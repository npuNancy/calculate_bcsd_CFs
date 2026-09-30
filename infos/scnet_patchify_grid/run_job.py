#!/usr/bin/env python3
"""Compute-node CF adapter: verify campaign context, run one unit, publish receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import pwd
import re
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
from create_jobs import MODELS, SCENARIOS, TECHS, catalog, unit_id


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require_env(name):
    value = os.environ.get(name)
    if not value:
        raise ValueError(f"missing environment variable {name}")
    return value


def stat_identity(path):
    path = Path(path).resolve(strict=True); s = path.stat()
    return {"path": str(path), "size": s.st_size, "mtime_ns": s.st_mtime_ns}


def inside(path, root):
    path = Path(path).resolve(strict=True)
    if not path.is_relative_to(Path(root).resolve(strict=True)):
        raise ValueError(f"path outside shared run root: {path}")
    return path


def preflight(args):
    workers, patches = catalog()
    user = pwd.getpwuid(os.getuid()).pw_name
    if user not in workers or args.patch not in patches:
        raise ValueError("unauthorized worker or patch")
    job_id = require_env("SLURM_JOB_ID")
    run_id = require_env("CF_RUN_ID")
    if not re.fullmatch(r"\d+", job_id) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", run_id):
        raise ValueError("invalid Slurm job ID/run ID")
    if os.environ.get("SLURM_ARRAY_TASK_ID") or int(require_env("SLURM_CPUS_PER_TASK")) < args.processes:
        raise ValueError("requires an individual job with enough CPUs")
    repo = Path(require_env("CF_REPO")).resolve(strict=True)
    if repo != ROOT:
        raise ValueError("CF_REPO differs from this checkout")
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()
    branch = subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip()
    if head != args.code_sha or dirty or branch != "develop-patch-grid":
        raise ValueError("requires clean develop-patch-grid at the pinned SHA")
    shared = Path(require_env("CF_SHARED_ROOT")).resolve(strict=True)
    home = (Path("/work/home") / "acjpoxgsdu").resolve(strict=True)
    if not shared.is_relative_to(home) or shared == home:
        raise ValueError("shared results must be inside aggregator home")
    release_path = inside(require_env("CF_RELEASE_FILE"), shared)
    release_hash = sha256(release_path)
    if release_hash != require_env("CF_RELEASE_SHA256"):
        raise ValueError("campaign release changed")
    release = json.loads(release_path.read_text())
    if (release["release"] != "production_v2" or release["run_id"] != run_id
            or release["code_sha"] != head or release["models"] != list(MODELS)
            or release["years"] != "2015-2060" or release["bcsd_spot_check_accepted"] is not True
            or Path(release["shared_root"]).resolve(strict=True) != shared):
        raise ValueError("release identity/scope mismatch")
    bcsd = Path(require_env("CF_BCSD_ROOT")).resolve(strict=True)
    if str(bcsd) != release["bcsd_root"]:
        raise ValueError("BCSD root differs from frozen release")
    land = Path(require_env("CF_LAND_PLAN")).resolve(strict=True)
    patch = Path(require_env("CF_PATCH_MANIFEST")).resolve(strict=True)
    if stat_identity(land) != release["land_plan"]:
        raise ValueError("land plan identity differs from release")
    if str(patch) != release["patch_manifest"]["path"] or sha256(patch) != release["patch_manifest"]["sha256"]:
        raise ValueError("patch manifest differs from release")
    if set(json.loads(patch.read_text())["patches"]) != set(patches):
        raise ValueError("patch manifest must contain the campaign's 47 active patches")
    units = release["input_units"]
    if set(units) != {"uas", "vas", "tas", "rsds"}:
        raise ValueError("release must declare all four CF input units")
    # Only frozen metadata is checked here. The compute entry validates this unit's blocks.
    return {"user": user, "job_id": job_id, "run_id": run_id, "shared": shared,
            "release_path": release_path, "release_hash": release_hash, "bcsd": bcsd,
            "land": land, "patch": patch, "units": units}


def compute_arguments(args, runtime):
    from patchify_grid_cf import parser as compute_parser
    variables = ("uas", "vas") if args.tech == "wind" else ("rsds", "tas", "uas", "vas")
    values = ["--bcsd-root", str(runtime["bcsd"]), "--model", args.model, "--scenario", args.scenario,
              "--tech", args.tech, "--patch", args.patch, "--patch-manifest", str(runtime["patch"]),
              "--land-plan", str(runtime["land"]), "--years", "2015-2060",
              "--output-root", str(runtime["shared"] / "outputs"), "--processes", str(args.processes),
              "--tile-shape", *map(str, args.tile_shape), "--time-chunk", str(args.time_chunk),
              "--compress-level", str(args.compress_level), "--input-units",
              *[f"{v}={runtime['units'][v]}" for v in variables]]
    if args.merge_final:
        values.append("--merge-final")
    return compute_parser().parse_args(values)


def receipt_evidence(manifest, args, shared):
    """Check only JSON/stat; scientific validation already ran on the compute node."""
    manifest = inside(manifest, shared)
    m = json.loads(manifest.read_text())
    p = m["provenance"]
    if m["status"] != "COMPLETED" or m["merge_final"] != args.merge_final or len(m["blocks"]) != 8:
        raise ValueError("CF unit incomplete")
    for key in ("model", "scenario", "tech", "patch"):
        if p[key] != getattr(args, key):
            raise ValueError(f"CF manifest {key} mismatch")
    if p["implementation"]["code_sha"] != args.code_sha:
        raise ValueError("CF code identity mismatch")
    artifacts = []
    for block in m["blocks"]:
        path = inside(block["path"], shared)
        side_path = inside(str(path) + ".json", shared)
        side = json.loads(side_path.read_text())
        if block["status"] != "COMPLETED" or side["status"] != "COMPLETED" or side["file"] != stat_identity(path):
            raise ValueError("CF block incomplete/changed")
        artifacts.append({"file": side["file"], "sidecar": stat_identity(side_path)})
    if args.merge_final:
        path = inside(m["final"]["path"], shared)
        side_path = inside(str(path) + ".json", shared)
        side = json.loads(side_path.read_text())
        if m["final"]["status"] != "COMPLETED" or side["status"] != "COMPLETED" or side["file"] != stat_identity(path):
            raise ValueError("CF final incomplete/changed")
        artifacts.append({"file": side["file"], "sidecar": stat_identity(side_path)})
    return {"manifest": stat_identity(manifest), "identity": m["identity"], "artifacts": artifacts}


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model", choices=MODELS, required=True)
    p.add_argument("--scenario", choices=SCENARIOS, required=True)
    p.add_argument("--tech", choices=TECHS, required=True)
    p.add_argument("--patch", required=True)
    p.add_argument("--code-sha", required=True)
    p.add_argument("--resource-profile", required=True)
    p.add_argument("--processes", type=int, default=8)
    p.add_argument("--tile-shape", type=int, nargs=2, default=[64, 64])
    p.add_argument("--time-chunk", type=int, default=240)
    p.add_argument("--compress-level", type=int, default=2)
    p.add_argument("--merge-final", action="store_true", default=False)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    runtime = preflight(args)
    from patchify_grid_cf import compute
    from grid_cf_io import atomic_json
    start = time.monotonic()
    manifest = compute(compute_arguments(args, runtime))
    evidence = receipt_evidence(manifest, args, runtime["shared"])
    if sha256(runtime["release_path"]) != runtime["release_hash"]:
        raise ValueError("release changed during computation")
    uid = unit_id(args.model, args.scenario, args.tech, args.patch)
    receipt = runtime["shared"] / "runtime/receipts" / uid / f"{runtime['job_id']}.json"
    atomic_json(receipt, {"status": "COMPLETED", "unit_id": uid, "run_id": runtime["run_id"],
                         "submit_username": runtime["user"], "slurm_job_id": runtime["job_id"],
                         "code_sha": args.code_sha, "release_sha256": runtime["release_hash"],
                         "resource_profile": args.resource_profile, "elapsed_seconds": time.monotonic()-start,
                         "merge_final": args.merge_final, **evidence})
    print(json.dumps({"unit_id": uid, "receipt": str(receipt)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
