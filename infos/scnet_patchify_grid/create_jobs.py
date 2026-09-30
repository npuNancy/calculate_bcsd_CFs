#!/usr/bin/env python3
"""Generate portable Slurm packs for all four production-V2 CF models; never submit."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import re
import shlex

HERE = Path(__file__).resolve().parent
MODELS = ("CANESM5", "MPI-ESM1-2-HR", "MRI-ESM2-0", "BCC-CSM2-MR")
SCENARIOS = ("ssp126", "ssp245", "ssp585")
TECHS = ("wind", "solar")


def catalog():
    with (HERE / "accounts.csv").open(encoding="utf-8-sig") as f:
        accounts = list(csv.DictReader(f))
    with (HERE / "作业分工/patch_assignment.csv").open(encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    workers = [r["username"] for r in accounts if r["role"] == "worker"]
    aggregators = [r["username"] for r in accounts if r["role"] == "aggregator"]
    if len(workers) != 14 or len(set(workers)) != 14 or aggregators != ["acjpoxgsdu"] or aggregators[0] in workers:
        raise ValueError("expected 14 unique workers and aggregator acjpoxgsdu")
    if len(rows) != 47 or len({r["patch_id"] for r in rows}) != 47:
        raise ValueError("expected 47 distinct patch assignments")
    for r in rows:
        if r["username"] not in workers or not re.fullmatch(r"R\d{2}C\d{2}", r["patch_id"]):
            raise ValueError("invalid patch assignment")
    return workers, {r["patch_id"]: r for r in rows}


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--jobs-dir", required=True)
    p.add_argument("--code-sha", required=True, help="Pinned full deployed Git SHA")
    p.add_argument("--patches", nargs="+", help="Optional local test subset; preparation omits this flag")
    p.add_argument("--partition", default="wzhctest")
    p.add_argument("--resource-profile", default="cf_grid_v2_v1")
    p.add_argument("--cpus-per-task", type=int, default=10)
    p.add_argument("--processes", type=int, default=8)
    p.add_argument("--time", default="24:00:00")
    p.add_argument("--tile-shape", nargs=2, type=int, default=[64, 64])
    p.add_argument("--time-chunk", type=int, default=240)
    p.add_argument("--compress-level", type=int, choices=range(10), default=2)
    p.add_argument("--merge-final", action="store_true", default=False)
    p.add_argument("--dry-run", action="store_true")
    return p


def unit_id(model, scenario, tech, patch):
    return f"cf-grid-v2/{model}/{scenario}/{tech}/{patch}"


def render(args, row, workers):
    command = ["python", "infos/scnet_patchify_grid/run_job.py", "--model", row["model"],
               "--scenario", row["scenario"], "--tech", row["tech"], "--patch", row["patch"],
               "--code-sha", args.code_sha, "--resource-profile", args.resource_profile,
               "--processes", str(args.processes), "--tile-shape", *map(str, args.tile_shape),
               "--time-chunk", str(args.time_chunk), "--compress-level", str(args.compress_level)]
    if args.merge_final:
        command.append("--merge-final")
    return "\n".join([
        "#!/usr/bin/env bash", f"#SBATCH --job-name={row['job_name']}",
        f"#SBATCH --partition={args.partition}", "#SBATCH --nodes=1", "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={args.cpus_per_task}", f"#SBATCH --time={args.time}",
        "#SBATCH --output=logs/%x-%j.out", "#SBATCH --error=logs/%x-%j.err",
        "set -eo pipefail", "module load apps/git/2.30.2", "git --version",
        ': "${CF_ENV_FILE:?set CF_ENV_FILE to the external account campaign.env}"',
        'source "$CF_ENV_FILE"',
        ': "${CF_CLIMATE_ACTIVATE:?set CF_CLIMATE_ACTIVATE}"',
        'source "$CF_CLIMATE_ACTIVATE" climate',
        "set -euo pipefail", "umask 0002",
        'run_user="$(id -un)"', 'case "$run_user" in',
        f"  {'|'.join(workers)}) ;;",
        '  *) echo "Not an authorized CF worker: $run_user" >&2; exit 2 ;;', "esac",
        ': "${SLURM_JOB_ID:?run on a Slurm compute node}"',
        ': "${CF_REPO:?set CF_REPO}"',
        'export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1',
        'cd "$CF_REPO"', shlex.join(command), "",
    ])


def build(args):
    workers, assignments = catalog()
    patches = args.patches if args.patches is not None else list(assignments)
    if not patches or len(set(patches)) != len(patches) or any(p not in assignments for p in patches):
        raise ValueError("invalid/duplicate patches")
    if not re.fullmatch(r"[0-9a-f]{40}", args.code_sha):
        raise ValueError("code-sha must be 40 lowercase hex characters")
    for name in ("partition", "resource_profile"):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", getattr(args, name)):
            raise ValueError(f"unsafe {name}")
    if not 1 <= args.processes <= min(8, args.cpus_per_task) or min(args.tile_shape) < 1 or args.time_chunk < 1:
        raise ValueError("require positive chunks and 1 <= processes <= min(8, cpus)")
    if not re.fullmatch(r"(?:\d+-)?\d{2,3}:[0-5]\d:[0-5]\d", args.time) or not any(int(v) for v in re.split("[-:]", args.time)):
        raise ValueError("invalid time limit")
    scripts, jobs = {}, []
    # Stable order is also the ready-queue priority; it is not a stage barrier.
    patches = sorted(patches, key=lambda p: (-int(assignments[p]["land_points_reference"]), p))
    for model, scenario, tech, patch in itertools.product(MODELS, SCENARIOS, TECHS, patches):
        name = f"cfg2_{model}_{scenario}_{tech}_{patch}"
        row = {"unit_id": unit_id(model, scenario, tech, patch), "model": model, "scenario": scenario,
               "tech": tech, "patch": patch, "years": "2015-2060", "job_name": name,
               "script": name + ".sh", "logical_owner": assignments[patch]["username"],
               "submit_username": None, "resource_profile": args.resource_profile,
               "cpus": args.cpus_per_task, "processes": args.processes,
               "merge_final": args.merge_final, "external_bc_variables": ["uas", "vas"] if tech == "wind" else ["rsds", "tas", "uas", "vas"],
               "expected_manifest": f"outputs/{model}/{scenario}/{patch}/{tech}/manifest.json"}
        script = render(args, row, workers)
        row["script_sha256"] = hashlib.sha256(script.encode()).hexdigest()
        scripts[row["script"]] = script; jobs.append(row)
    manifest = {"schema": "cf-grid-v2-job-pack-v1", "code_sha": args.code_sha,
                "resource_profile": args.resource_profile, "models": list(MODELS),
                "scenarios": list(SCENARIOS), "techs": list(TECHS), "patches": patches,
                "workers": workers, "aggregator": "acjpoxgsdu", "count": len(jobs),
                "years": "2015-2060", "merge_final": args.merge_final, "jobs": jobs}
    return manifest, scripts


def main(argv=None):
    p = parser(); args = p.parse_args(argv)
    try:
        manifest, scripts = build(args)
        target = Path(args.jobs_dir).expanduser().resolve()
        if not args.dry_run:
            if target.is_relative_to(HERE.parents[1]):
                raise ValueError("generated packs must be outside the checkout")
            for name in [*scripts, "manifest.json"]:
                path = target / name
                if path.exists() or path.is_symlink():
                    raise FileExistsError(f"immutable pack already exists: {path}; choose a new directory")
            target.mkdir(parents=True, exist_ok=True)
            for name, script in scripts.items():
                with (target / name).open("x", encoding="utf-8") as f:
                    f.write(script)
                (target / name).chmod(0o750)
            with (target / "manifest.json").open("x", encoding="utf-8") as f:
                json.dump(manifest, f, ensure_ascii=False, indent=2); f.write("\n")
    except (OSError, ValueError) as exc:
        p.exit(2, f"error: {exc}\n")
    print(json.dumps({"dry_run": args.dry_run, "jobs": manifest["count"], "models": list(MODELS)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
