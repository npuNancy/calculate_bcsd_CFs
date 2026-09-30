"""Portable full-campaign packs and the compute-node receipt boundary."""
import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

HERE = Path(__file__).resolve().parents[1] / "infos/scnet_patchify_grid"


def load(name):
    spec = importlib.util.spec_from_file_location("cf_test_" + name, HERE / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


jobs = load("create_jobs")
runner = load("run_job")
SHA = "a" * 40


def arguments(target, *extra):
    return jobs.parser().parse_args(["--jobs-dir", str(target), "--code-sha", SHA, *extra])


def test_full_campaign_packs_are_portable(tmp_path):
    workers, patches = jobs.catalog()
    expected, scripts = jobs.build(arguments(tmp_path / "first"))
    assert len(patches) == 47 and len(workers) == 14
    assert expected["aggregator"] not in workers
    assert expected["count"] == len(scripts) == 1128
    assert len({r["unit_id"] for r in expected["jobs"]}) == 1128
    counts = Counter((r["model"], r["scenario"]) for r in expected["jobs"])
    assert len(counts) == 12 and set(counts.values()) == {94}
    for row in expected["jobs"]:
        assert row["logical_owner"] in workers
        assert row["processes"] == 8 and row["merge_final"] is False
        assert row["script_sha256"] == hashlib.sha256(scripts[row["script"]].encode()).hexdigest()
    # Materialize every worker's entire inventory; account-specific paths must not leak in.
    for worker in workers:
        target = tmp_path / worker / "pack"
        jobs.main(["--jobs-dir", str(target), "--code-sha", SHA])
        assert json.loads((target / "manifest.json").read_text()) == expected
        actual = {p.name: p.read_text() for p in target.glob("*.sh")}
        assert actual == scripts
    for script in scripts.values():
        assert "/work/home/" not in script
        assert "--merge-final" not in script
        assert "#SBATCH --account" not in script
        subprocess.run(["bash", "-n"], input=script, text=True, check=True, capture_output=True)


def test_dry_run_and_immutable_pack(tmp_path):
    target = tmp_path / "pack"
    base = ["--jobs-dir", str(target), "--code-sha", SHA]
    jobs.main([*base, "--dry-run"])
    assert not target.exists()
    jobs.main(base)
    before = (target / "manifest.json").read_bytes()
    with pytest.raises(SystemExit):
        jobs.main(base)
    assert (target / "manifest.json").read_bytes() == before
    with pytest.raises(SystemExit):
        jobs.main(["--jobs-dir", str(HERE / "generated"), "--code-sha", SHA])


@pytest.mark.parametrize("extra", [
    ["--code-sha", "HEAD"], ["--processes", "9"], ["--cpus-per-task", "4"],
    ["--partition", "x\n#SBATCH --account=other"], ["--time", "00:00:00"],
    ["--tile-shape", "0", "64"], ["--patches", "R02C09", "R02C09"],
])
def test_invalid_pack_rejected(tmp_path, extra):
    with pytest.raises(ValueError):
        jobs.build(arguments(tmp_path, *extra))


@pytest.mark.parametrize("tech", ["wind", "solar"])
@pytest.mark.parametrize("merge", [False, True])
def test_compute_contract(tmp_path, tech, merge):
    args = runner.parser().parse_args([
        "--model", "CANESM5", "--scenario", "ssp126", "--tech", tech,
        "--patch", "R02C09", "--code-sha", SHA, "--resource-profile", "test",
        *(["--merge-final"] if merge else []),
    ])
    runtime = {"bcsd": tmp_path / "bcsd", "shared": tmp_path / "shared",
               "land": tmp_path / "land.nc", "patch": tmp_path / "patch.json",
               "units": {"uas": "m/s", "vas": "m/s", "tas": "K", "rsds": "W/m2"}}
    computed = runner.compute_arguments(args, runtime)
    assert computed.processes == 8 and computed.merge_final is merge
    assert computed.years == "2015-2060"
    assert Path(computed.output_root) == runtime["shared"] / "outputs"
    assert computed.overwrite is False and computed.parts_root is None
    assert len(computed.input_units) == (2 if tech == "wind" else 4)
    _, scripts = jobs.build(arguments(tmp_path, *(["--merge-final"] if merge else [])))
    assert all(("--merge-final" in script) is merge for script in scripts.values())


def receipt_fixture(root, merge=False):
    args = SimpleNamespace(model="CANESM5", scenario="ssp126", tech="wind", patch="R02C09",
                           code_sha=SHA, merge_final=merge)
    state = {"identity": "unit-identity", "status": "COMPLETED", "merge_final": merge,
             "provenance": {k: getattr(args, k) for k in ("model", "scenario", "tech", "patch")},
             "blocks": []}
    state["provenance"]["implementation"] = {"code_sha": SHA}
    for i in range(8 + int(merge)):
        path = root / f"cf_{i}.nc"
        path.write_bytes(b"receipt stat fixture")
        side = {"status": "COMPLETED", "file": runner.stat_identity(path)}
        Path(str(path) + ".json").write_text(json.dumps(side))
        item = {"status": "COMPLETED", "path": str(path)}
        if i < 8:
            state["blocks"].append(item)
        else:
            state["final"] = item
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps(state))
    return args, manifest, state


@pytest.mark.parametrize("merge", [False, True])
def test_receipt_validates_products(tmp_path, merge):
    args, manifest, _ = receipt_fixture(tmp_path, merge)
    evidence = runner.receipt_evidence(manifest, args, tmp_path)
    assert len(evidence["artifacts"]) == 8 + int(merge)
    (tmp_path / "cf_0.nc").write_bytes(b"changed")
    with pytest.raises(ValueError, match="changed"):
        runner.receipt_evidence(manifest, args, tmp_path)


@pytest.mark.parametrize("defect", ["missing", "incomplete", "wrong_model", "wrong_sha", "outside", "missing_final"])
def test_receipt_rejects_bad_evidence(tmp_path, defect):
    root = tmp_path / "shared"
    root.mkdir()
    args, manifest, state = receipt_fixture(root, defect == "missing_final")
    if defect == "missing":
        (root / "cf_0.nc.json").unlink()
    elif defect == "incomplete":
        state["blocks"].pop()
    elif defect == "wrong_model":
        state["provenance"]["model"] = "other"
    elif defect == "wrong_sha":
        state["provenance"]["implementation"]["code_sha"] = "b" * 40
    elif defect == "outside":
        outside = tmp_path / "outside.nc"
        outside.write_bytes(b"foreign")
        state["blocks"][0]["path"] = str(outside)
    else:
        state["final"]["status"] = "pending"
    manifest.write_text(json.dumps(state))
    with pytest.raises((ValueError, FileNotFoundError)):
        runner.receipt_evidence(manifest, args, root)


def test_aggregator_cannot_run(monkeypatch):
    monkeypatch.setattr(runner.pwd, "getpwuid", lambda _: SimpleNamespace(pw_name="acjpoxgsdu"))
    with pytest.raises(ValueError, match="unauthorized"):
        runner.preflight(SimpleNamespace(patch="R02C09"))
