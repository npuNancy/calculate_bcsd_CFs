"""Metadata-only selection of one BCSD input layout for a CF unit."""
import json
from pathlib import Path

VARS = {"wind": ("uas", "vas"), "solar": ("rsds", "tas", "uas", "vas")}


def relocated_path(value, root, suffix):
    parts = Path(value).parts
    if tuple(parts[-len(suffix)-1:-1]) != tuple(suffix):
        raise ValueError(f"unexpected BCSD path identity: {value}")
    path = root.joinpath(*suffix, parts[-1]).resolve()
    if not path.is_relative_to((root / suffix[0]).resolve()):
        raise ValueError(f"input outside production {suffix[0]}: {path}")
    return path


def missing(path):
    # Permission and I/O failures must not masquerade as missing inputs.
    try:
        path.stat()
    except FileNotFoundError:
        return True
    return False


def select_sources(root, model, scenario, tech, patch, mode="auto"):
    if mode not in ("auto", "blocks", "final"):
        raise ValueError(f"invalid input mode: {mode}")
    root = Path(root).expanduser().resolve()
    records, absent = {}, []
    for variable in VARS[tech]:
        manifest = root / "manifests" / f"{model}__{scenario}__{variable}__{patch}.json"
        meta = json.loads(manifest.read_text())
        if meta.get("kind") != "global-bcsd-three-stage":
            raise ValueError(f"invalid BCSD manifest: {manifest}")
        for key, value in (("model", model), ("scenario", scenario), ("variable", variable)):
            if key in meta and meta[key] != value:
                raise ValueError(f"manifest {key} mismatch: {manifest}")
        tasks = [t for t in meta["final_tasks"] if t["patch_id"] == patch and t["variable"] == variable]
        if len(tasks) != 1 or len(tasks[0]["block_files"]) != len(meta["time_block_contract"]["blocks"]):
            raise ValueError("manifest final task/block count mismatch")
        task = tasks[0]
        paths = [relocated_path(v, root, ("blocks", model, scenario, variable, patch)) for v in task["block_files"]]
        sides = [p.with_name("block_" + p.stem + ".json") for p in paths]
        absent.extend(p for p in [*paths, *sides] if missing(p))
        records[variable] = {"manifest": manifest, "meta": meta, "task": task,
                             "blocks": paths, "sidecars": sides}
    selected = ("final" if absent else "blocks") if mode == "auto" else mode
    if selected == "blocks" and absent:
        raise FileNotFoundError(f"missing BCSD block input: {absent[0]}")
    if selected == "final":
        for variable, record in records.items():
            if "output" not in record["task"]:
                raise FileNotFoundError(f"manifest has no final input for {variable}")
            record["final"] = relocated_path(record["task"]["output"], root, ("outputs", model, scenario, variable))
    return selected, records, absent
