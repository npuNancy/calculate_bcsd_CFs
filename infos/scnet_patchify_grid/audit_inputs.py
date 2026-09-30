#!/usr/bin/env python3
"""Classify all CF units from manifests and file metadata, without opening NC arrays."""
import argparse
import csv
from datetime import datetime
import io
import json
from pathlib import Path
import sys
from zoneinfo import ZoneInfo

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from grid_cf_sources import VARS, missing, relocated_path, select_sources
from create_jobs import MODELS, SCENARIOS, TECHS, catalog, unit_id


def audit(root, patches):
    checked = datetime.now(ZoneInfo("Asia/Shanghai")).isoformat(timespec="seconds")
    rows = []
    for model in MODELS:
        for scenario in SCENARIOS:
            for tech in TECHS:
                for patch in sorted(patches):
                    mode, sources, absent = select_sources(root, model, scenario, tech, patch)
                    final_absent, invalid = [], []
                    for variable, source in sources.items():
                        if len(source["blocks"]) != 8:
                            raise ValueError("production_v2 inventory requires eight year blocks")
                        path = relocated_path(source["task"]["output"], root, ("outputs", model, scenario, variable))
                        side = Path(str(path)+".json")
                        for p in (path, side):
                            if missing(p):
                                final_absent.append(str(p))
                            elif mode == "final" and p.stat().st_size == 0:
                                invalid.append(str(p))
                        if mode == "final" and not missing(side) and side.stat().st_size:
                            meta = json.loads(side.read_text())
                            if (meta.get("kind") != "global-final-patch" or meta.get("variable") != variable
                                    or meta.get("patch_id") != patch
                                    or relocated_path(meta["output"], root, ("outputs", model, scenario, variable)) != path):
                                invalid.append(str(side))
                        if mode == "blocks":
                            invalid.extend(str(p) for p in source["blocks"]+source["sidecars"] if p.stat().st_size == 0)
                    rows.append({"unit_id": unit_id(model, scenario, tech, patch), "model": model,
                                 "scenario": scenario, "tech": tech, "patch": patch, "input_mode": mode,
                                 "required_block_nc": len(VARS[tech])*8,
                                 "missing_block_nc": sum(p.suffix == ".nc" for p in absent),
                                 "missing_block_sidecars": sum(p.suffix == ".json" for p in absent),
                                 "missing_final_files": len(final_absent), "invalid_files": len(invalid),
                                 "files_ready": not invalid and (mode == "blocks" or not final_absent),
                                 "checked_at": checked, "bcsd_root": str(root)})
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bcsd-root", required=True)
    p.add_argument("--output", help="CSV destination; stdout when omitted")
    args = p.parse_args()
    _, patches = catalog()
    rows = audit(Path(args.bcsd_root).resolve(), patches)
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
    if args.output:
        target = Path(args.output)
        with target.open("x", encoding="utf-8", newline="") as f:
            f.write(stream.getvalue())
    else:
        print(stream.getvalue(), end="")


if __name__ == "__main__":
    main()
