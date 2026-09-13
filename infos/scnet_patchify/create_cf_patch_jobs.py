#!/usr/bin/env python3
"""Generate one Slurm script per model/scenario/patch/technology CF unit."""
from __future__ import annotations
import argparse, itertools, re, shlex
from pathlib import Path

MODELS=("CANESM5","MPI-ESM1-2-HR","MRI-ESM2-0","BCC-CSM2-MR")
SCENARIOS=("ssp126","ssp245","ssp585")
TECHS=("wind","solar")
PATCHES=()

def _token(x):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", x): raise ValueError(f"unsafe token {x!r}")
    return x

def build_parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--models", nargs="+", default=list(MODELS)); p.add_argument("--scenarios", nargs="+", default=list(SCENARIOS)); p.add_argument("--patches", nargs="+", required=True, metavar="PATCH"); p.add_argument("--techs", nargs="+", choices=TECHS, default=list(TECHS)); p.add_argument("--years", default="2015-2060"); p.add_argument("--bcsd-root", required=True); p.add_argument("--patch-manifest", required=True); p.add_argument("--stations-csv", required=True); p.add_argument("--output-root", required=True); p.add_argument("--project-dir", default=str(Path(__file__).resolve().parents[2])); p.add_argument("--jobs-dir", required=True); p.add_argument("--logs-dir", required=True); p.add_argument("--partition", default="wzhctest"); p.add_argument("--account", default=None); p.add_argument("--cpus-per-task", type=int, default=10); p.add_argument("--time", default=None); p.add_argument("--overwrite", action="store_true"); p.add_argument("--dry-run", action="store_true"); return p

def render(a, model, scenario, patch, tech):
    jid=f"cfp_{model}_{scenario}_{patch}_{tech}"; q=shlex.quote
    lines=["#!/bin/bash", "#SBATCH --job-name="+jid, f"#SBATCH --partition={a.partition}", f"#SBATCH --cpus-per-task={a.cpus_per_task}", f"#SBATCH --output={a.logs_dir}/{jid}_%j.out", f"#SBATCH --error={a.logs_dir}/{jid}_%j.out"]
    if a.account: lines.append(f"#SBATCH --account={a.account}")
    if a.time: lines.append(f"#SBATCH --time={a.time}")
    cmd=["python", str(Path(a.project_dir)/"patchify_station_cf.py"), "--bcsd-root", a.bcsd_root, "--model", model, "--scenario", scenario, "--patch", patch, "--patch-manifest", a.patch_manifest, "--stations-csv", a.stations_csv, "--tech", tech, "--years", a.years, "--output-root", a.output_root]
    if a.overwrite: cmd.append("--overwrite")
    # The SCNet activation scripts reference optional variables (for example
    # MAGPLUS_HOME).  Activate before enabling nounset so a valid job does not
    # fail during environment setup.
    lines += ["source /work/home/acbpgywfpz/miniconda3/bin/activate climate", "set -euo pipefail", f"mkdir -p {q(a.logs_dir)} {q(a.output_root)}", f"cd {q(a.project_dir)}", " ".join(q(x) for x in cmd)]
    return jid, "\n".join(lines)+"\n"

def main(argv=None):
    a=build_parser().parse_args(argv); out=Path(a.jobs_dir).expanduser();
    if not a.dry_run: out.mkdir(parents=True,exist_ok=True); Path(a.logs_dir).expanduser().mkdir(parents=True,exist_ok=True)
    rows=[]
    for model,scenario,patch,tech in itertools.product(a.models,a.scenarios,a.patches,a.techs):
        for x in (model,scenario,patch): _token(x)
        jid,body=render(a,model,scenario,patch,tech)
        if a.dry_run:
            rows.append({"unit_id":jid,"model":model,"scenario":scenario,"patch":patch,"tech":tech}); continue
        path=out/(jid+".sh");
        if path.exists() and not a.overwrite: raise FileExistsError(path)
        path.write_text(body,encoding="utf-8"); path.chmod(0o750); rows.append({"unit_id":jid,"model":model,"scenario":scenario,"patch":patch,"tech":tech,"script":str(path)})
    if not a.dry_run:
        import json; (out/"manifest.json").write_text(json.dumps(rows,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(f"generated {len(rows)} jobs in {out} (dry-run={a.dry_run})")
if __name__ == "__main__": main()
