#!/usr/bin/env python3
"""Generate identical complete station-CF Slurm packs; never submit jobs."""
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
MODELS = ('CANESM5', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'BCC-CSM2-MR')
SCENARIOS = ('ssp126', 'ssp245', 'ssp585')
TECHS = ('wind', 'solar')
AGGREGATOR = 'acjpoxgsdu'


def catalog():
    with (HERE/'accounts.csv').open(encoding='utf-8-sig', newline='') as f:
        accounts = list(csv.DictReader(f))
    with (HERE/'作业分工/patch_assignment.csv').open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    workers = [a['username'] for a in accounts if a['role'] == 'worker']
    if (len(accounts) != 15 or len(workers) != 14 or len(set(workers)) != 14
            or [a['username'] for a in accounts if a['role'] == 'aggregator'] != [AGGREGATOR]
            or AGGREGATOR in workers):
        raise ValueError('expected fourteen workers and aggregator acjpoxgsdu')
    by_user = {a['username']: a for a in accounts}
    if len(rows) != 47 or len({r['patch_id'] for r in rows}) != 47:
        raise ValueError('expected exactly 47 unique patches')
    for r in rows:
        if (r['username'] not in workers or not re.fullmatch(r'R\d{2}C\d{2}', r['patch_id'])
                or int(r['land_points_reference']) < 0 or r['host'] != by_user[r['username']]['host']):
            raise ValueError('invalid patch assignment')
    return workers, {r['patch_id']: r for r in rows}


def unit_id(model, climate, station, tech, patch):
    return f'station-cf-v1/{model}/climate-{climate}/station-{station}/{tech}/{patch}'


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--jobs-dir', required=True)
    p.add_argument('--code-sha', required=True)
    p.add_argument('--config-sha256', required=True, help='SHA256 of frozen inputs/campaign_config.json')
    p.add_argument('--resource-profile', default='station_cf_v1')
    p.add_argument('--partition', default='wzhctest')
    p.add_argument('--cpus-per-task', type=int, default=10)
    p.add_argument('--prepare-cpus', type=int, default=10)
    p.add_argument('--publish-cpus', type=int, default=10)
    p.add_argument('--processes', type=int, default=8)
    p.add_argument('--time', default='24:00:00')
    p.add_argument('--time-chunk', type=int, default=240)
    p.add_argument('--station-chunk', type=int, default=1024)
    p.add_argument('--compress-level', type=int, choices=range(10), default=2)
    p.add_argument('--dry-run', action='store_true')
    return p


def render(args, row, workers):
    command = ['python', 'infos/scnet_patchify_stations/run_job.py', '--stage', row['stage'],
               '--code-sha', args.code_sha, '--config-sha256', args.config_sha256,
               '--resource-profile', args.resource_profile]
    if row['stage'] == 'extract':
        for key in ('model', 'climate_scenario', 'station_scenario', 'tech', 'patch'):
            command.extend(['--'+key.replace('_', '-'), row[key]])
        for key in ('processes', 'time_chunk', 'station_chunk', 'compress_level'):
            command.extend(['--'+key.replace('_', '-'), str(getattr(args, key))])
    return '\n'.join([
        '#!/usr/bin/env bash', f'#SBATCH --job-name={row["job_name"]}',
        f'#SBATCH --partition={args.partition}', '#SBATCH --nodes=1', '#SBATCH --ntasks=1',
        f'#SBATCH --cpus-per-task={row["cpus"]}', f'#SBATCH --time={args.time}',
        '#SBATCH --output=logs/%x-%j.out', '#SBATCH --error=logs/%x-%j.err',
        'set -eo pipefail', 'module load apps/git/2.30.2', 'git --version',
        ': "${SCF_ENV_FILE:?set SCF_ENV_FILE to the external campaign.env}"',
        'source "$SCF_ENV_FILE"', ': "${SCF_CLIMATE_ACTIVATE:?set SCF_CLIMATE_ACTIVATE}"',
        'source "$SCF_CLIMATE_ACTIVATE" climate', 'set -euo pipefail', 'umask 0002',
        'run_user="$(id -un)"', 'case "$run_user" in',
        f'  {"|".join(workers)}) ;;',
        '  *) echo "Not an authorized station-CF worker: $run_user" >&2; exit 2 ;;', 'esac',
        ': "${SLURM_JOB_ID:?requires a Slurm compute job}"',
        ': "${SCF_REPO:?set SCF_REPO}"',
        'export SCF_JOB_SCRIPT="$(realpath "${BASH_SOURCE[0]}")"',
        'export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1',
        'cd "$SCF_REPO"', shlex.join(command), '',
    ])


def build(args):
    workers, assignments = catalog()
    if not re.fullmatch('[0-9a-f]{40}', args.code_sha) or not re.fullmatch('[0-9a-f]{64}', args.config_sha256):
        raise ValueError('require full lowercase code SHA and configuration SHA256')
    for key in ('resource_profile', 'partition'):
        if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*', getattr(args, key)):
            raise ValueError('unsafe '+key)
    if (not 1 <= args.processes <= min(8, args.cpus_per_task) or min(args.prepare_cpus, args.publish_cpus) < 1
            or min(args.time_chunk, args.station_chunk) < 1):
        raise ValueError('invalid CPU/process/chunk configuration')
    if not re.fullmatch(r'(?:\d+-)?\d{2,3}:[0-5]\d:[0-5]\d', args.time) or not any(int(x) for x in re.split('[-:]', args.time)):
        raise ValueError('invalid walltime')
    patches = sorted(assignments, key=lambda p: (-int(assignments[p]['land_points_reference']), p))
    rows = [{'stage': 'prepare', 'unit_id': 'station-cf-v1/prepare', 'job_name': 'stcf_prepare',
             'logical_owner': workers[0], 'cpus': args.prepare_cpus, 'processes': 1,
             'expected_manifest': 'prepared.json', 'depends_on': []}]
    for m, c, s, t, patch in itertools.product(MODELS, SCENARIOS, SCENARIOS, TECHS, patches):
        rows.append({'stage': 'extract', 'unit_id': unit_id(m, c, s, t, patch),
                     'model': m, 'climate_scenario': c, 'station_scenario': s, 'tech': t, 'patch': patch,
                     'job_name': f'stcf_{m}_c{c[3:]}_s{s[3:]}_{t}_{patch}',
                     'logical_owner': assignments[patch]['username'], 'cpus': args.cpus_per_task,
                     'processes': args.processes, 'time_chunk': args.time_chunk,
                     'station_chunk': args.station_chunk, 'compress_level': args.compress_level,
                     'years': '2015-2060', 'depends_on': ['station-cf-v1/prepare'],
                     'expected_manifest': f'outputs/{m}/climate_{c}/station_{s}/{patch}/{t}/manifest.json'})
    rows.append({'stage': 'publish', 'unit_id': 'station-cf-v1/publish', 'job_name': 'stcf_publish',
                 'logical_owner': workers[0], 'cpus': args.publish_cpus, 'processes': 1,
                 'expected_manifest': 'index/authoritative_index.json', 'depends_on': 'all_extract_succeeded'})
    scripts = {}
    for row in rows:
        row.update(resource_profile=args.resource_profile, submit_username=None, script=row['job_name']+'.sh')
        script = render(args, row, workers)
        if row['script'] in scripts:
            raise ValueError('script filename collision')
        scripts[row['script']] = script
        row['script_sha256'] = hashlib.sha256(script.encode()).hexdigest()
    manifest = {'schema': 'station-cf-slurm-v1', 'code_sha': args.code_sha, 'config_sha256': args.config_sha256,
                'resource_profile': args.resource_profile, 'models': list(MODELS),
                'climate_scenarios': list(SCENARIOS), 'station_scenarios': list(SCENARIOS),
                'techs': list(TECHS), 'patches': patches, 'workers': workers, 'aggregator': AGGREGATOR,
                'years': '2015-2060', 'extract_count': 3384, 'count': len(rows),
                'stage_counts': {'prepare': 1, 'extract': 3384, 'publish': 1}, 'jobs': rows}
    return manifest, scripts


def main(argv=None):
    p = parser(); args = p.parse_args(argv)
    try:
        manifest, scripts = build(args)
        target = Path(args.jobs_dir).expanduser().resolve()
        if not args.dry_run:
            if target.is_relative_to(HERE.parents[1]):
                raise ValueError('generated packs must be outside checkout')
            content = {**scripts, 'manifest.json': json.dumps(manifest, ensure_ascii=False, indent=2)+'\n'}
            for name in content:
                if (target/name).exists() or (target/name).is_symlink():
                    raise FileExistsError('immutable pack exists; use another directory')
            target.mkdir(parents=True, exist_ok=True)
            for name, value in content.items():
                with (target/name).open('x', encoding='utf-8') as f:
                    f.write(value)
                (target/name).chmod(0o750 if name.endswith('.sh') else 0o640)
    except (ValueError, OSError) as exc:
        p.exit(2, f'error: {exc}\n')
    print(json.dumps({'dry_run': args.dry_run, 'stages': manifest['stage_counts'], 'scripts': manifest['count']}))


if __name__ == '__main__':
    main()
