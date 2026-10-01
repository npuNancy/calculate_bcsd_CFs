#!/usr/bin/env python3
"""Compute-node adapter for preparation, station CF extraction and publication."""
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
from infos.scnet_patchify_stations.create_jobs import MODELS, SCENARIOS, TECHS, AGGREGATOR, catalog, unit_id


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def env(name):
    value = os.environ.get(name)
    if not value:
        raise ValueError('missing environment variable '+name)
    return value


def stat(path):
    p = Path(path).resolve(strict=True); s = p.stat()
    return {'path': str(p), 'size': s.st_size, 'mtime_ns': s.st_mtime_ns}


def inside(path, root):
    p = Path(path).resolve(strict=True)
    if not p.is_relative_to(Path(root).resolve(strict=True)):
        raise ValueError('artifact outside shared root')
    return p


def job_unit(args):
    if args.stage == 'extract':
        if any(getattr(args, k) is None for k in ('model', 'climate_scenario', 'station_scenario', 'tech', 'patch')):
            raise ValueError('extract requires model, both scenarios, tech and patch')
        return unit_id(args.model, args.climate_scenario, args.station_scenario, args.tech, args.patch)
    if any(getattr(args, k) is not None for k in ('model', 'climate_scenario', 'station_scenario', 'tech', 'patch')):
        raise ValueError('global stages do not take extract selectors')
    return 'station-cf-v1/'+args.stage


def prepared_contract(path, shared, code_sha):
    """JSON contract only; science entry points perform metadata/array validation."""
    prep = read(path); _, patches = catalog()
    if (prep['sample'] is not False or prep['identity'] != digest({k:v for k,v in prep.items() if k != 'identity'})
            or Path(prep['output_root']).resolve() != Path(shared).resolve()
            or prep['implementation']['code_sha'] != code_sha or prep['capacity_semantics'] != 'snapshot_total'):
        raise ValueError('prepared production identity mismatch')
    expected = {unit_id(m,c,s,t,p) for m in MODELS for c in SCENARIOS for s in SCENARIOS for t in TECHS for p in patches}
    if len(prep['tasks']) != 3384 or {t['unit_id'] for t in prep['tasks']} != expected:
        raise ValueError('prepared must contain all 3384 unique units')
    for t in prep['tasks']:
        if (t['unit_id'] != unit_id(t['model'],t['climate_scenario'],t['station_scenario'],t['tech'],t['source_patch'])
                or t['years'] != '2015-2060'):
            raise ValueError('prepared unit fields mismatch')
    return prep


def preflight(args):
    workers, patches = catalog(); user = pwd.getpwuid(os.getuid()).pw_name
    uid = job_unit(args)
    if user not in workers or (args.stage == 'extract' and args.patch not in patches):
        raise ValueError('unauthorized worker or patch')
    job_id, run_id = env('SLURM_JOB_ID'), env('SCF_RUN_ID')
    if not re.fullmatch(r'\d+', job_id) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', run_id):
        raise ValueError('invalid job or run identity')
    if os.environ.get('SLURM_ARRAY_TASK_ID') or os.environ.get('SLURM_JOB_ACCOUNT', user) != user:
        raise ValueError('requires individual job billed to actual worker')
    if Path(env('SCF_REPO')).resolve(strict=True) != ROOT:
        raise ValueError('SCF_REPO differs from checkout')
    head = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    dirty = subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    branch = subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip()
    if head != args.code_sha or dirty or branch != 'develop-patch-grid':
        raise ValueError('requires clean develop-patch-grid at pinned SHA')
    shared = Path(env('SCF_SHARED_ROOT')).resolve(strict=True)
    home = (Path('/work/home')/AGGREGATOR).resolve(strict=True)
    if shared == home or not shared.is_relative_to(home):
        raise ValueError('shared root must be physically inside aggregator home')
    release_path = inside(env('SCF_RELEASE_FILE'), shared); release_hash = sha(release_path)
    if release_hash != env('SCF_RELEASE_SHA256'):
        raise ValueError('release changed')
    release = read(release_path)
    if (release['schema'] != 'station-cf-campaign-v1' or release['run_id'] != run_id
            or release['code_sha'] != head or Path(release['shared_root']).resolve() != shared
            or release['models'] != list(MODELS) or release['climate_scenarios'] != list(SCENARIOS)
            or release['station_scenarios'] != list(SCENARIOS) or release['years'] != '2015-2060'
            or release['capacity_semantics'] != 'snapshot_total' or release['bcsd_spot_check_accepted'] is not True):
        raise ValueError('campaign release scope mismatch')
    config_path = inside(env('SCF_CONFIG'), shared)
    if (config_path == shared/'inputs/station_cf_config.json' or str(config_path) != release['config']['path']
            or sha(config_path) != release['config']['sha256'] or sha(config_path) != args.config_sha256):
        raise ValueError('configuration changed or collides with generated config')
    pack_path = Path(env('SCF_PACK_MANIFEST')).resolve(strict=True); pack = read(pack_path)
    if (pack['code_sha'] != head or pack['config_sha256'] != args.config_sha256
            or pack['resource_profile'] != args.resource_profile or pack['count'] != 3386):
        raise ValueError('wrong job pack')
    rows = [r for r in pack['jobs'] if r['unit_id'] == uid]
    if len(rows) != 1:
        raise ValueError('unit absent or duplicated in pack')
    row = rows[0]
    if sha(env('SCF_JOB_SCRIPT')) != row['script_sha256'] or row['stage'] != args.stage:
        raise ValueError('executed script hash/stage differs from pack')
    if int(env('SLURM_CPUS_PER_TASK')) < row['cpus']:
        raise ValueError('Slurm CPU allocation below prepared resource profile')
    if args.stage == 'extract':
        for key in ('model','climate_scenario','station_scenario','tech','patch','processes','time_chunk','station_chunk','compress_level'):
            if getattr(args,key) != row[key]:
                raise ValueError('extract arguments differ from pack')
        for key in ('time_chunk','station_chunk','compress_level'):
            if getattr(args,key) != release['encoding'][key]:
                raise ValueError('encoding changed; use a new scientific run')
    prepared = shared/'prepared.json'; prepared_hash = None
    if args.stage != 'prepare':
        prepared_hash = sha(prepared)
        if prepared_hash != env('SCF_PREPARED_SHA256'):
            raise ValueError('prepared file changed')
        prep = prepared_contract(prepared, shared, head)
        if prep['config_hash'] != digest(read(config_path)):
            raise ValueError('prepared configuration differs from frozen campaign config')
    return {'user': user, 'job_id': job_id, 'run_id': run_id, 'unit_id': uid, 'shared': shared,
            'release_path': release_path, 'release_hash': release_hash, 'config': config_path,
            'prepared': prepared, 'prepared_hash': prepared_hash, 'pack_hash': sha(pack_path), 'row': row}


def extract_evidence(path, args, runtime):
    shared = runtime['shared']; path = inside(path, shared); m = read(path)
    uid = job_unit(args); prepared = read(runtime['prepared'])
    tasks = [t for t in prepared['tasks'] if t['unit_id'] == uid]
    if len(tasks) != 1:
        raise ValueError('unknown prepared unit')
    task = tasks[0]; source = prepared['sources'][task['source_key']]
    context = dict(model=args.model,climate_scenario=args.climate_scenario,station_scenario=args.station_scenario,
                   tech=args.tech,source_patch=args.patch)
    if (m['unit_id'] != uid or m['prepared_identity'] != prepared['identity'] or m['context'] != context
            or m['station_count'] != task['station_count'] or m['implementation']['code_sha'] != args.code_sha
            or m['profile'] != {k:getattr(args,k) for k in ('time_chunk','station_chunk','compress_level')}):
        raise ValueError('unit manifest context/profile mismatch')
    audit_path = inside(path.parent/'audit.json', shared); audit = read(audit_path)
    if audit['identity'] != m['identity'] or audit['status'] != m['status']:
        raise ValueError('unit audit mismatch')
    artifacts = []
    if m['status'] == 'EMPTY_NO_STATIONS':
        if m['station_count'] != 0 or m['blocks']:
            raise ValueError('invalid empty-unit marker')
    elif m['status'] == 'COMPLETED':
        if m['station_count'] <= 0 or len(m['blocks']) != 8 or len(source['blocks']) != 8:
            raise ValueError('unit year blocks incomplete')
        for i, (b, src) in enumerate(zip(m['blocks'], source['blocks'])):
            nc_path = inside(b['path'], shared); side_path = inside(str(nc_path)+'.json', shared); side = read(side_path)
            expected = digest({'unit':m['identity'],'source':src['identity'],'time':src['time'],'index':i})
            if (b['index'] != i or b['status'] != 'COMPLETED' or side['status'] != 'COMPLETED'
                    or b['identity'] != expected or side['identity'] != expected or side['file'] != stat(nc_path)
                    or side['source_block_identity'] != src['identity'] or side['station_count'] != m['station_count']
                    or side['time_count'] != src['time']['count']):
                raise ValueError('year output identity/stat mismatch')
            artifacts.append({'file':stat(nc_path),'sidecar':{**stat(side_path),'sha256':sha(side_path)}})
    else:
        raise ValueError('unit incomplete')
    return {'scientific_status':m['status'],'identity':m['identity'],'station_count':m['station_count'],
            'manifest':{**stat(path),'sha256':sha(path)}, 'audit':{**stat(audit_path),'sha256':sha(audit_path)},
            'artifacts':artifacts}


def execute(args, runtime):
    if args.stage == 'prepare':
        from prepare_station_cf import prepare
        path = prepare(runtime['config'], runtime['shared'])
        prep = prepared_contract(path, runtime['shared'], args.code_sha)
        return {'scientific_status':'COMPLETED','prepared_identity':prep['identity'],
                'prepared':{**stat(path),'sha256':sha(path)},'units':len(prep['tasks']),
                'storage_estimate':prep['storage_estimate']}
    if args.stage == 'extract':
        from extract_station_cf import extract_unit
        values = {k:getattr(args,k) for k in ('model','climate_scenario','station_scenario','tech','patch',
                                             'processes','time_chunk','station_chunk','compress_level')}
        path = extract_unit(runtime['prepared'], **values)
        return extract_evidence(path,args,runtime)
    from extract_station_cf import publish
    path = publish(runtime['prepared']); index = read(path); summary = index['summary']
    if (index['status'] != 'COMPLETED' or summary['scope'] != 'production' or summary['units'] != 3384
            or summary['nonempty_units']+summary['empty_units'] != 3384
            or summary['year_nc_files'] != 8*summary['nonempty_units']):
        raise ValueError('global publication incomplete')
    files = [path, *(runtime['shared']/'index'/name for name in
                    ('validation_summary.json','coverage_summary.csv','unmapped_stations.csv.gz'))]
    return {'scientific_status':'COMPLETED','index_identity':index['identity'],'summary':summary,
            'artifacts':[{**stat(inside(p,runtime['shared'])),'sha256':sha(p)} for p in files]}


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=('prepare','extract','publish'),required=True)
    p.add_argument('--code-sha',required=True); p.add_argument('--config-sha256',required=True)
    p.add_argument('--resource-profile',required=True)
    p.add_argument('--model',choices=MODELS); p.add_argument('--climate-scenario',choices=SCENARIOS)
    p.add_argument('--station-scenario',choices=SCENARIOS); p.add_argument('--tech',choices=TECHS)
    p.add_argument('--patch'); p.add_argument('--processes',type=int,default=8)
    p.add_argument('--time-chunk',type=int,default=240); p.add_argument('--station-chunk',type=int,default=1024)
    p.add_argument('--compress-level',type=int,default=2)
    return p


def main(argv=None):
    args = parser().parse_args(argv); runtime = preflight(args); start = time.monotonic()
    receipt = runtime['shared']/'runtime/receipts'/runtime['unit_id']/(runtime['job_id']+'.json')
    from grid_cf_io import atomic_json, output_lock
    with output_lock(str(receipt)+'.lock'):
        if receipt.exists():
            raise ValueError('JobID receipt already exists; do not overwrite')
        evidence = execute(args,runtime)
        if sha(runtime['release_path']) != runtime['release_hash'] or sha(runtime['config']) != args.config_sha256:
            raise ValueError('release/config changed during execution')
        if runtime['prepared_hash'] is not None and sha(runtime['prepared']) != runtime['prepared_hash']:
            raise ValueError('prepared changed during execution')
        atomic_json(receipt, {'status':'COMPLETED','stage':args.stage,'unit_id':runtime['unit_id'],
                             'run_id':runtime['run_id'],'submit_username':runtime['user'],'slurm_job_id':runtime['job_id'],
                             'code_sha':args.code_sha,'config_sha256':args.config_sha256,
                             'release_sha256':runtime['release_hash'],'prepared_sha256':sha(runtime['prepared']),
                             'resource_profile':args.resource_profile,'pack_sha256':runtime['pack_hash'],
                             'script_sha256':runtime['row']['script_sha256'],'elapsed_seconds':time.monotonic()-start,
                             **evidence})
    print(json.dumps({'unit_id':runtime['unit_id'],'receipt':str(receipt)}))


if __name__ == '__main__':
    main()
