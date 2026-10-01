"""Extract one Climate × Station unit, then optionally publish an audited index."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import os
from pathlib import Path
import uuid

from grid_cf_io import atomic_json, digest, output_lock, read_json
from station_cf_catalog import scenario, sha, write_csv
import station_cf_io as io
from station_cf_mapping import read_mapping, STATUS


def task_context(prep, task, *, verify_catalog=True):
    s = task['station_scenario']; catalog = prep['catalogs'][s]
    if verify_catalog:
        for path, expected in catalog['files'].items():
            if sha(path) != expected:
                raise ValueError('catalog changed')
    return {**task, 'source': prep['sources'][task['source_key']],
            'mapping': prep['mappings'][task['mapping_key']], 'catalog_identity': catalog['identity'],
            'source_run_id': prep['source_run_id'], 'code_sha': prep['implementation']['code_sha'],
            'implementation_sha256': digest(prep['implementation'])}


def unit_directory(prep, task):
    return (Path(prep['output_root'])/'outputs'/task['model']/('climate_'+task['climate_scenario'])/
            ('station_'+task['station_scenario'])/task['source_patch']/task['tech'])


def unit_identity(task, profile, implementation):
    return digest({'task': task, 'profile': profile, 'implementation': implementation})


def extract_unit(prepared, *, model, climate_scenario, station_scenario, tech, patch,
                 processes=8, time_chunk=240, station_chunk=1024, compress_level=2):
    if processes <= 0 or time_chunk <= 0 or station_chunk <= 0 or not 0 <= compress_level <= 9:
        raise ValueError('invalid process/chunk/compression configuration')
    climate_scenario = scenario(climate_scenario); station_scenario = scenario(station_scenario, station=True)
    prep = io.load_prepared(prepared)
    impl = io.implementation()
    if prep['implementation'] != impl:
        raise ValueError('implementation changed after preparation')
    matches = [t for t in prep['tasks'] if (t['model'], t['climate_scenario'], t['station_scenario'], t['tech'], t['source_patch']) ==
               (model, climate_scenario, station_scenario, tech, patch)]
    if len(matches) != 1:
        raise ValueError('unit not uniquely defined in prepared scope')
    task = task_context(prep, matches[0]); io.check_source(task['source'])
    directory = unit_directory(prep, task); directory.mkdir(parents=True, exist_ok=True)
    profile = {'time_chunk': time_chunk, 'station_chunk': station_chunk, 'compress_level': compress_level}
    identity = unit_identity(task, profile, impl)
    manifest = directory/'manifest.json'
    with output_lock(directory/'.unit.lock'):
        old = read_json(manifest) if manifest.exists() else None
        if old and old['identity'] != identity:
            raise ValueError('unit identity/profile changed; use a new output root')
        frame = io.sort_mapping(io.selected_mapping(task), task['source']['blocks'][0]['path'], tech)
        if len(frame) != task['station_count']:
            raise ValueError('prepared station count mismatch')
        state = {'unit_id': task['unit_id'], 'identity': identity, 'prepared_identity': prep['identity'],
                 'status': 'RUNNING', 'profile': profile, 'implementation': impl, 'station_count': len(frame),
                 'context': {k: task[k] for k in ('model', 'climate_scenario', 'station_scenario', 'tech', 'source_patch')},
                 'blocks': []}
        if frame.empty:
            io.check_source(task['source'])
            state['status'] = 'EMPTY_NO_STATIONS'
            atomic_json(directory/'audit.json', {'status': state['status'], 'identity': identity, 'station_count': 0})
            atomic_json(manifest, state)
            return manifest
        payloads = []
        for b in task['source']['blocks']:
            prior = next((x for x in (old or {}).get('blocks', []) if x['index'] == b['index']), None)
            path = Path(prior['path']) if prior else directory/'blocks'/identity/f'cf_{b["start_year"]}-{b["end_year"]}.nc'
            bid = io.block_id(identity, b)
            reusable = io.completed(path, task, b, frame, bid)
            if not reusable and (path.exists() or Path(str(path)+'.json').exists()):
                path = directory/'blocks'/identity/('attempt_'+uuid.uuid4().hex)/path.name
            state['blocks'].append({'index': b['index'], 'path': str(path), 'identity': bid,
                                    'status': 'COMPLETED' if reusable else 'pending'})
            if not reusable:
                payloads.append((task, b, profile, identity, str(path)))
        atomic_json(manifest, state)
        try:
            count = min(processes, len(payloads), 8)
            if count <= 1:
                for payload in payloads:
                    row = io.extract_block(payload); state['blocks'][row['index']] = row; atomic_json(manifest, state)
            else:
                for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
                    os.environ[name] = '1'
                with ProcessPoolExecutor(max_workers=count, mp_context=multiprocessing.get_context('spawn')) as pool:
                    futures = [pool.submit(io.extract_block, p) for p in payloads]
                    try:
                        for future in as_completed(futures):
                            row = future.result(); state['blocks'][row['index']] = row; atomic_json(manifest, state)
                    except BaseException:
                        for future in futures:
                            future.cancel()
                        raise
            for b, row in zip(task['source']['blocks'], state['blocks']):
                if not io.completed(row['path'], task, b, frame, io.block_id(identity, b)):
                    raise ValueError('year block verification failed')
                io.validate_station(row['path'], task, b, frame, row['identity'], sample=True)
            io.check_source(task['source'])
            audit = {'status': 'COMPLETED', 'identity': identity, 'station_count': len(frame),
                     'year_files': len(state['blocks']), 'blocks': [io.hashed_file(str(r['path'])+'.json') for r in state['blocks']],
                     'independent_value_audit': 'fixed-seed source samples per block',
                     'source_files': task['source']['files']}
            atomic_json(directory/'audit.json', audit)
            state['status'] = 'COMPLETED'; atomic_json(manifest, state)
        except BaseException as exc:
            state.update(status='FAILED', error=str(exc)); atomic_json(manifest, state)
            raise
    return manifest


def publish(prepared):
    """Serial complete-scope audit and index publication; no task execution."""
    import pandas as pd
    prep = io.load_prepared(prepared); root = Path(prep['output_root'])
    with output_lock(root/'.publish.lock'):
        entries, coverage, unmapped = [], [], []
        checked_catalogs = set(); cached_mapping_key = None; cached_mapping_frame = None
        for base in prep['tasks']:
            task = task_context(prep, base, verify_catalog=base['station_scenario'] not in checked_catalogs)
            checked_catalogs.add(base['station_scenario'])
            io.check_source(task['source'])
            if cached_mapping_key != base['mapping_key']:
                if sha(task['mapping']['path']) != task['mapping']['sha256']:
                    raise ValueError('mapping changed')
                cached_mapping_frame = read_mapping(task['mapping']['path'])
                cached_mapping_key = base['mapping_key']
            directory = unit_directory(prep, task); path = directory/'manifest.json'
            if not path.exists():
                raise ValueError(f'incomplete unit: {base["unit_id"]}')
            state = read_json(path)
            expected = unit_identity(task, state['profile'], prep['implementation'])
            if state['identity'] != expected or state['prepared_identity'] != prep['identity']:
                raise ValueError('unit manifest identity mismatch')
            frame = io.sort_mapping(io.selected_mapping(task, cached_mapping_frame), task['source']['blocks'][0]['path'], task['tech'])
            if len(frame) != base['station_count'] or state['station_count'] != len(frame):
                raise ValueError('station count mismatch')
            if len(frame):
                if state['status'] != 'COMPLETED' or len(state['blocks']) != len(task['source']['blocks']):
                    raise ValueError(f'incomplete unit: {base["unit_id"]}')
                for b, row in zip(task['source']['blocks'], state['blocks']):
                    if row['index'] != b['index'] or not io.completed(row['path'], task, b, frame, io.block_id(expected, b)):
                        raise ValueError('incomplete year output')
                    io.validate_station(row['path'], task, b, frame, io.block_id(expected, b), sample=True)
            elif state['status'] != 'EMPTY_NO_STATIONS' or state['blocks']:
                raise ValueError('invalid empty unit')
            audit = read_json(directory/'audit.json')
            if audit['status'] != state['status'] or audit['identity'] != expected:
                raise ValueError('missing unit audit')
            valid_count = sum(read_json(str(r['path'])+'.json')['valid_count'] for r in state['blocks'])
            missing_count = sum(read_json(str(r['path'])+'.json')['missing_count'] for r in state['blocks'])
            entries.append({**base, 'valid_count': valid_count, 'missing_count': missing_count, 'manifest': io.hashed_file(path), 'status': state['status'],
                            'blocks': [{**r, 'file': io.file_identity(r['path']),
                                        'sidecar': io.hashed_file(str(r['path'])+'.json')} for r in state['blocks']]})
        seen = set()
        for t in prep['tasks']:
            key = (t['model'], t['climate_scenario'], t['station_scenario'], t['tech'])
            if key in seen:
                continue
            seen.add(key); mapping = prep['mappings'][t['mapping_key']]
            if sha(mapping['path']) != mapping['sha256']:
                raise ValueError('mapping changed')
            frame = read_mapping(mapping['path']); counts = {k: int((frame.mapping_status == v).sum()) for k, v in STATUS.items()}
            if sum(counts.values()) != len(frame) or len(frame) != prep['catalogs'][t['station_scenario']]['catalogs'][t['tech']]['count']:
                raise ValueError('catalog coverage conservation failed')
            assigned = sum(e['station_count'] for e in entries if tuple(e[k] for k in ('model','climate_scenario','station_scenario','tech')) == key)
            if assigned != counts['MATCHED'] + counts['OUTSIDE_DOMAIN']:
                raise ValueError('patch coverage conservation failed')
            context = dict(zip(('model','climate_scenario','station_scenario','tech'), key))
            selected_entries = [e for e in entries if tuple(e[k] for k in ('model','climate_scenario','station_scenario','tech')) == key]
            nt = sum(b['time']['count'] for b in prep['sources'][t['source_key']]['blocks'])
            valid_count = sum(e['valid_count'] for e in selected_entries)
            missing_count = sum(e['missing_count'] for e in selected_entries) + (counts['TOO_FAR']+counts['NO_SOURCE_PATCH'])*nt
            if valid_count + missing_count != len(frame)*nt:
                raise ValueError('global time sample conservation failed')
            coverage.append({**context, 'catalog_count': len(frame), **counts, 'valid_count': valid_count,
                             'missing_count': missing_count, 'valid_fraction': valid_count/(len(frame)*nt) if len(frame) else None})
            unmapped.append(frame.loc[frame.mapping_status.isin([2, 3])].assign(**context))
        write_csv(pd.DataFrame(coverage), root/'index/coverage_summary.csv')
        write_csv(pd.concat(unmapped, ignore_index=True), root/'index/unmapped_stations.csv.gz')
        summary = {'status': 'COMPLETED', 'scope': 'sample' if prep['sample'] else 'production',
                   'units': len(entries), 'nonempty_units': sum(e['status'] == 'COMPLETED' for e in entries),
                   'empty_units': sum(e['status'] == 'EMPTY_NO_STATIONS' for e in entries),
                   'year_nc_files': sum(len(e['blocks']) for e in entries),
                   'year_nc_bytes': sum(b['file']['size'] for e in entries for b in e['blocks'])}
        atomic_json(root/'index/validation_summary.json', summary)
        index = {'schema': 'station-cf-v1', 'status': 'COMPLETED', 'prepared': io.hashed_file(prepared),
                 'prepared_identity': prep['identity'], 'entries': entries, 'summary': summary}
        index['identity'] = digest(index)
        atomic_json(root/'index/authoritative_index.json', index)
        return root/'index/authoritative_index.json'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', required=True); parser.add_argument('--publish', action='store_true')
    for k in ('model', 'climate-scenario', 'station-scenario', 'tech', 'patch'):
        parser.add_argument('--'+k)
    parser.add_argument('--processes', type=int, default=8); parser.add_argument('--time-chunk', type=int, default=240)
    parser.add_argument('--station-chunk', type=int, default=1024); parser.add_argument('--compress-level', type=int, default=2)
    args = vars(parser.parse_args()); publication = args.pop('publish')
    if publication:
        if any(args[k] for k in ('model','climate_scenario','station_scenario','tech','patch')):
            parser.error('--publish cannot be combined with unit selectors')
        print(publish(args['prepared']))
    else:
        if any(not args[k] for k in ('model','climate_scenario','station_scenario','tech','patch')):
            parser.error('both scenarios, model, tech and patch are required')
        print(extract_unit(**args))


if __name__ == '__main__':
    main()
