"""Freeze grid-CF inputs, snapshot catalogs and global nearest mappings."""
from __future__ import annotations

import argparse
import itertools
import json
import time
from concurrent.futures import ProcessPoolExecutor
from contextlib import nullcontext
from multiprocessing import get_context
from pathlib import Path

from grid_cf_io import atomic_json, digest, output_lock, read_json
from station_cf_catalog import SCENARIOS, TECHS, build_catalog, scenario, sha
from station_cf_mapping import build_mapping, grid_fingerprints
import station_cf_io as io

MODELS = ('CANESM5', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'BCC-CSM2-MR')


def source_key(model, climate, tech, patch):
    return '/'.join((model, climate, tech, patch))


def task_id(model, climate, station, tech, patch):
    return f'station-cf-v1/{model}/climate-{climate}/station-{station}/{tech}/{patch}'


def validate_config(config, patches, *, sample=False):
    config = dict(config)
    for key in ('climate_scenarios', 'station_scenarios'):
        values = [scenario(x, station=key == 'station_scenarios') for x in config[key]]
        if len(set(values)) != len(values) or not values:
            raise ValueError(f'duplicate/empty {key}')
        config[key] = values
    for key in ('models', 'techs'):
        if not config[key] or len(set(config[key])) != len(config[key]):
            raise ValueError(f'duplicate/empty {key}')
    if not set(config['techs']) <= set(TECHS):
        raise ValueError('unsupported tech')
    for name in list(config['models']) + list(patches):
        if not name or '/' in name or '\\' in name or name in ('.', '..'):
            raise ValueError('unsafe model/patch identifier')
    if config.get('spatial_method', 'nearest_regular') != 'nearest_regular':
        raise ValueError('only nearest_regular is supported')
    if not sample and (set(config['models']) != set(MODELS) or len(patches) != 47
                       or set(config['climate_scenarios']) != set(SCENARIOS)
                       or set(config['station_scenarios']) != set(SCENARIOS)
                       or set(config['techs']) != set(TECHS) or config['years'] != '2015-2060'):
        raise ValueError('production scope must be 4 GCM × 3 Climate × 3 Station × 2 tech × 47 patches, 2015-2060')
    for key in ('station_files', 'station_reference_files'):
        normalized = {}
        for k, v in config[key].items():
            canonical = scenario(k, station=True)
            if canonical in normalized:
                raise ValueError('duplicate Station file alias')
            normalized[canonical] = v
        config[key] = normalized
        if not set(config['station_scenarios']) <= set(normalized):
            raise ValueError(f'missing {key}')
    return config


def enumerate_tasks(config, patches):
    return [{'unit_id': task_id(m, c, s, t, p), 'model': m, 'climate_scenario': c,
             'station_scenario': s, 'tech': t, 'source_patch': p,
             'source_key': source_key(m, c, t, p), 'years': config['years']}
            for m, c, s, t, p in itertools.product(config['models'], config['climate_scenarios'],
                                                  config['station_scenarios'], config['techs'], sorted(patches))]


_TIME_CACHE = {}


def _scan_request(request):
    path, kwargs = request
    return io.scan_source(path, **kwargs, time_cache=_TIME_CACHE)


def _progress(stage, completed, total, started):
    elapsed = time.monotonic() - started
    print(json.dumps({'stage': stage, 'completed': completed, 'total': total,
                      'elapsed_seconds': round(elapsed, 2),
                      'remaining_estimate_seconds': round(elapsed * (total-completed) / completed, 2)
                      if completed else None}), flush=True)


def scan_sources(config, grid_root, patches, release, *, sample, processes):
    sources, groups, time_contracts = {}, {}, {}
    total = len(config['models'])*len(config['climate_scenarios'])*len(config['techs'])*len(patches)
    started = time.monotonic()
    executor = (ProcessPoolExecutor(max_workers=processes, mp_context=get_context('spawn'))
                if processes > 1 else nullcontext(None))
    with executor as pool:
        for m, c, t in itertools.product(config['models'], config['climate_scenarios'], config['techs']):
            grids = {}
            ordered = sorted(patches)
            requests = [(grid_root/'outputs'/m/c/p/t/'manifest.json',
                         dict(model=m, climate_scenario=c, tech=t, patch=p, years=config['years'],
                              expected_blocks=None if sample else 8)) for p in ordered]
            results = pool.map(_scan_request, requests, chunksize=1) if pool else map(_scan_request, requests)
            for p, (source, grid) in zip(ordered, results):
                if not sample:
                    provenance = read_json(source['manifest'])['provenance']
                    if provenance['implementation']['code_sha'] != release['code_sha']:
                        raise ValueError('source code SHA differs from frozen release')
                sources[source_key(m, c, t, p)] = source; grids[p] = grid
                times = [b['time'] for b in source['blocks']]
                key = '/'.join((m, c, t))
                if key in time_contracts and time_contracts[key] != times:
                    raise ValueError('cross-patch time contract mismatch')
                time_contracts[key] = times
                if len(sources) == 1 or len(sources) % len(patches) == 0:
                    _progress('source_index', len(sources), total, started)
            fingerprint = digest(grid_fingerprints(grids))
            groups.setdefault(fingerprint, grids)
            for p in patches:
                sources[source_key(m, c, t, p)]['spatial_group'] = fingerprint
    return sources, groups


def prepare(config_path, output_root, *, sample=False, processes=1):
    if not 1 <= processes <= 16:
        raise ValueError("prepare processes must be between 1 and 16")
    config_path = Path(config_path).resolve(); root = Path(output_root).resolve()
    raw = read_json(config_path); grid_root = io.resolve(raw['grid_cf_root'], config_path.parent)
    if root == grid_root or root.is_relative_to(grid_root) or grid_root.is_relative_to(root):
        raise ValueError('source and station output roots must be separate')
    patch_path = grid_root/'inputs/patch_manifest.json'
    patches = read_json(patch_path)['patches']; config = validate_config(raw, patches, sample=sample)
    root.mkdir(parents=True, exist_ok=True)
    with output_lock(root/'.prepare.lock'):
        prepared_path = root/'prepared.json'
        if prepared_path.exists():
            old = io.load_prepared(prepared_path)
            if old['config_hash'] != digest(config) or old['sample'] != sample or old['implementation'] != io.implementation():
                raise ValueError('prepared run differs; use a new output root')
            for catalog in old['catalogs'].values():
                for path, expected in catalog['files'].items():
                    if sha(path) != expected:
                        raise ValueError('catalog changed')
                for key in ('source', 'reference'):
                    io.check_file(catalog[key])
                    if sha(catalog[key]['path']) != catalog[key+'_sha256']:
                        raise ValueError('station source/reference changed')
            for mapping in old['mappings'].values():
                if sha(mapping['path']) != mapping['sha256']:
                    raise ValueError('mapping changed')
            for source in old['sources'].values():
                io.check_source(source)
            return prepared_path
        completion = read_json(grid_root/'runtime/completion.json')
        release_path = grid_root/'inputs/release.json'
        release = read_json(release_path)
        if completion['status'] != 'COMPLETED' or completion['run_id'] != release['run_id']:
            raise ValueError('grid CF run is not completed')
        if sha(release_path) != completion['release_sha256']:
            raise ValueError('grid release hash mismatch')
        if not sample and (completion['succeeded_units'] != 1128 or completion['year_nc_files'] != 9024):
            raise ValueError('grid CF production scope incomplete')
        if not sample and sha(patch_path) != release['patch_manifest']['sha256']:
            raise ValueError('source patch manifest changed from release')
        assets = [io.hashed_file(config_path), io.hashed_file(patch_path), io.hashed_file(release_path),
                  io.hashed_file(grid_root/'runtime/completion.json')]
        sources, groups = scan_sources(config, grid_root, patches, release, sample=sample, processes=processes)
        catalogs = {}; stage_started = time.monotonic()
        station_root = io.resolve(config['stations_root'], config_path.parent)
        for s in config['station_scenarios']:
            catalogs[s] = build_catalog(io.resolve(config['station_files'][s], station_root), s, root,
                                        reference_csv=io.resolve(config['station_reference_files'][s], config_path.parent))
            _progress('catalogs', len(catalogs), len(config['station_scenarios']), stage_started)
        mappings = {}; stage_started = time.monotonic()
        for fingerprint, grids in groups.items():
            for s, t in itertools.product(config['station_scenarios'], config['techs']):
                mappings[f'{fingerprint}/{s}/{t}'] = build_mapping(
                    catalogs[s], t, grids, patches, root, config.get('max_distance_deg', 0.15))
                _progress('mappings', len(mappings), len(groups)*len(config['station_scenarios'])*len(config['techs']), stage_started)
        tasks = enumerate_tasks(config, patches)
        for task in tasks:
            group = sources[task['source_key']]['spatial_group']
            task['mapping_key'] = f'{group}/{task["station_scenario"]}/{task["tech"]}'
            task['station_count'] = mappings[task['mapping_key']]['patch_counts'][task['source_patch']]
            task['unique_grid_count'] = mappings[task['mapping_key']]['unique_grid_counts'][task['source_patch']]
        atomic_json(root/'inputs/station_cf_config.json', config)
        atomic_json(root/'inputs/grid_cf_index.json', sources)
        input_release = {'source_run_id': completion['run_id'], 'release': release, 'assets': assets,
                         'capacity_semantics': 'snapshot_total', 'implementation': io.implementation()}
        atomic_json(root/'inputs/input_release.json', input_release)
        atomic_json(root/'index/units.json', {'tasks': tasks})
        for path in ('inputs/station_cf_config.json', 'inputs/grid_cf_index.json', 'inputs/input_release.json', 'index/units.json'):
            assets.append(io.hashed_file(root/path))
        result = {'schema': 'station-cf-v1', 'output_root': str(root), 'config_hash': digest(config),
                  'sample': sample, 'source_run_id': completion['run_id'], 'implementation': io.implementation(),
                  'assets': assets, 'sources': sources, 'catalogs': catalogs, 'mappings': mappings,
                  'tasks': tasks, 'capacity_semantics': 'snapshot_total'}
        result['storage_estimate'] = {
            'raw_cf_bytes': sum(4*t['station_count']*sum(b['time']['count'] for b in sources[t['source_key']]['blocks']) for t in tasks),
            'year_nc_files': sum(len(sources[t['source_key']]['blocks']) for t in tasks if t['station_count']),
            'note': 'float32 logical bytes; compression and temporary files are not predicted'}
        result['identity'] = digest(result)
        for source in sources.values():
            io.check_source(source)
        atomic_json(prepared_path, result)
        return prepared_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True); parser.add_argument('--output-root', required=True)
    parser.add_argument('--sample', action='store_true', help='explicit small-fixture scope; never a full production release')
    parser.add_argument("--processes", type=int, default=16, help="parallel source-index workers (1-16)")
    args = parser.parse_args()
    print(prepare(args.config, args.output_root, sample=args.sample, processes=args.processes))


if __name__ == '__main__':
    main()
