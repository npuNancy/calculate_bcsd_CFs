"""Publish an index from previously verified extraction receipts."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from grid_cf_io import atomic_json, digest, output_lock, read_json
from station_cf_catalog import sha, write_csv
from station_cf_mapping import read_mapping
import station_cf_io as io


def verified_entry(payload):
    prep, task, row = payload
    if row['status'] not in ('succeeded', 'succeeded_empty') or not row.get('verified_at'):
        raise ValueError('extraction has not been verified')
    receipt_path = Path(row['receipt_path'])
    if sha(receipt_path) != row['receipt_sha256']:
        raise ValueError('verified receipt changed')
    receipt = read_json(receipt_path)
    for key in ('run_id', 'unit_id', 'stage', 'submit_username', 'slurm_job_id', 'code_sha',
                'config_sha256', 'release_sha256', 'prepared_sha256', 'resource_profile',
                'pack_sha256', 'script_sha256'):
        if receipt[key] != row[key]:
            raise ValueError(f'receipt identity mismatch: {key}')
    if receipt['status'] != 'COMPLETED' or receipt['station_count'] != task['station_count']:
        raise ValueError('receipt completion mismatch')
    for jid in (row['slurm_job_id'], row['slurm_job_id']+'.batch'):
        if not any(a['job_id']==jid and a['state']=='COMPLETED' and a['exit']=='0:0'
                   for a in row['accounting']):
            raise ValueError('missing successful Slurm evidence')
    empty = task['station_count']==0
    if receipt['scientific_status'] != ('EMPTY_NO_STATIONS' if empty else 'COMPLETED'):
        raise ValueError('invalid empty-unit receipt')
    artifacts = receipt['artifacts']; source = prep['sources'][task['source_key']]
    if len(artifacts) != (0 if empty else len(source['blocks'])):
        raise ValueError('receipt year count mismatch')
    root = Path(prep['output_root']); blocks = []
    if Path(receipt['manifest']['path']) != root/row['expected_manifest']:
        raise ValueError('receipt manifest path mismatch')
    for i, artifact in enumerate(artifacts):
        path = Path(artifact['file']['path'])
        if not path.is_relative_to(root) or artifact['sidecar']['path'] != str(path)+'.json':
            raise ValueError('receipt output path mismatch')
        source_block = source['blocks'][i]
        identity = digest(dict(unit=receipt['identity'], source=source_block['identity'],
                               time=source_block['time'], index=i))
        blocks.append(dict(index=i, path=str(path), identity=identity, status='COMPLETED', **artifact))
    return {**task, 'manifest':receipt['manifest'], 'status':receipt['scientific_status'],
            'blocks':blocks}


def publish(prepared, ledger_path, *, processes=8):
    """Read receipt metadata; never open or stat extraction output files."""
    import pandas as pd
    if processes < 1:
        raise ValueError('publish processes must be positive')
    prep = read_json(prepared); root = Path(prep['output_root'])
    if prep['identity'] != digest({k:v for k,v in prep.items() if k!='identity'}):
        raise ValueError('prepared identity mismatch')
    prepared_hash = sha(prepared); rows = read_json(ledger_path)['units']
    tasks = {t['unit_id']:t for t in prep['tasks']}
    extracts = {u:r for u,r in rows.items() if r['stage']=='extract'}
    if set(tasks) != set(extracts):
        raise ValueError('verified ledger scope mismatch')
    for row in extracts.values():
        if row['prepared_sha256'] != prepared_hash or row['code_sha'] != prep['implementation']['code_sha']:
            raise ValueError('verified ledger scientific identity mismatch')
    with output_lock(root/'.publish.lock'):
        payloads = [(prep, task, extracts[uid]) for uid,task in tasks.items()]
        with ThreadPoolExecutor(max_workers=processes) as pool:
            entries = list(pool.map(verified_entry, payloads))
        print(f'Publish loaded {len(entries)} verified receipts', flush=True)
        coverage, unmapped, frames = [], [], {}
        seen = set()
        for task in tasks.values():
            key = tuple(task[k] for k in ('model','climate_scenario','station_scenario','tech'))
            if key in seen:
                continue
            seen.add(key); mapping = prep['mappings'][task['mapping_key']]
            counts = mapping['counts']; count = mapping['catalog_count']
            if sum(counts.values()) != count or count != prep['catalogs'][task['station_scenario']]['catalogs'][task['tech']]['count']:
                raise ValueError('catalog coverage conservation failed')
            assigned = sum(e['station_count'] for e in entries if tuple(e[k] for k in ('model','climate_scenario','station_scenario','tech'))==key)
            if assigned != counts['MATCHED']+counts['OUTSIDE_DOMAIN']:
                raise ValueError('patch coverage conservation failed')
            context = dict(zip(('model','climate_scenario','station_scenario','tech'),key))
            coverage.append({**context, 'catalog_count':count, **counts})
            if task['mapping_key'] not in frames:
                if sha(mapping['path']) != mapping['sha256']:
                    raise ValueError('frozen mapping changed')
                frame = read_mapping(mapping['path'])
                frames[task['mapping_key']] = frame.loc[frame.mapping_status.isin([2,3])]
            unmapped.append(frames[task['mapping_key']].assign(**context))
        write_csv(pd.DataFrame(coverage),root/'index/coverage_summary.csv')
        write_csv(pd.concat(unmapped,ignore_index=True),root/'index/unmapped_stations.csv.gz')
        summary = dict(status='COMPLETED',scope='sample' if prep['sample'] else 'production',
                       units=len(entries),nonempty_units=sum(e['status']=='COMPLETED' for e in entries),
                       empty_units=sum(e['status']=='EMPTY_NO_STATIONS' for e in entries),
                       year_nc_files=sum(len(e['blocks']) for e in entries),
                       year_nc_bytes=sum(b['file']['size'] for e in entries for b in e['blocks']),
                       verification='previously_verified_extraction_receipts',output_files_reopened=0,
                       value_counts_recomputed=False)
        atomic_json(root/'index/validation_summary.json',summary)
        index = dict(schema='station-cf-v1',status='COMPLETED',prepared=io.hashed_file(prepared),
                     prepared_identity=prep['identity'],entries=entries,summary=summary)
        index['identity']=digest(index)
        atomic_json(root/'index/authoritative_index.json',index)
        return root/'index/authoritative_index.json'
