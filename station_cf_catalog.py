"""Validated Station catalogs and cumulative capacity snapshots."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from grid_cf_io import atomic_json, digest, file_identity, read_json, partial_path

SCENARIOS = ('ssp126', 'ssp245', 'ssp585')
TECHS = ('wind', 'solar')
LABELS = dict(zip(SCENARIOS, ('SSP1-2.6', 'SSP2-4.5', 'SSP5-6.0')))
SNAPSHOTS = (2030, 2040, 2050)


def scenario(value, *, station=False):
    value = 'ssp585' if station and value == 'ssp560' else value
    if value not in SCENARIOS:
        raise ValueError(f'unsupported {"Station" if station else "Climate"}: {value}')
    return value


def lon180(values):
    return (np.asarray(values, dtype='f8') + 180) % 360 - 180


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_csv(frame, path):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = partial_path(path)
    try:
        frame.to_csv(tmp, index=False, float_format='%.17g',
                     compression={'method': 'gzip', 'mtime': 0} if path.suffix == '.gz' else None)
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def read_csv(path):
    return pd.read_csv(path, keep_default_na=False, float_precision='round_trip',
                       dtype={'station_id': str})


def validated_rows(path, report_dir):
    frame = pd.read_csv(path, float_precision='round_trip')
    columns = ['year', 'type', 'lon', 'lat', 'capacity_gw']
    if not set(columns) <= set(frame):
        raise ValueError(f'missing station columns: {set(columns) - set(frame)}')
    frame = frame[columns].copy()
    frame['source_row'] = np.arange(len(frame)) + 2
    frame['type'] = frame['type'].astype(str).str.strip().str.lower()
    for name in ('year', 'lon', 'lat', 'capacity_gw'):
        frame[name] = pd.to_numeric(frame[name], errors='coerce')
    bad = (~np.isfinite(frame[['year', 'lon', 'lat', 'capacity_gw']])).any(axis=1)
    bad |= (~frame.type.isin(TECHS) | ~frame.lat.between(-90, 90)
            | ~frame.lon.between(-180, 360) | (frame.capacity_gw < 0)
            | ~frame.year.isin(SNAPSHOTS))
    if bad.any():
        write_csv(frame.loc[bad], Path(report_dir) / 'invalid_rows.csv.gz')
        raise ValueError('invalid station rows; see invalid_rows.csv.gz')
    frame['lon'] = lon180(frame.lon)
    frame['year'] = frame.year.astype('int16')
    duplicate = frame.duplicated(['year', 'type', 'lon', 'lat'], keep=False)
    if duplicate.any():
        write_csv(frame.loc[duplicate], Path(report_dir) / 'duplicate_rows.csv.gz')
        raise ValueError('duplicate snapshot position')
    if set(frame.year) != set(SNAPSHOTS):
        raise ValueError('missing capacity snapshot year')
    return frame


def build_catalog(csv_path, station_scenario, output_root, *, reference_csv):
    """Compare against an upstream snapshot export, retaining every capacity row."""
    station_scenario = scenario(station_scenario, station=True)
    csv_path, reference_csv = Path(csv_path).resolve(), Path(reference_csv).resolve()
    if csv_path == reference_csv:
        raise ValueError('reference_csv must be an independent upstream snapshot export')
    source_hash, reference_hash = sha(csv_path), sha(reference_csv)
    identity = digest({'schema': 'station-catalog-v1', 'station': station_scenario,
                       'source_sha256': source_hash, 'reference_sha256': reference_hash,
                       'implementation': sha(__file__), 'capacity_semantics': 'snapshot_total'})
    directory = Path(output_root) / 'catalogs' / station_scenario / identity
    directory.mkdir(parents=True, exist_ok=True)
    done = directory / 'catalog.json'
    if done.exists():
        result = read_json(done)
        if all(sha(p) == h for p, h in result['files'].items()):
            return result
        raise ValueError('published catalog changed')
    frame = validated_rows(csv_path, directory)
    reference = validated_rows(reference_csv, directory / 'reference')
    keys = ['year', 'type', 'lon', 'lat']
    left, right = [x.sort_values(keys).reset_index(drop=True) for x in (frame, reference)]
    if not left[keys].equals(right[keys]) or not np.array_equal(
            left.capacity_gw.to_numpy(), right.capacity_gw.to_numpy()):
        raise ValueError('station CSV differs from upstream snapshot export')
    result = {'identity': identity, 'station_scenario': station_scenario,
              'station_source_label': LABELS[station_scenario],
              'capacity_semantics': 'snapshot_total', 'snapshot_years': list(SNAPSHOTS),
              'source': file_identity(csv_path), 'source_sha256': source_hash,
              'reference': file_identity(reference_csv), 'reference_sha256': reference_hash,
              'raw_rows': len(frame), 'catalogs': {}, 'files': {}}
    capacities = []
    for tech in TECHS:
        raw = frame.loc[frame.type == tech]
        sites = raw.groupby(['lon', 'lat'], sort=True, as_index=False).agg(first_snapshot_year=('year', 'min'))
        sites['station_id'] = [hashlib.sha1(f'{station_scenario}|{tech}|{x:.4f}|{y:.4f}'.encode()).hexdigest()[:20]
                               for x, y in zip(sites.lon, sites.lat)]
        if sites.station_id.duplicated().any():
            raise ValueError('station_id coordinate rounding collision')
        sites = sites.sort_values('station_id').reset_index(drop=True)
        path = directory / tech / 'stations.csv.gz'
        write_csv(sites, path)
        result['catalogs'][tech] = {'path': str(path.resolve()), 'count': len(sites), 'sha256': sha(path)}
        result['files'][str(path.resolve())] = sha(path)
        rows = raw.merge(sites[['lon', 'lat', 'station_id']], on=['lon', 'lat'], validate='many_to_one')
        rows['station_scenario'] = station_scenario
        rows = rows.rename(columns={'type': 'tech'})
        capacities.append(rows)
    path = directory / 'capacity_rows.csv.gz'
    write_csv(pd.concat(capacities).sort_values('source_row'), path)
    result['capacity_path'] = str(path.resolve()); result['files'][str(path.resolve())] = sha(path)
    if sha(csv_path) != source_hash or sha(reference_csv) != reference_hash:
        raise ValueError('station source changed during preparation')
    atomic_json(done, result)
    return result


def capacity_snapshot(catalog, tech, snapshot_year, station_ids=None):
    """Return GW for exactly one published snapshot, with catalog-known absences zero."""
    if catalog['capacity_semantics'] != 'snapshot_total' or snapshot_year not in catalog['snapshot_years']:
        raise ValueError('unknown capacity snapshot year or semantics')
    if tech not in TECHS:
        raise ValueError('unknown tech')
    for path, expected in catalog['files'].items():
        if sha(path) != expected:
            raise ValueError('catalog file changed')
    sites = read_csv(catalog['catalogs'][tech]['path'])
    ids = sites.station_id.tolist() if station_ids is None else list(station_ids)
    if len(set(ids)) != len(ids) or not set(ids) <= set(sites.station_id):
        raise ValueError('duplicate/unknown station IDs')
    rows = read_csv(catalog['capacity_path'])
    rows = rows.loc[(rows.tech == tech) & (rows.year == snapshot_year)]
    return rows.set_index('station_id').capacity_gw.reindex(ids, fill_value=0.0)
