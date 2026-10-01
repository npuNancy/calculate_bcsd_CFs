"""Read-only grid CF contracts and bounded station CF extraction."""
from __future__ import annotations

import os
from pathlib import Path
import resource
import subprocess
import time

import netCDF4 as nc
import numpy as np

from grid_cf_io import (FILL, array_digest, atomic_json, digest, file_identity,
                        output_lock, partial_path, read_json)
from station_cf_catalog import LABELS, sha
from station_cf_mapping import read_mapping

MODULES = ('grid_cf_io.py', 'station_cf_catalog.py', 'station_cf_mapping.py', 'station_cf_io.py',
           'prepare_station_cf.py', 'extract_station_cf.py', 'station_cf_reader.py')


def implementation():
    root = Path(__file__).resolve().parent
    try:
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True,
                                         stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = 'unavailable'
    return {'code_sha': commit, 'source_sha256': {name: sha(root/name) for name in MODULES}}


def check_file(record):
    if file_identity(record['path']) != {k: record[k] for k in ('path', 'size', 'mtime_ns')}:
        raise ValueError(f'file changed: {record["path"]}')
    if 'sha256' in record and sha(record['path']) != record['sha256']:
        raise ValueError(f'file hash changed: {record["path"]}')


def hashed_file(path):
    return {**file_identity(path), 'sha256': sha(path)}


def resolve(value, root):
    path = Path(value)
    return (path if path.is_absolute() else Path(root)/path).resolve()


def axis(ds, name):
    v = ds[name]; values = np.asarray(v[:])
    if v.dimensions != (name,) or not len(values) or not np.isfinite(values).all():
        raise ValueError(f'invalid {name} axis')
    if len(values) > 1 and not (np.all(np.diff(values) > 0) or np.all(np.diff(values) < 0)):
        raise ValueError(f'nonmonotonic {name}')
    return values


def time_meta(ds):
    v = ds['time']; raw = axis(ds, 'time')
    units, calendar = v.units, v.calendar
    dates = nc.num2date(raw, units, calendar, only_use_cftime_datetimes=True)
    hours = np.asarray(nc.date2num(dates, 'hours since 1970-01-01', calendar), dtype='f8')
    if len(hours) > 1 and not np.allclose(np.diff(hours), 3, rtol=0, atol=1e-7):
        raise ValueError('source CF is not a continuous 3-hour time axis')
    return {'sha256': array_digest(raw), 'dtype': raw.dtype.str, 'count': len(raw),
            'units': units, 'calendar': calendar, 'first_hour': float(hours[0]), 'last_hour': float(hours[-1]),
            'start_year': dates[0].year, 'end_year': dates[-1].year}


def source_schema(ds, tech):
    v = ds[tech+'_cf']
    if v.dimensions != ('time', 'lat', 'lon') or v.dtype != np.dtype('f4') or v.units != '1':
        raise ValueError('source CF variable schema mismatch')
    if any(key in v.ncattrs() for key in ('scale_factor', 'add_offset')) or v._FillValue != FILL:
        raise ValueError('unexpected source CF encoding')
    lat, lon = axis(ds, 'lat'), axis(ds, 'lon')
    mask = np.asarray(ds['domain_mask'][:])
    if ds['domain_mask'].dimensions != ('lat', 'lon') or mask.shape != (len(lat), len(lon)) or not np.isin(mask, [0, 1]).all():
        raise ValueError('source mask schema mismatch')
    return {'lat': lat, 'lon': lon, 'mask': mask}


def scan_source(manifest_path, *, model, climate_scenario, tech, patch, years, expected_blocks=None):
    manifest_path = Path(manifest_path).resolve(); meta = read_json(manifest_path)
    prov = meta['provenance']
    if meta['status'] != 'COMPLETED' or digest(prov) != meta['identity']:
        raise ValueError('source manifest incomplete/identity mismatch')
    for key, value in dict(model=model, scenario=climate_scenario, tech=tech, patch=patch, years=years).items():
        if prov[key] != value:
            raise ValueError(f'source manifest {key} mismatch')
    blocks = meta['blocks']
    if not blocks or (expected_blocks is not None and len(blocks) != expected_blocks):
        raise ValueError('wrong number of CF year blocks')
    files = [hashed_file(manifest_path)]; entries = []; spatial = None; previous = None; spatial_chunk = None
    for i, block in enumerate(blocks):
        if block['index'] != i or block['status'] != 'COMPLETED':
            raise ValueError('source blocks incomplete/unordered')
        path = resolve(block['path'], manifest_path.parent)
        side_path = Path(str(path)+'.json'); side = read_json(side_path)
        identity = digest({'unit': meta['identity'], 'block_index': i})
        if side['status'] != 'COMPLETED' or side['identity'] != identity or side['file'] != file_identity(path):
            raise ValueError('source CF sidecar/stat mismatch')
        if side['provenance'] != prov or side['block_index'] != i:
            raise ValueError('source CF provenance mismatch')
        with nc.Dataset(path) as ds:
            ds.set_auto_maskandscale(False)
            for k, val in dict(model=model, scenario=climate_scenario, tech=tech, patch_id=patch, identity=identity).items():
                if ds.getncattr(k) != val:
                    raise ValueError(f'source NC {k} mismatch')
            grid = source_schema(ds, tech)
            fingerprint = array_digest(grid['lat'], grid['lon'], grid['mask'])
            if ds.grid_fingerprint != fingerprint or prov['grid_fingerprint'] != fingerprint:
                raise ValueError('source declared grid fingerprint mismatch')
            chunks = ds[tech+'_cf'].chunking()
            current_chunk = chunks[1:] if isinstance(chunks, list) else [64, 64]
            if spatial_chunk is not None and current_chunk != spatial_chunk:
                raise ValueError('source spatial chunk layout changes between years')
            spatial_chunk = current_chunk
            if spatial is None:
                spatial = grid
            elif fingerprint != array_digest(spatial['lat'], spatial['lon'], spatial['mask']):
                raise ValueError('source grid/mask changes between years')
            tm = time_meta(ds)
            for k in ('start_year', 'end_year'):
                if tm[k] != block[k] or side[k] != block[k]:
                    raise ValueError('year block contract mismatch')
            if side['time_count'] != tm['count']:
                raise ValueError('source sidecar time count mismatch')
            if previous and (tm['start_year'] != previous['end_year'] + 1 or
                             tm['calendar'] != previous['calendar'] or
                             not np.isclose(tm['first_hour'] - previous['last_hour'], 3, atol=1e-7, rtol=0)):
                raise ValueError('source time blocks gap/overlap/calendar mismatch')
            previous = tm
        record = file_identity(path); files.extend([record, hashed_file(side_path)])
        entries.append({'index': i, 'path': str(path), 'identity': identity, 'time': tm,
                        'start_year': block['start_year'], 'end_year': block['end_year'], 'file': record})
    if years != f'{entries[0]["start_year"]}-{entries[-1]["end_year"]}':
        raise ValueError('analysis years mismatch')
    for record in files:
        check_file(record)
    return {'model': model, 'climate_scenario': climate_scenario, 'tech': tech, 'source_patch': patch,
            'identity': meta['identity'], 'manifest': str(manifest_path), 'files': files,
            'grid_fingerprint': fingerprint, 'spatial_chunk': spatial_chunk, 'blocks': entries}, spatial


def load_prepared(path):
    prep = read_json(path)
    identity = prep['identity']
    if identity != digest({k: v for k, v in prep.items() if k != 'identity'}):
        raise ValueError('prepared identity mismatch')
    for record in prep['assets']:
        check_file(record)
    return prep


def check_source(source):
    for record in source['files']:
        check_file(record)


def selected_mapping(task, mapping_frame=None):
    mapping = task['mapping']
    if mapping_frame is None:
        if sha(mapping['path']) != mapping['sha256']:
            raise ValueError('mapping changed')
        mapping_frame = read_mapping(mapping['path'])
    frame = mapping_frame
    frame = frame.loc[(frame.source_patch == task['source_patch']) & frame.mapping_status.isin([0, 1])]
    return frame.sort_values(['source_iy', 'source_ix', 'station_id']).reset_index(drop=True)


def sort_mapping(frame, source_path, tech):
    with nc.Dataset(source_path) as ds:
        chunk = ds[tech+'_cf'].chunking()
    cy, cx = chunk[1:] if isinstance(chunk, list) else (64, 64)
    frame = frame.assign(chunk_y=frame.source_iy//cy, chunk_x=frame.source_ix//cx)
    return frame.sort_values(['chunk_y', 'chunk_x', 'source_iy', 'source_ix', 'station_id']).reset_index(drop=True)


def clean_values(values, variable):
    out = np.asarray(values, dtype='f4').copy()
    missing = ~np.isfinite(out)
    for key in ('_FillValue', 'missing_value'):
        for value in np.atleast_1d(getattr(variable, key, FILL)):
            missing |= out == value
    if np.any(((out < 0) | (out > 1)) & ~missing):
        raise ValueError('source CF outside [0,1]')
    out[missing] = FILL
    return out


def block_id(identity, block):
    return digest({'unit': identity, 'source': block['identity'], 'time': block['time'], 'index': block['index']})


def create_station_output(path, source_ds, frame, task, identity, unit_identity, profile):
    ds = nc.Dataset(path, 'w', format='NETCDF4')
    try:
        n = len(frame); nt = len(source_ds.dimensions['time'])
        ds.createDimension('time', nt); ds.createDimension('station', n)
        t = ds.createVariable('time', source_ds['time'].dtype, ('time',))
        t[:] = source_ds['time'][:]; t.units = source_ds['time'].units; t.calendar = source_ds['time'].calendar
        t.standard_name = 'time'
        ds.createVariable('station', 'i4', ('station',))[:] = np.arange(n)
        for name in frame:
            if name in ('chunk_y', 'chunk_x', 'source_patch', 'station_patch'):
                continue
            values = frame[name].to_numpy()
            dtype = (str if values.dtype.kind in 'OUS' else 'i1' if name in ('mapping_status', 'domain_mask')
                     else 'i2' if name == 'first_snapshot_year' else 'i4' if values.dtype.kind in 'iu'
                     else 'f4' if name.startswith('match_dist_') else 'f8')
            v = ds.createVariable(name, dtype, ('station',)); v[:] = values.astype(object) if dtype is str else values
        ds['station_id'].cf_role = 'timeseries_id'; ds.featureType = 'timeSeries'
        for k in ('lon', 'source_grid_lon'):
            ds[k].units = 'degrees_east'
        for k in ('lat', 'source_grid_lat'):
            ds[k].units = 'degrees_north'
        ds['mapping_status'].flag_values = np.array([0, 1], dtype='i1')
        ds['mapping_status'].flag_meanings = 'MATCHED OUTSIDE_DOMAIN'
        ds['match_dist_deg'].units = 'degrees'; ds['match_dist_km'].units = 'km'
        v = ds.createVariable(task['tech']+'_cf', 'f4', ('time', 'station'), fill_value=FILL,
                              zlib=profile['compress_level'] > 0, complevel=profile['compress_level'], shuffle=True,
                              chunksizes=(min(profile['time_chunk'], nt), min(profile['station_chunk'], n)))
        v.units = '1'; v.valid_range = np.array([0, 1], dtype='f4'); v.coordinates = 'station_id lon lat'
        ds.createVariable('valid_count', 'i8', ('station',))
        ds.createVariable('missing_count', 'i8', ('station',))
        attrs = {k: task[k] for k in ('model', 'climate_scenario', 'station_scenario', 'tech', 'source_patch')}
        attrs.update(code_sha=task['code_sha'], source_file=source_ds.filepath(),
                     implementation_sha256=task['implementation_sha256'],
                     time_chunk=profile['time_chunk'], station_chunk=profile['station_chunk'],
                     compress_level=profile['compress_level'], schema_version='station-cf-v1', identity=identity, block_identity=identity,
                     unit_identity=unit_identity, source_run_id=task['source_run_id'],
                     source_cf_identity=task['source']['identity'],
                     source_block_identity=source_ds.identity, station_source_label=LABELS[task['station_scenario']],
                     mapping_identity=task['mapping']['identity'], catalog_identity=task['catalog_identity'],
                     match_method='nearest_regular', max_distance_deg=task['mapping']['max_distance_deg'],
                     capacity_semantics='snapshot_total', capacity_time_mask='off', years=task['years'])
        ds.setncatts(attrs)
        return ds
    except BaseException:
        ds.close(); raise


def validate_station(path, task, block, frame, identity, *, sample=False):
    with nc.Dataset(path) as ds, nc.Dataset(block['path']) as source:
        ds.set_auto_maskandscale(False); source.set_auto_maskandscale(False)
        if ds.identity != identity or ds.schema_version != 'station-cf-v1':
            raise ValueError('station output identity mismatch')
        for k in ('model', 'climate_scenario', 'station_scenario', 'tech', 'source_patch'):
            if ds.getncattr(k) != task[k]:
                raise ValueError('station output context mismatch')
        if time_meta(ds) != block['time']:
            raise ValueError('station time mismatch')
        for k in frame:
            if k in ('chunk_y', 'chunk_x', 'source_patch', 'station_patch'):
                continue
            if not np.array_equal(ds[k][:], frame[k].to_numpy()):
                raise ValueError(f'station coordinate mismatch: {k}')
        v = ds[task['tech']+'_cf']
        if v.dimensions != ('time', 'station') or v.dtype != np.dtype('f4') or v._FillValue != FILL:
            raise ValueError('station CF schema mismatch')
        valid, missing = ds['valid_count'][:], ds['missing_count'][:]
        if np.any(valid < 0) or np.any(missing < 0) or not np.all(valid+missing == len(ds.dimensions['time'])):
            raise ValueError('invalid sample counts')
        if np.any(valid[frame.domain_mask.to_numpy() == 0] != 0):
            raise ValueError('domain-outside station has valid values')
        if sample and len(frame):
            rng = np.random.default_rng(20261001)
            stations = np.unique(np.r_[0, len(frame)-1, rng.integers(len(frame), size=min(32, len(frame)))])
            times = np.unique(np.r_[0, v.shape[0]-1, rng.integers(v.shape[0], size=min(32, v.shape[0]))])
            for s in stations:
                row = frame.iloc[s]
                expected = (clean_values(source[task['tech']+'_cf'][times, int(row.source_iy), int(row.source_ix)], source[task['tech']+'_cf'])
                            if row.domain_mask == 1 else np.full(len(times), FILL, dtype='f4'))
                if not np.array_equal(v[times, int(s)], expected):
                    raise ValueError('independent source sample mismatch')
        return {'valid_count': int(valid.sum()), 'missing_count': int(missing.sum()),
                'station_count': len(frame), 'time_count': len(ds.dimensions['time'])}


def completed(path, task, block, frame, identity):
    try:
        side = read_json(str(path)+'.json')
        if side['status'] != 'COMPLETED' or side['identity'] != identity or side['file'] != file_identity(path):
            return False
        counts = validate_station(path, task, block, frame, identity)
        if any(side[key] != value for key, value in counts.items()):
            return False
        if side['source_block_identity'] != block['identity']:
            return False
        return True
    except (OSError, ValueError, KeyError, AttributeError, RuntimeError, IndexError):
        return False


def extract_block(payload):
    task, block, profile, unit_identity, filename = payload
    begin = time.monotonic(); path = Path(filename); path.parent.mkdir(parents=True, exist_ok=True)
    identity = block_id(unit_identity, block)
    frame = sort_mapping(selected_mapping(task), task['source']['blocks'][0]['path'], task['tech'])
    check_source(task['source'])
    with output_lock(str(path)+'.lock'):
        if completed(path, task, block, frame, identity):
            return {'index': block['index'], 'path': str(path), 'identity': identity, 'status': 'COMPLETED'}
        if path.exists() or Path(str(path)+'.json').exists():
            raise ValueError('incomplete output requires a new attempt path')
        tmp = partial_path(path); timing = {'read': 0.0, 'gather_write': 0.0}
        valid_counts = np.zeros(len(frame), dtype='i8'); minimum, maximum = None, None
        read_bytes = 0
        try:
            with nc.Dataset(block['path']) as src:
                src.set_auto_maskandscale(False)
                grid = source_schema(src, task['tech'])
                if array_digest(grid['lat'], grid['lon'], grid['mask']) != task['source']['grid_fingerprint']:
                    raise ValueError('source grid fingerprint mismatch')
                if time_meta(src) != block['time'] or src.identity != block['identity']:
                    raise ValueError('source block metadata mismatch')
                if len(frame) and (np.any(frame.source_iy < 0) or np.any(frame.source_iy >= len(grid['lat'])) or
                                   np.any(frame.source_ix < 0) or np.any(frame.source_ix >= len(grid['lon']))):
                    raise ValueError('mapping indices out of bounds')
                yi = frame.source_iy.to_numpy(dtype='i8'); xi = frame.source_ix.to_numpy(dtype='i8')
                from station_cf_catalog import lon180
                if (not np.array_equal(grid['mask'][yi, xi], frame.domain_mask.to_numpy()) or
                        not np.array_equal(grid['lat'][yi], frame.source_grid_lat.to_numpy()) or
                        not np.array_equal(lon180(grid['lon'][xi]), frame.source_grid_lon.to_numpy())):
                    raise ValueError('mapping does not match source geometry/domain')
                v = src[task['tech']+'_cf']; chunk = v.chunking()
                cy, cx = chunk[1:] if isinstance(chunk, list) else (64, 64)
                v.set_var_chunk_cache(size=16*1024*1024, nelems=1009, preemption=0.75)
                groups = frame.groupby([frame.source_iy//cy, frame.source_ix//cx], sort=True).indices
                with create_station_output(tmp, src, frame, task, identity, unit_identity, profile) as dst:
                    out = dst[task['tech']+'_cf']; out.set_var_chunk_cache(size=16*1024*1024, nelems=1009, preemption=0.75)
                    for start in range(0, len(src.dimensions['time']), profile['time_chunk']):
                        stop = min(start+profile['time_chunk'], len(src.dimensions['time']))
                        for (iyc, ixc), positions in groups.items():
                            rows = frame.iloc[positions]
                            tile = None
                            if (rows.domain_mask == 1).any():
                                mark = time.monotonic(); y0, x0 = int(iyc)*cy, int(ixc)*cx
                                tile = v[start:stop, y0:y0+cy, x0:x0+cx]
                                read_bytes += tile.nbytes; timing['read'] += time.monotonic()-mark
                            for off in range(0, len(positions), profile['station_chunk']):
                                mark = time.monotonic(); pos = positions[off:off+profile['station_chunk']]
                                subset = frame.iloc[pos]; values = np.full((stop-start, len(pos)), FILL, dtype='f4')
                                keep = subset.domain_mask.to_numpy() == 1
                                if keep.any():
                                    coordinates = np.column_stack((subset.source_iy.to_numpy()[keep]-y0, subset.source_ix.to_numpy()[keep]-x0))
                                    unique, inv = np.unique(coordinates, axis=0, return_inverse=True)
                                    gathered = clean_values(tile[:, unique[:, 0], unique[:, 1]], v)
                                    values[:, keep] = gathered[:, inv]
                                good = values != FILL; valid_counts[pos] += good.sum(axis=0)
                                if good.any():
                                    lo, hi = float(values[good].min()), float(values[good].max())
                                    minimum = lo if minimum is None else min(minimum, lo)
                                    maximum = hi if maximum is None else max(maximum, hi)
                                # The mapping order makes every source tile a contiguous station run.
                                out[start:stop, int(pos[0]):int(pos[-1])+1] = values
                                timing['gather_write'] += time.monotonic()-mark
                    dst['valid_count'][:] = valid_counts
                    dst['missing_count'][:] = len(src.dimensions['time'])-valid_counts
            validate_station(tmp, task, block, frame, identity, sample=True)
            check_source(task['source'])
            tmp.replace(path)
            side = {'status': 'COMPLETED', 'identity': identity, 'file': file_identity(path),
                    'source_block_identity': block['identity'], 'station_count': len(frame),
                    'time_count': block['time']['count'], 'valid_count': int(valid_counts.sum()),
                    'missing_count': int(len(frame)*block['time']['count']-valid_counts.sum()),
                    'minimum': minimum, 'maximum': maximum, 'read_array_bytes': read_bytes,
                    'timing_seconds': {**timing, 'total': time.monotonic()-begin},
                    'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    'audit': 'all streamed values range checked; independent fixed-seed source samples checked'}
            atomic_json(str(path)+'.json', side)
        finally:
            tmp.unlink(missing_ok=True)
    return {'index': block['index'], 'path': str(path), 'identity': identity, 'status': 'COMPLETED'}
