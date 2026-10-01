"""Global regular-axis nearest mapping with explicit patch and domain coverage."""
from __future__ import annotations

from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd

from grid_cf_io import array_digest, atomic_json, digest, partial_path
from station_cf_catalog import lon180, read_csv, sha

STATUS = {'MATCHED': 0, 'OUTSIDE_DOMAIN': 1, 'TOO_FAR': 2, 'NO_SOURCE_PATCH': 3}
TIE_TOL = 1e-8


def nearest(axis, points, *, periodic=False):
    """Return original indices; ties use ascending normalized coordinate."""
    axis = lon180(axis) if periodic else np.asarray(axis, dtype='f8')
    order = np.argsort(axis); sorted_axis = axis[order]
    if len(axis) == 0 or np.any(np.diff(sorted_axis) <= 0):
        raise ValueError('empty or duplicate coordinate axis')
    points = lon180(points) if periodic else np.asarray(points, dtype='f8')
    at = np.searchsorted(sorted_axis, points)
    low = (at - 1) % len(axis) if periodic else np.clip(at - 1, 0, len(axis)-1)
    high = at % len(axis) if periodic else np.clip(at, 0, len(axis)-1)
    def distance(v):
        return abs((v - points + 180) % 360 - 180) if periodic else abs(v - points)
    dl, dh = distance(sorted_axis[low]), distance(sorted_axis[high])
    choose_low = (dl < dh - TIE_TOL) | ((abs(dl-dh) <= TIE_TOL) & (sorted_axis[low] <= sorted_axis[high]))
    return order[np.where(choose_low, low, high)]


def theoretical_axis(grids, name):
    raw_arrays = [np.asarray(g[name]) for g in grids.values()]
    arrays = [a.astype('f8') for a in raw_arrays]
    tolerance = max([1e-7] + [4*np.finfo(a.dtype).eps*max(1., float(abs(a).max()))
                             for a in raw_arrays if a.dtype.kind == 'f'])
    differences = [abs(np.diff(a)) for a in arrays if len(a) > 1]
    if not differences:
        raise ValueError(f'cannot establish regular {name} spacing')
    steps = np.concatenate(differences)
    # Fit the whole axis so float32 coordinate quantization does not accumulate globally.
    step = float(np.median([abs(np.polynomial.polynomial.polyfit(np.arange(len(a)), a, 1)[1])
                            for a in arrays if len(a) > 1]))
    if step <= 0 or not np.allclose(steps, step, rtol=1e-5, atol=tolerance):
        raise ValueError(f'nonregular {name} spacing')
    origin = float(arrays[0][0])
    for a in arrays:
        if not np.allclose(a-origin, np.rint((a-origin)/step)*step, atol=tolerance, rtol=0):
            raise ValueError(f'inconsistent {name} lattice origin')
    lower, upper = (-90, 90) if name == 'lat' else (-180, 180)
    if name == 'lon' and not np.isclose(360/step, round(360/step), atol=1e-5):
        raise ValueError('longitude lattice is not periodic')
    first = int(np.ceil((lower-origin)/step - 1e-7))
    last = int(np.floor((upper-origin)/step + 1e-7))
    axis = origin + np.arange(first, last+1) * step
    if name == 'lon':
        axis = axis[axis < 180 - 1e-7]
    return axis


def inside(lon, lat, bbox):
    west, east, south, north = bbox
    return (((np.asarray(lon) - west) % 360 < east-west)
            & (np.asarray(lat) >= south)
            & ((np.asarray(lat) <= north) if north == 90 else (np.asarray(lat) < north)))


def validate_grids(grids, patches):
    for patch, g in grids.items():
        if patch not in patches:
            raise ValueError('unknown source patch')
        y, x, mask = np.asarray(g['lat']), np.asarray(g['lon']), np.asarray(g['mask'])
        if mask.shape != (len(y), len(x)) or not np.isin(mask, [0, 1]).all():
            raise ValueError('invalid domain mask')
        if not inside(x, np.full(len(x), y[0]), patches[patch]['core_bbox_360']).all() or not inside(
                np.full(len(y), x[0]), y, patches[patch]['core_bbox_360']).all():
            raise ValueError('source coordinates outside patch core')
    keys = sorted(grids)
    for i, key in enumerate(keys):
        a = grids[key]
        for other in keys[i+1:]:
            b = grids[other]
            if np.intersect1d(a['lat'], b['lat']).size and np.intersect1d(lon180(a['lon']), lon180(b['lon'])).size:
                raise ValueError('duplicate core grid ownership')


def map_stations(sites, grids, patches, max_distance_deg=0.15):
    if not np.isfinite(max_distance_deg) or max_distance_deg < 0:
        raise ValueError('invalid max distance')
    validate_grids(grids, patches)
    lat_axis, lon_axis = (theoretical_axis(grids, k) for k in ('lat', 'lon'))
    lat, lon = sites.lat.to_numpy(), sites.lon.to_numpy()
    gy = lat_axis[nearest(lat_axis, lat)]; gx = lon_axis[nearest(lon_axis, lon, periodic=True)]
    dx = abs((gx-lon+180) % 360-180); dy = abs(gy-lat)
    result = sites.copy()
    result['source_grid_lat'] = gy; result['source_grid_lon'] = lon180(gx)
    result['source_patch'] = ''; result['station_patch'] = ''
    result['source_iy'] = np.int32(-1); result['source_ix'] = np.int32(-1)
    result['domain_mask'] = np.int8(-1)
    result['mapping_status'] = np.int8(STATUS['NO_SOURCE_PATCH'])
    result['match_dist_deg'] = np.maximum(dx, dy)
    a = np.sin(np.radians(gy-lat)/2)**2 + np.cos(np.radians(gy))*np.cos(np.radians(lat))*np.sin(np.radians(dx)/2)**2
    result['match_dist_km'] = 6371.0088 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
    station_owned = np.zeros(len(sites), dtype=bool)
    for patch in sorted(grids):
        g = grids[patch]; bbox = patches[patch]['core_bbox_360']
        s = inside(lon, lat, bbox)
        if (station_owned & s).any():
            raise ValueError('overlapping station patch cores')
        station_owned |= s; result.loc[s, 'station_patch'] = patch
        # Expand only candidate lookup, not ownership: actual coordinate equality below
        # resolves roundoff at an exactly shared patch boundary.
        west, east, south, north = bbox
        delta = (gx-west) % 360
        candidate = np.flatnonzero(((delta <= east-west+2e-5) | (delta >= 360-2e-5))
                                   & (gy >= south-2e-5) & (gy <= north+2e-5))
        if not len(candidate):
            continue
        iy = nearest(g['lat'], gy[candidate]); ix = nearest(g['lon'], gx[candidate], periodic=True)
        exact = ((abs(np.asarray(g['lat'])[iy]-gy[candidate]) < 2e-5)
                 & (abs((np.asarray(g['lon'])[ix]-gx[candidate]+180) % 360-180) < 2e-5))
        candidate, iy, ix = candidate[exact], iy[exact], ix[exact]
        if (result.loc[candidate, 'source_patch'] != '').any():
            raise ValueError('ambiguous grid owner')
        mask = np.asarray(g['mask'])[iy, ix]
        result.loc[candidate, 'source_patch'] = patch
        result.loc[candidate, 'source_iy'] = iy.astype('i4'); result.loc[candidate, 'source_ix'] = ix.astype('i4')
        result.loc[candidate, 'source_grid_lat'] = np.asarray(g['lat'])[iy]
        result.loc[candidate, 'source_grid_lon'] = lon180(np.asarray(g['lon'])[ix])
        result.loc[candidate, 'domain_mask'] = mask
        result.loc[candidate, 'mapping_status'] = np.where(mask == 1, STATUS['MATCHED'], STATUS['OUTSIDE_DOMAIN']).astype('i1')
    result.loc[result.match_dist_deg > max_distance_deg + TIE_TOL, 'mapping_status'] = STATUS['TOO_FAR']
    return result


def grid_fingerprints(grids):
    return {p: array_digest(np.asarray(g['lat']), np.asarray(g['lon']), np.asarray(g['mask'])) for p, g in sorted(grids.items())}


def write_mapping(frame, path, identity):
    tmp = partial_path(path)
    try:
        with nc.Dataset(tmp, 'w') as ds:
            ds.createDimension('station', len(frame))
            ds.identity = identity
            for column in frame:
                values = frame[column].to_numpy()
                if values.dtype.kind in 'OUS':
                    v = ds.createVariable(column, str, ('station',)); v[:] = values.astype(object)
                else:
                    dtype = ('i1' if column in ('mapping_status', 'domain_mask') else 'i2' if column == 'first_snapshot_year'
                             else 'i4' if values.dtype.kind in 'iu' else 'f4' if column.startswith('match_dist_') else 'f8')
                    v = ds.createVariable(column, dtype, ('station',), zlib=True); v[:] = values
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def read_mapping(path):
    with nc.Dataset(path) as ds:
        ds.set_auto_mask(False)
        return pd.DataFrame({name: ds[name][:] for name in ds.variables})


def build_mapping(catalog, tech, grids, patches, output_root, max_distance_deg=0.15):
    entry = catalog['catalogs'][tech]
    if sha(entry['path']) != entry['sha256']:
        raise ValueError('catalog changed')
    identity = digest({'catalog': catalog['identity'], 'tech': tech, 'grids': grid_fingerprints(grids),
                       'patches': patches, 'method': 'nearest_regular', 'max_distance_deg': max_distance_deg,
                       'tie_tol': TIE_TOL, 'implementation': sha(__file__)})
    directory = Path(output_root) / 'mappings' / identity
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / 'mapping.nc'
    if path.exists():
        from grid_cf_io import read_json
        old = read_json(directory / 'coverage.json')
        if sha(path) != old['sha256']:
            raise ValueError('published mapping changed')
        return old
    frame = map_stations(read_csv(entry['path']), grids, patches, max_distance_deg)
    write_mapping(frame, path, identity)
    counts = {k: int((frame.mapping_status == v).sum()) for k, v in STATUS.items()}
    report = {'identity': identity, 'path': str(path.resolve()), 'sha256': sha(path),
              'catalog_identity': catalog['identity'], 'catalog_count': len(frame), 'counts': counts,
              'max_distance_deg': max_distance_deg,
              'cross_patch_count': int(((frame.source_patch != frame.station_patch) & (frame.source_patch != '')).sum()),
              'distance_bins': {str(limit): int((frame.match_dist_deg <= limit).sum()) for limit in (1e-6, 0.05, 0.15)},
              'patch_counts': {p: int(((frame.source_patch == p) & frame.mapping_status.isin([0, 1])).sum()) for p in grids},
              'unique_grid_counts': {p: int(len(frame.loc[(frame.source_patch == p) & frame.mapping_status.isin([0, 1]),
                                                          ['source_iy', 'source_ix']].drop_duplicates())) for p in grids}}
    atomic_json(directory / 'coverage.json', report)
    return report
