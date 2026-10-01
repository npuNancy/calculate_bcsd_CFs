"""Lazy Climate × Station access and bounded iteration over published CF."""
from __future__ import annotations

import netCDF4 as nc
import numpy as np
import xarray as xr
from xarray.backends import BackendArray
from xarray.core import indexing

from grid_cf_io import FILL, digest, read_json
from station_cf_catalog import scenario, sha
from station_cf_mapping import read_mapping
import station_cf_io as io


class StationCFArray(BackendArray):
    """Open files only for an explicit outer time/station selection."""
    def __init__(self, entries, source_blocks, ids, time_indices, tech):
        self.dtype = np.dtype('f4'); self.shape = (len(time_indices), len(ids))
        self.entries = entries; self.time_indices = np.asarray(time_indices)
        self.tech = tech; self.station_positions = {}
        self.starts = np.r_[0, np.cumsum([b['time']['count'] for b in source_blocks])]
        positions = {sid: i for i, sid in enumerate(ids)}
        for entry in entries:
            if not entry['blocks']:
                continue
            io.check_file(entry['manifest'])
            first = entry['blocks'][0]; io.check_file(first['file']); io.check_file(first['sidecar'])
            with nc.Dataset(first['path']) as ds:
                file_ids = list(ds['station_id'][:])
            chosen = [(positions[sid], i) for i, sid in enumerate(file_ids) if sid in positions]
            for global_pos, local_pos in chosen:
                if global_pos in self.station_positions:
                    raise ValueError('station appears in multiple source patches')
                self.station_positions[global_pos] = (entry['source_patch'], local_pos)

    def __getitem__(self, key):
        return indexing.explicit_indexing_adapter(key, self.shape, indexing.IndexingSupport.OUTER, self._read)

    def _read(self, key):
        time_key, station_key = key
        tids = np.atleast_1d(np.arange(self.shape[0])[time_key])
        sids = np.atleast_1d(np.arange(self.shape[1])[station_key])
        raw_times = self.time_indices[tids]
        out = np.full((len(tids), len(sids)), np.nan, dtype='f4')
        for entry in self.entries:
            matches = [(j, self.station_positions[int(s)][1]) for j, s in enumerate(sids)
                       if int(s) in self.station_positions and self.station_positions[int(s)][0] == entry['source_patch']]
            if not matches:
                continue
            target_stations, local_stations = map(np.asarray, zip(*matches))
            # netCDF4 outer indexing accepts ordered integer vectors. Restore arbitrary caller order.
            unique_stations, station_inverse = np.unique(local_stations, return_inverse=True)
            for b, block in enumerate(entry['blocks']):
                where = np.flatnonzero((raw_times >= self.starts[b]) & (raw_times < self.starts[b+1]))
                if not len(where):
                    continue
                io.check_file(block['file']); io.check_file(block['sidecar'])
                local_times = raw_times[where]-self.starts[b]
                unique_times, time_inverse = np.unique(local_times, return_inverse=True)
                with nc.Dataset(block['path']) as ds:
                    ds.set_auto_maskandscale(False)
                    values = np.asarray(ds[self.tech+'_cf'][unique_times, unique_stations], dtype='f4')
                values = values[time_inverse][:, station_inverse]
                values[values == FILL] = np.nan
                out[np.ix_(where, target_stations)] = values
        if np.isscalar(time_key):
            out = out[0]
        if np.isscalar(station_key):
            out = out[..., 0]
        return out


def open_station_cf(index_path, *, model, climate_scenario, station_scenario, tech, station_ids=None, years=None):
    climate_scenario = scenario(climate_scenario); station_scenario = scenario(station_scenario, station=True)
    index = read_json(index_path)
    if index['status'] != 'COMPLETED' or index['identity'] != digest({k:v for k,v in index.items() if k != 'identity'}):
        raise ValueError('invalid authoritative index')
    io.check_file(index['prepared']); prep = io.load_prepared(index['prepared']['path'])
    if index['prepared_identity'] != prep['identity']:
        raise ValueError('index/prepared identity mismatch')
    keys = ('model', 'climate_scenario', 'station_scenario', 'tech')
    wanted = (model, climate_scenario, station_scenario, tech)
    entries = [e for e in index['entries'] if tuple(e[k] for k in keys) == wanted]
    if not entries:
        raise ValueError('unknown Climate × Station combination')
    mapping = prep['mappings'][entries[0]['mapping_key']]
    if sha(mapping['path']) != mapping['sha256']:
        raise ValueError('mapping changed')
    frame = read_mapping(mapping['path']).set_index('station_id', drop=False)
    ids = frame.index.tolist() if station_ids is None else list(station_ids)
    if len(set(ids)) != len(ids) or not set(ids) <= set(frame.index):
        raise ValueError('unknown/duplicate station IDs in selected Station catalog')
    frame = frame.loc[ids]
    source_blocks = prep['sources'][entries[0]['source_key']]['blocks']
    dates = []
    for b in source_blocks:
        io.check_file(b['file'])
        with nc.Dataset(b['path']) as ds:
            ds.set_auto_maskandscale(False)
            if io.time_meta(ds) != b['time']:
                raise ValueError('reference time axis changed')
            dates.extend(nc.num2date(ds['time'][:], ds['time'].units, ds['time'].calendar,
                                    only_use_cftime_datetimes=True))
    dates = np.asarray(dates, dtype=object); time_indices = np.arange(len(dates))
    if years is not None:
        if len(years) != 2 or years[0] > years[1]:
            raise ValueError('invalid requested year range')
        chosen = np.array([years[0] <= t.year <= years[1] for t in dates])
        dates, time_indices = dates[chosen], time_indices[chosen]
    backend = StationCFArray(entries, source_blocks, ids, time_indices, tech)
    covered = set(np.flatnonzero(frame.mapping_status.isin([0, 1]).to_numpy()))
    if set(backend.station_positions) != covered:
        raise ValueError('published station coverage mismatch')
    variable = xr.Variable(('time', 'station'), indexing.LazilyIndexedArray(backend),
                           attrs={'units': '1', 'coordinates': 'station_id lon lat'})
    coordinates = {'time': dates, 'station': np.arange(len(ids))}
    coordinates.update({k: ('station', frame[k].to_numpy()) for k in frame})
    result = xr.Dataset({tech+'_cf': variable}, coords=coordinates,
                        attrs={**dict(zip(keys, wanted)), 'capacity_semantics': 'snapshot_total'})
    result.time.encoding.update(units=source_blocks[0]['time']['units'], calendar=source_blocks[0]['time']['calendar'])
    return result


def iter_station_cf(index_path, *, model, climate_scenario, station_scenario, tech,
                    station_ids=None, years=None, time_chunk=240, station_chunk=1024):
    """Yield loaded bounded arrays; uncovered locations form their own patch group."""
    if time_chunk <= 0 or station_chunk <= 0:
        raise ValueError('chunk sizes must be positive')
    ds = open_station_cf(index_path, model=model, climate_scenario=climate_scenario,
                         station_scenario=station_scenario, tech=tech, station_ids=station_ids, years=years)
    patches = ds.source_patch.values
    for patch in sorted(set(patches)):
        positions = np.flatnonzero(patches == patch)
        for s in range(0, len(positions), station_chunk):
            for t in range(0, ds.sizes['time'], time_chunk):
                yield ds.isel(station=positions[s:s+station_chunk], time=slice(t, t+time_chunk)).load()
