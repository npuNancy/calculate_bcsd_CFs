"""Parallel preparation and content-keyed time conversion preserve contracts."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import gzip
import netCDF4 as nc
import numpy as np
import pytest

from grid_cf_io import read_json
from prepare_station_cf import prepare
import station_cf_io as io
from station_cf_fixtures import build


def time_file(path, calendar='noleap', offset=0):
    with nc.Dataset(path, 'w') as ds:
        ds.createDimension('time', 100)
        v=ds.createVariable('time','f8',('time',))
        v.units='hours since 2015-01-01';v.calendar=calendar
        v[:]=np.arange(100)*3+offset
    return path


@pytest.mark.parametrize('calendar', ['noleap','gregorian','360_day'])
@pytest.mark.parametrize('offset', [0,1.5])
def test_time_cache_exact_and_isolated(tmp_path, monkeypatch, calendar, offset):
    path=time_file(tmp_path/'time.nc',calendar,offset)
    with nc.Dataset(path) as ds: expected=io.time_meta(ds)
    original=io.nc.num2date;calls=[]
    def convert(*args,**kwargs):
        calls.append(1);return original(*args,**kwargs)
    monkeypatch.setattr(io.nc,'num2date',convert)
    cache={}
    with nc.Dataset(path) as ds:
        first=io.time_meta(ds,cache=cache)
        assert first==expected
        first['count']=-1
        assert io.time_meta(ds,cache=cache)==expected
    assert len(calls)==1
    with nc.Dataset(path,'r+') as ds:ds['time'][:]=ds['time'][:]+3
    with nc.Dataset(path) as ds:
        changed=io.time_meta(ds,cache=cache)
        assert changed['first_hour']==expected['first_hour']+3
    assert len(calls)==2
    # A content change is checked even if all surrounding metadata is unchanged.
    with nc.Dataset(path,'r+') as ds:ds['time'][50]+=0.25
    with nc.Dataset(path) as ds:
        with pytest.raises(ValueError,match='continuous 3-hour'):
            io.time_meta(ds,cache=cache)


def test_cache_units_and_calendar_are_part_of_identity(tmp_path):
    path=time_file(tmp_path/'time.nc');cache={}
    with nc.Dataset(path) as ds:first=io.time_meta(ds,cache=cache)
    with nc.Dataset(path,'r+') as ds:ds['time'].calendar='gregorian'
    with nc.Dataset(path) as ds:second=io.time_meta(ds,cache=cache)
    assert first['calendar']!=second['calendar'] and len(cache)==2
    with nc.Dataset(path,'r+') as ds:ds['time'].units='days since 2015-01-01'
    with nc.Dataset(path) as ds:
        with pytest.raises(ValueError,match='continuous 3-hour'):
            io.time_meta(ds,cache=cache)


def test_parallel_preparation_matches_serial(tmp_path):
    cfg=build(tmp_path,climates=('ssp126','ssp585'),techs=('wind','solar'))
    serial=read_json(prepare(cfg,tmp_path/'serial',sample=True,processes=1))
    parallel=read_json(prepare(cfg,tmp_path/'parallel',sample=True,processes=2))
    for key in ('sources','tasks','storage_estimate','capacity_semantics','config_hash'):
        assert serial[key]==parallel[key],key
    assert serial['mappings'].keys()==parallel['mappings'].keys()
    for key in serial['mappings']:
        assert serial['mappings'][key]['sha256']==parallel['mappings'][key]['sha256']
    for key in serial['catalogs']:
        a={Path(p).name:gzip.decompress(Path(p).read_bytes()) for p in serial['catalogs'][key]['files']}
        b={Path(p).name:gzip.decompress(Path(p).read_bytes()) for p in parallel['catalogs'][key]['files']}
        assert a==b


def test_parallel_prepare_rejects_changed_source(tmp_path):
    cfg=build(tmp_path)
    source=next((tmp_path/'grid/outputs').rglob('*.nc'))
    with source.open('ab') as f:f.write(b'changed')
    with pytest.raises(ValueError,match='sidecar/stat mismatch'):
        prepare(cfg,tmp_path/'parallel',sample=True,processes=2)


@pytest.mark.parametrize('processes', [0,17])
def test_invalid_process_count(tmp_path, processes):
    with pytest.raises(ValueError,match='processes'):
        prepare(tmp_path/'absent.json',tmp_path/'output',processes=processes)
