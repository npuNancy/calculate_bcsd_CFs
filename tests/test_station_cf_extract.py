from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import json
import subprocess
import netCDF4 as nc
import numpy as np
import pytest

from station_cf_fixtures import build
from prepare_station_cf import prepare,enumerate_tasks,MODELS
from station_cf_catalog import SCENARIOS
from extract_station_cf import extract_unit,publish
from grid_cf_io import read_json,atomic_json,FILL


def run(prepared,station='ssp585',**kwargs):
    return extract_unit(prepared,model='M',climate_scenario='ssp126',station_scenario=station,
                        tech='wind',patch='A',**kwargs)


def test_production_combinations():
    config={'models':MODELS,'climate_scenarios':SCENARIOS,'station_scenarios':SCENARIOS,'techs':['wind','solar'],'years':'2015-2060'}
    tasks=enumerate_tasks(config,{f'P{i}':{} for i in range(47)})
    assert len(tasks)==len({t['unit_id'] for t in tasks})==3384
    assert sum(t['climate_scenario']==t['station_scenario'] for t in tasks)==1128


def test_extraction_resume_and_orphan(tmp_path):
    cfg=build(tmp_path);prepared=prepare(cfg,tmp_path/'out',sample=True)
    path=run(prepared,processes=1,time_chunk=37,station_chunk=2)
    state=read_json(path);block=state['blocks'][0];first=Path(block['path']);before=first.stat().st_mtime_ns
    with nc.Dataset(first) as ds:
        ds.set_auto_maskandscale(False);values=ds['wind_cf'][:]
        assert np.all(values[:,ds['domain_mask'][:]==0]==FILL)
        ids=np.flatnonzero(ds['domain_mask'][:]==1)
        assert values[0,ids[0]]==0 and values[1,ids[0]]==FILL
        assert np.array_equal(values[:,ids[0]],values[:,ids[1]])
        assert np.all(ds['valid_count'][:]+ds['missing_count'][:]==values.shape[0])
        assert ds.climate_scenario=='ssp126' and ds.station_scenario=='ssp585'
    run(prepared,processes=2,time_chunk=37,station_chunk=2)
    assert first.stat().st_mtime_ns==before
    Path(str(first)+'.json').unlink()
    run(prepared,processes=1,time_chunk=37,station_chunk=2)
    assert first.exists() and read_json(path)['blocks'][0]['path']!=str(first)
    with pytest.raises(ValueError,match='profile'):run(prepared,processes=1,time_chunk=38)


@pytest.mark.parametrize('processes',[1,2,4,8])
def test_real_cli_spawn(tmp_path,processes):
    cfg=build(tmp_path,blocks=8,stations=('ssp585',));prepared=prepare(cfg,tmp_path/'out',sample=True)
    subprocess.run([sys.executable,'extract_station_cf.py','--prepared',str(prepared),
                    '--model','M','--climate-scenario','ssp126','--station-scenario','ssp560',
                    '--tech','wind','--patch','A','--processes',str(processes)],check=True,capture_output=True,text=True)
    path=tmp_path/'out/outputs/M/climate_ssp126/station_ssp585/A/wind/manifest.json'
    state=read_json(path);assert state['status']=='COMPLETED' and len(state['blocks'])==8
    for b in state['blocks']:
        with nc.Dataset(b['path']) as ds:
            ds.set_auto_maskandscale(False)
            good=np.flatnonzero(ds['domain_mask'][:]==1)
            np.testing.assert_array_equal(ds['wind_cf'][2,good],np.full(len(good),.1,dtype='f4'))


def test_frozen_source_change_and_scope(tmp_path):
    cfg=build(tmp_path)
    with pytest.raises(ValueError,match='production scope'):prepare(cfg,tmp_path/'bad')
    prepared=prepare(cfg,tmp_path/'out',sample=True)
    with pytest.raises(ValueError):publish(prepared)
    p=Path(next(iter(read_json(prepared)['sources'].values()))['blocks'][0]['path'])
    with nc.Dataset(p,'a') as ds:ds.changed='yes'
    with pytest.raises(ValueError,match='changed'):run(prepared,processes=1)


def test_solar_empty_unit_and_direct_oracle(tmp_path):
    import pandas as pd
    cfg=build(tmp_path,techs=('solar',),calendar='proleptic_gregorian',stations=('ssp585',))
    config=read_json(cfg)
    # Keep every snapshot, but remove all B-owned stations so that B is an empty unit.
    for p in (Path(config['stations_root'])/'ssp585.csv',Path(config['station_reference_files']['ssp585'])):
        frame=pd.read_csv(p);frame=frame.loc[frame.lon!=1.25];frame.to_csv(p,index=False)
    prepared=prepare(cfg,tmp_path/'out',sample=True)
    prep=read_json(prepared)
    for t in prep['tasks']:
        path=extract_unit(prepared,model='M',climate_scenario='ssp126',station_scenario='ssp585',
                          tech='solar',patch=t['source_patch'],processes=2)
        m=read_json(path)
        if t['source_patch']=='B':
            assert m['status']=='EMPTY_NO_STATIONS' and not m['blocks']
            continue
        source=prep['sources'][t['source_key']]
        for b,out in zip(source['blocks'],m['blocks']):
            with nc.Dataset(b['path']) as src,nc.Dataset(out['path']) as dst:
                src.set_auto_maskandscale(False);dst.set_auto_maskandscale(False)
                np.testing.assert_array_equal(src['time'][:],dst['time'][:])
                assert dst['time'][0]%3==1.5
                for i,(y,x) in enumerate(zip(dst['source_iy'][:],dst['source_ix'][:])):
                    np.testing.assert_array_equal(dst['solar_cf'][:,i],src['solar_cf'][:,int(y),int(x)])
    index=read_json(publish(prepared,processes=1))
    assert read_json(publish(prepared,processes=2))==index
    assert index['summary']['empty_units']==1 and index['summary']['year_nc_files']==2


def test_invalid_values_fail_and_output_lock(tmp_path):
    from grid_cf_io import output_lock,file_identity
    from extract_station_cf import unit_directory
    cfg=build(tmp_path,stations=('ssp585',))
    config=read_json(cfg)
    manifest=Path(config['grid_cf_root'])/'outputs/M/ssp126/A/wind/manifest.json'
    source=read_json(manifest);p=source['blocks'][0]['path']
    with nc.Dataset(p,'a') as ds:ds['wind_cf'][4,0,0]=2.0
    side=read_json(p+'.json');side['file']=file_identity(p);atomic_json(p+'.json',side)
    prepared=prepare(cfg,tmp_path/'out',sample=True);prep=read_json(prepared)
    directory=unit_directory(prep,prep['tasks'][0])
    with output_lock(directory/'.unit.lock'):
        with pytest.raises(RuntimeError,match='locked'):run(prepared,processes=1)
    with pytest.raises(ValueError,match='outside'):run(prepared,processes=1)
    assert read_json(directory/'manifest.json')['status']=='FAILED'


def test_preparation_rejects_time_and_sidecar_identity(tmp_path):
    from grid_cf_io import file_identity
    cfg=build(tmp_path,stations=('ssp585',));config=read_json(cfg)
    p=Path(config['grid_cf_root'])/'outputs/M/ssp126/B/wind/manifest.json'
    meta=read_json(p);block=Path(meta['blocks'][0]['path'])
    side=read_json(str(block)+'.json');side['identity']='wrong';atomic_json(str(block)+'.json',side)
    with pytest.raises(ValueError,match='sidecar'):prepare(cfg,tmp_path/'out',sample=True)


def test_cross_patch_time_mismatch(tmp_path):
    from grid_cf_io import file_identity
    cfg=build(tmp_path,stations=('ssp585',));config=read_json(cfg)
    m=read_json(Path(config['grid_cf_root'])/'outputs/M/ssp126/B/wind/manifest.json')
    for b in m['blocks']:
        p=b['path']
        with nc.Dataset(p,'a') as ds:ds['time'][:]=ds['time'][:]+1
        side=read_json(p+'.json');side['file']=file_identity(p);atomic_json(p+'.json',side)
    with pytest.raises(ValueError,match='cross-patch time'):prepare(cfg,tmp_path/'out',sample=True)


def test_parallel_publish_failure_preserves_index(tmp_path):
    cfg=build(tmp_path,stations=('ssp585',));prepared=prepare(cfg,tmp_path/'out',sample=True)
    prep=read_json(prepared)
    for t in prep['tasks']:
        extract_unit(prepared,model='M',climate_scenario='ssp126',station_scenario='ssp585',
                     tech=t['tech'],patch=t['source_patch'],processes=1)
    path=publish(prepared,processes=2);before=path.read_bytes()
    entry=read_json(path)['entries'][0]
    sidecar=Path(entry['blocks'][0]['sidecar']['path'])
    side=read_json(sidecar);side['identity']='corrupt';atomic_json(sidecar,side)
    with pytest.raises(ValueError,match='incomplete year output'):
        publish(prepared,processes=2)
    assert path.read_bytes()==before
    with pytest.raises(ValueError,match='positive'):
        publish(prepared,processes=0)
