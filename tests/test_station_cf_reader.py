from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import numpy as np
import pytest
from station_cf_fixtures import build
from prepare_station_cf import prepare
from extract_station_cf import extract_unit,publish
from station_cf_reader import open_station_cf,iter_station_cf
from grid_cf_io import read_json


def test_counterfactual_reader_and_iteration(tmp_path):
    cfg=build(tmp_path,climates=('ssp126','ssp245'),calendar='proleptic_gregorian')
    prepared=prepare(cfg,tmp_path/'out',sample=True);prep=read_json(prepared)
    for task in prep['tasks']:
        extract_unit(prepared,model=task['model'],climate_scenario=task['climate_scenario'],
                     station_scenario=task['station_scenario'],tech=task['tech'],patch=task['source_patch'],processes=1)
    index=publish(prepared)
    opts=dict(model='M',climate_scenario='ssp126',station_scenario='ssp585',tech='wind')
    ds=open_station_cf(index,**opts)
    assert ds.sizes['station']==6 and not isinstance(ds.wind_cf.variable._data,np.ndarray)
    missing=np.flatnonzero(ds.mapping_status.values==3)
    assert np.isnan(ds.wind_cf.isel(time=slice(0,3),station=missing).values).all()
    order=[5,0,3,1];ids=ds.station_id.values[order].tolist()
    part=open_station_cf(index,**opts,station_ids=ids,years=(2016,2016))
    assert part.station_id.values.tolist()==ids and part.sizes['time']==2928
    expected=ds.wind_cf.isel(time=slice(2920,None),station=order).values
    np.testing.assert_array_equal(part.wind_cf.values,expected)
    np.testing.assert_array_equal(ds.wind_cf.isel(time=[2,0,1],station=[1,0]).values,
                                  ds.wind_cf.values[np.ix_([2,0,1],[1,0])])
    other=open_station_cf(index,**{**opts,'station_scenario':'ssp126'})
    # Different Station IDs, identical values at identical catalog locations.
    for i in range(ds.sizes['station']):
        j=np.flatnonzero((other.lon.values==ds.lon.values[i]) & (other.lat.values==ds.lat.values[i]))[0]
        np.testing.assert_array_equal(ds.wind_cf.isel(station=i).values,other.wind_cf.isel(station=j).values)
    changed=open_station_cf(index,**{**opts,'climate_scenario':'ssp245'})
    assert changed.station_id.values.tolist()==ds.station_id.values.tolist()
    point=np.flatnonzero((ds.lon.values==.25)&(ds.lat.values==-.5))[0]
    assert ds.wind_cf.isel(time=2,station=point).item()!=changed.wind_cf.isel(time=2,station=point).item()
    with pytest.raises(ValueError,match='unknown'):open_station_cf(index,**opts,station_ids=['bad'])
    chunks=list(iter_station_cf(index,**opts,years=(2015,2015),time_chunk=1000,station_chunk=2))
    assert sum(x.wind_cf.size for x in chunks)==2920*6
    assert max(x.sizes['station'] for x in chunks)<=2
