from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from station_cf_catalog import build_catalog, capacity_snapshot, read_csv, scenario


def source(tmp_path, change=None):
    rows=[dict(year=y,type='wind',lon=180.,lat=1.,capacity_gw=c) for y,c in ((2030,1.),(2040,1.5),(2050,2.))]
    frame=pd.DataFrame(rows)
    if change:frame=change(frame)
    path=tmp_path/'stations.csv';ref=tmp_path/'reference.csv'
    frame.to_csv(path,index=False);frame.to_csv(ref,index=False)
    return path,ref


def test_snapshot_identity_and_alias(tmp_path):
    p,r=source(tmp_path);catalog=build_catalog(p,'ssp560',tmp_path/'out',reference_csv=r)
    sites=read_csv(catalog['catalogs']['wind']['path'])
    assert len(sites)==1 and sites.lon.iloc[0]==-180
    assert sites.first_snapshot_year.iloc[0]==2030 and 'activation_year' not in sites
    assert capacity_snapshot(catalog,'wind',2040).iloc[0]==1.5
    assert catalog['station_scenario']=='ssp585'
    with pytest.raises(ValueError):capacity_snapshot(catalog,'wind',2035)
    with pytest.raises(ValueError):capacity_snapshot(catalog,'wind',2030,['unknown'])
    with pytest.raises(ValueError):scenario('ssp560')
    assert build_catalog(p,'ssp585',tmp_path/'out',reference_csv=r)==catalog


@pytest.mark.parametrize('change,match',[
    (lambda f: pd.concat([f,f.iloc[:1]]),'duplicate'),
    (lambda f: f.assign(capacity_gw=-1),'invalid'),
    (lambda f: f.loc[f.year!=2040],'missing'),
    (lambda f: pd.concat([f,f.assign(lon=180.000001)]),'collision'),
])
def test_catalog_rejects_invalid(tmp_path,change,match):
    p,r=source(tmp_path,change)
    with pytest.raises(ValueError,match=match):build_catalog(p,'ssp126',tmp_path/'out',reference_csv=r)


def test_reference_mismatch(tmp_path):
    p,r=source(tmp_path);f=pd.read_csv(r);f.loc[0,'capacity_gw']=3;f.to_csv(r,index=False)
    with pytest.raises(ValueError,match='differs'):build_catalog(p,'ssp126',tmp_path/'out',reference_csv=r)
