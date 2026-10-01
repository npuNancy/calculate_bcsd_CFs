from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest
from station_cf_mapping import nearest,map_stations,STATUS


def grids():
    return {'A':{'lat':np.array([-.5,.5]),'lon':np.array([.25,.75]),'mask':np.array([[1,0],[1,1]],dtype='i1')},
            'B':{'lat':np.array([.5,-.5]),'lon':np.array([1.75,1.25]),'mask':np.ones((2,2),dtype='i1')}}, {
            'A':{'core_bbox_360':[0,1,-1,1]},'B':{'core_bbox_360':[1,2,-1,1]}}


def test_nearest_ties_and_dateline():
    a=np.array([179.95,-179.95])
    assert a[nearest(a,[180],periodic=True)][0]==-179.95
    assert nearest([1,0], [.5]).tolist()==[1]
    assert nearest([1,0],[.5001]).tolist()==[0]


def test_global_owner_missing_domain_and_pole():
    g,p=grids();s=pd.DataFrame({'lon':[.25,1.,1.25,4.25,.25], 'lat':[-.5,-.5,.5,-.5,90.]})
    result=map_stations(s,g,p,.6)
    assert result.mapping_status.tolist()==[0,1,0,3,3]
    assert result.iloc[1].source_patch=='A' and result.iloc[1].station_patch=='B'
    assert result.iloc[2].source_iy==0 and result.iloc[2].source_ix==1
    assert result.iloc[4].source_grid_lat==pytest.approx(89.5)
    assert map_stations(s,g,p,.01).iloc[1].mapping_status==STATUS['TOO_FAR']


def test_lattice_and_duplicate_ownership():
    g,p=grids();g['B']['lon'][1]=1.3
    with pytest.raises(ValueError,match='regular|lattice'):map_stations(pd.DataFrame({'lon':[0.],'lat':[0.]}),g,p)


def test_float32_axis_quantization():
    from station_cf_mapping import theoretical_axis
    g={'P':{'lat':np.linspace(30.05,59.95,300,dtype='f4')}}
    axis=theoretical_axis(g,'lat')
    assert len(axis)==1800
    np.testing.assert_allclose(axis[[0,-1]],[-89.95,89.95],atol=2e-5)


def test_exact_grid_coordinate_on_patch_boundary():
    g,p=grids()
    g['A']['lon']=np.array([0.,.5]);g['B']['lon']=np.array([1.5,1.])
    result=map_stations(pd.DataFrame({'lon':[1.],'lat':[.5]}),g,p,.15)
    assert result.iloc[0].mapping_status==STATUS['MATCHED']
    assert result.iloc[0].source_patch=='B' and result.iloc[0].source_ix==1
