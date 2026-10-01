"""Small native-grid CF products, independent of BCSD readers/physics."""
from pathlib import Path
import json

import cftime
import netCDF4 as nc
import numpy as np
import pandas as pd

import grid_cf_io as gio
from station_cf_catalog import sha


def build(root, *, blocks=2, climates=('ssp126',), stations=('ssp126','ssp585'), techs=('wind',), calendar='noleap'):
    root = Path(root).resolve(); grid = root/'grid'; grid.mkdir(parents=True)
    patch_meta = {'patches': {'A': {'core_bbox_360': [0, 1, -1, 1]}, 'B': {'core_bbox_360': [1, 2, -1, 1]}}}
    gio.atomic_json(grid/'inputs/patch_manifest.json', patch_meta)
    gio.atomic_json(grid/'inputs/release.json', {'run_id': 'sample-grid'})
    gio.atomic_json(grid/'runtime/completion.json', {'status': 'COMPLETED', 'run_id': 'sample-grid',
                                                   'release_sha256': sha(grid/'inputs/release.json')})
    lat = np.array([-.5, .5], dtype='f8')
    for ci, climate in enumerate(climates):
        for tech in techs:
            for pi, patch in enumerate(('A','B')):
                lon = np.array([.25, .75], dtype='f8') + pi
                mask = np.ones((2,2), dtype='i1')
                if patch == 'A': mask[0,1] = 0
                prov = {'model':'M','scenario':climate,'tech':tech,'patch':patch,'years':f'2015-{2014+blocks}','grid_fingerprint':gio.array_digest(lat,lon,mask)}
                identity = gio.digest(prov)
                entries = []
                for i in range(blocks):
                    year=2015+i; units='hours since 2015-01-01'; offset=1.5 if tech=='solar' else 0
                    low=nc.date2num(cftime.datetime(year,1,1,calendar=calendar),units,calendar)
                    high=nc.date2num(cftime.datetime(year+1,1,1,calendar=calendar),units,calendar)
                    times=np.arange(low,high,3,dtype='f8')+offset
                    values=np.broadcast_to(np.array([[0,.2],[.7,1]],dtype='f4'),(len(times),2,2)).copy()
                    values[:,mask==0]=gio.FILL
                    values[1,0,0]=gio.FILL
                    values[2,0,0]=np.float32(.1+ci*.1+pi*.01)
                    path=grid/'outputs'/'M'/climate/patch/tech/'blocks'/identity/f'cf_{year}-{year}.nc'
                    path.parent.mkdir(parents=True,exist_ok=True)
                    bid=gio.digest({'unit':identity,'block_index':i})
                    with nc.Dataset(path,'w') as ds:
                        for name,arr in (('time',times),('lat',lat),('lon',lon)):
                            ds.createDimension(name,len(arr));ds.createVariable(name,arr.dtype,(name,))[:]=arr
                        ds['time'].units=units;ds['time'].calendar=calendar
                        ds.createVariable('domain_mask','i1',('lat','lon'))[:]=mask
                        v=ds.createVariable(tech+'_cf','f4',('time','lat','lon'),fill_value=gio.FILL,zlib=True,chunksizes=(240,2,2))
                        v.units='1';v[:]=values
                        ds.setncatts({'identity':bid,'model':'M','scenario':climate,'patch_id':patch,'tech':tech,'grid_fingerprint':prov['grid_fingerprint']})
                    gio.atomic_json(str(path)+'.json',{'status':'COMPLETED','identity':bid,'file':gio.file_identity(path),
                                                       'provenance':prov,'block_index':i,'start_year':year,'end_year':year,'time_count':len(times)})
                    entries.append({'index':i,'start_year':year,'end_year':year,'path':str(path),'status':'COMPLETED'})
                gio.atomic_json(grid/'outputs'/'M'/climate/patch/tech/'manifest.json',
                                {'identity':identity,'provenance':prov,'status':'COMPLETED','blocks':entries})
    station_root=root/'stations';station_root.mkdir(); refs=root/'reference';refs.mkdir()
    files={};references={}
    for station in stations:
        rows=[]
        for year,capacity in ((2030,1.),(2040,1.5),(2050,2.)):
            for tech in techs:
                for x,y in ((.25,-.5),(.251,-.5),(.75,-.5),(1.,-.5),(1.25,.5),(4.25,-.5)):
                    rows.append(dict(year=year,type=tech,lon=x,lat=y,capacity_gw=capacity))
        name=f'{station}.csv';pd.DataFrame(rows).to_csv(station_root/name,index=False)
        (refs/name).write_bytes((station_root/name).read_bytes())
        files[station]=name;references[station]=str(refs/name)
    config={'grid_cf_root':str(grid),'stations_root':str(station_root),'station_files':files,
            'station_reference_files':references,'models':['M'],'climate_scenarios':list(climates),
            'station_scenarios':list(stations),'techs':list(techs),'years':f'2015-{2014+blocks}',
            'spatial_method':'nearest_regular','max_distance_deg':.6}
    path=root/'config.json';gio.atomic_json(path,config)
    return path
