"""Native-grid CF contracts, numerical equivalence, multiprocessing and recovery."""
import json
from pathlib import Path
import subprocess
import sys

import cftime
import netCDF4 as nc
import numpy as np
import pytest
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import grid_cf_io as io
import patchify_grid_cf as grid
from patchify_station_cf import _solar_points_block, _wind_points, _doy_hour


def build(root, calendar="noleap", years=2, reverse=False, shape=None):
    root = Path(root)
    lat = np.array([70., 80., 90.], dtype="f8")
    if reverse:
        lat = lat[::-1].copy()
    lon = np.array([0., 90., 180., 270.], dtype="f8")
    mask = np.array([[1, 1, 0, 0], [1, 0, 1, 0], [1, 1, 0, 0]], dtype="i1")
    if shape is not None:
        lat = np.linspace(70, 90, shape[0], dtype="f8")
        lon = np.linspace(0, 360, shape[1], endpoint=False, dtype="f8")
        mask = np.ones(shape, dtype="i1")
        mask[:, -1] = 0
    yy, xx = np.where(mask)
    land = root / "land.nc"
    root.mkdir(parents=True, exist_ok=True)
    xr.Dataset({"land_mask": (("lat", "lon"), mask), "point_y": ("point", yy), "point_x": ("point", xx)},
               coords={"lat": lat, "lon": lon}).to_netcdf(land)
    pm = root / "patches.json"
    pm.write_text(json.dumps({"patches": {"P": {"core_bbox_360": [0, 360, 70, 90]}}}))
    contracts = [{"start_year": 2024+i, "end_year": 2024+i} for i in range(years)]
    units = "hours since 2024-01-01 00:00:00"
    rng = np.random.default_rng(42)
    source = {}
    for variable in ("uas", "vas", "rsds", "tas"):
        paths, times_all, values_all = [], [], []
        for index, contract in enumerate(contracts):
            year = contract["start_year"]
            lo = nc.date2num(cftime.datetime(year, 1, 1, calendar=calendar), units, calendar)
            hi = nc.date2num(cftime.datetime(year+1, 1, 1, calendar=calendar), units, calendar)
            times = np.arange(lo, hi, 3.) + (1.5 if variable == "rsds" else 0)
            vals = rng.normal(size=(len(times), len(yy))).astype("f4")
            vals = {"uas": lambda v: v*4+8, "vas": lambda v: v*3,
                    "tas": lambda v: v*8+280, "rsds": lambda v: np.abs(v)*300}[variable](vals)
            vals[3, 0] = np.nan; vals[8, 1] = np.inf
            if variable == "uas":
                vals[15:18, 0] = [0, 25/(10**(1/7)), 30]
            if variable == "vas":
                vals[15:18, 0] = 0
            if variable == "rsds":
                vals[20, 0] = 0
            p = root / "blocks" / "M" / "s" / variable / "P" / f"{variable}_{index}.nc"
            p.parent.mkdir(parents=True, exist_ok=True)
            ds = xr.Dataset({variable: (("time", "point"), vals), "y_index": ("point", yy),
                             "x_index": ("point", xx), "flat_index": ("point", yy*len(lon)+xx),
                             "global_point": ("point", np.arange(len(yy))),
                             "lat": ("point", lat[yy]), "lon": ("point", lon[xx])},
                            coords={"time": times}, attrs={"patch_id": "P"})
            ds.time.attrs.update(units=units, calendar=calendar)
            ds[variable].attrs["units"] = {"tas": "K", "rsds": "W m-2"}.get(variable, "m s-1")
            ds.to_netcdf(p, encoding={variable: {"zlib": True, "chunksizes": (120, min(128, len(yy)) if shape else 3), "_FillValue": io.FILL}})
            side = {"kind": "global-step6-block", "variable": variable, "patch_ids": ["P"], "outputs": {"P": str(p)}}
            p.with_name("block_"+p.stem+".json").write_text(json.dumps(side))
            paths.append(str(p)); times_all.append(times); values_all.append(vals)
        manifest = root / "manifests" / f"M__s__{variable}__P.json"
        manifest.parent.mkdir(exist_ok=True)
        manifest.write_text(json.dumps({"kind": "global-bcsd-three-stage", "land_plan": str(land),
                                       "time_block_contract": {"blocks": contracts},
                                       "final_tasks": [{"patch_id": "P", "variable": variable, "block_files": paths}]}))
        source[variable] = (np.concatenate(times_all), np.concatenate(values_all))
    return {"root": root, "lat": lat, "lon": lon, "mask": mask, "yy": yy, "xx": xx,
            "source": source, "units": units, "calendar": calendar, "years": f"2024-{2023+years}"}


def arguments(data, output, tech="solar", **changes):
    a = grid.parser().parse_args(["--bcsd-root", str(data["root"]), "--model", "M", "--scenario", "s",
                                 "--patch", "P", "--patch-manifest", str(data["root"]/"patches.json"),
                                 "--land-plan", str(data["root"]/"land.nc"), "--tech", tech,
                                 "--years", data["years"], "--output-root", str(output)])
    a.processes = 1; a.tile_shape = [2, 2]; a.time_chunk = 480
    for key, value in changes.items():
        setattr(a, key, value)
    return a


def outputs(manifest):
    m = json.loads(Path(manifest).read_text())
    pieces = []
    for row in m["blocks"]:
        with xr.open_dataset(row["path"]) as ds:
            pieces.append(ds.load())
    return m, xr.concat(pieces, dim="time")


def reference(data, tech):
    src = data["source"]
    times = src[io.VARS[tech][0]][0]
    values = {}
    for name in io.VARS[tech]:
        t, v = src[name]
        a = xr.DataArray(np.where(np.isfinite(v), v, np.nan), dims=("time", "point"), coords={"time": t})
        arr = (a.values if np.array_equal(t, times) else a.interp(time=times).values).astype("f4")
        if name == "tas": arr -= np.float32(273.15)
        if name == "rsds": arr /= 1000
        values[name] = arr
    valid = np.logical_and.reduce([np.isfinite(v) for v in values.values()])
    if tech == "wind":
        out = _wind_points(values["uas"], values["vas"])
    else:
        dates = nc.num2date(times, data["units"], data["calendar"])
        doy, hour = _doy_hour(dates)
        out = _solar_points_block(*(values[v] for v in io.VARS[tech]), doy, hour,
                                  data["lat"][data["yy"]], ((data["lon"][data["xx"]]+180)%360)-180)
    result = np.full((len(times), *data["mask"].shape), np.nan, dtype="f4")
    result[:, data["yy"], data["xx"]] = np.where(valid, out, np.nan)
    return result


@pytest.mark.parametrize("calendar,reverse", [("noleap", False), ("proleptic_gregorian", True)])
@pytest.mark.parametrize("tech", ["wind", "solar"])
def test_grid_reference_and_encoding(tmp_path, calendar, reverse, tech):
    data = build(tmp_path/"input", calendar, reverse=reverse)
    args = arguments(data, tmp_path/"out", tech)
    m, got = outputs(grid.compute(args))
    np.testing.assert_allclose(got[tech+"_cf"].values, reference(data, tech), rtol=1e-6, atol=1e-6)
    assert np.array_equal(got.lat, data["lat"]) and np.array_equal(got.lon, data["lon"])
    assert m["status"] == "COMPLETED" and m["final"]["status"] == "not_requested"
    assert not (tmp_path/"out/M/s/P"/(tech+".nc")).exists()
    for row in m["blocks"]:
        with nc.Dataset(row["path"]) as ds:
            assert ds[tech+"_cf"].dtype == np.dtype("f4")
            assert ds[tech+"_cf"]._FillValue == io.FILL
            assert ds["time"].units == data["units"] and ds["time"].calendar == calendar
    assert np.isnan(got[tech+"_cf"].values[:, data["mask"] == 0]).all()


def test_merge_reuse_and_corrupt_block(tmp_path):
    data = build(tmp_path/"input")
    args = arguments(data, tmp_path/"out")
    m, original = outputs(grid.compute(args))
    paths = [Path(r["path"]) for r in m["blocks"]]
    stamps = [p.stat().st_mtime_ns for p in paths]
    args.merge_final = True
    merged = io.read_json(grid.compute(args))
    assert [p.stat().st_mtime_ns for p in paths] == stamps
    assert all(r["reused"] for r in merged["blocks"])
    with xr.open_dataset(merged["final"]["path"]) as ds:
        np.testing.assert_array_equal(ds.solar_cf.values, original.solar_cf.values)
    Path(str(paths[0])+".json").write_text('{bad')
    grid.compute(args)
    assert paths[0].stat().st_mtime_ns != stamps[0]
    assert paths[1].stat().st_mtime_ns == stamps[1]
    args.time_chunk += 1
    with pytest.raises(ValueError, match="different inputs/configuration"):
        grid.compute(args)


@pytest.mark.parametrize("workers", [2, 4, 8])
def test_parallel_and_chunk_equivalence(tmp_path, workers):
    data = build(tmp_path/"input", years=8)
    a = arguments(data, tmp_path/"serial", tech="wind", time_chunk=1000)
    _, reference_ds = outputs(grid.compute(a))
    a.output_root = str(tmp_path/"parallel"); a.processes = workers; a.tile_shape = [3, 4]; a.time_chunk = 2000
    _, got = outputs(grid.compute(a))
    np.testing.assert_array_equal(got.wind_cf, reference_ds.wind_cf)


@pytest.mark.parametrize("fault", ["units", "time", "indices", "calendar", "sidecar", "domain", "years"])
def test_reject_bad_inputs(tmp_path, fault):
    data = build(tmp_path/"input")
    a = arguments(data, tmp_path/"out")
    p = next((data["root"]/"blocks/M/s/uas/P").glob('uas_0.nc'))
    if fault == "sidecar":
        p.with_name("block_"+p.stem+".json").unlink()
    elif fault == "years":
        a.years = "2024-2024"
    elif fault == "domain":
        with nc.Dataset(data["root"]/"land.nc", "a") as ds: ds["land_mask"][0, 2] = 1
    else:
        with nc.Dataset(p, "a") as ds:
            if fault == "units": ds["uas"].units = ""
            if fault == "time": ds["time"][10] += 1
            if fault == "indices": ds["global_point"][0] = 2
            if fault == "calendar": ds["time"].calendar = "360_day"
    with pytest.raises((ValueError, FileNotFoundError)):
        grid.compute(a)


def test_locks_and_merge_recovery(tmp_path, monkeypatch):
    data = build(tmp_path/"input")
    a = arguments(data, tmp_path/"out", merge_final=True)
    merge = grid._merge
    def fail(*args): raise OSError("interrupted merge")
    monkeypatch.setattr(grid, "_merge", fail)
    with pytest.raises(OSError, match="interrupted merge"):
        grid.compute(a)
    m = io.read_json(tmp_path/"out/M/s/P/solar/manifest.json")
    assert m["status"] == "FAILED" and all(r["status"] == "COMPLETED" for r in m["blocks"])
    monkeypatch.setattr(grid, "_merge", merge)
    m = io.read_json(grid.compute(a))
    assert all(r["reused"] for r in m["blocks"])
    with io.output_lock(tmp_path/"out/M/s/P/.solar.lock"):
        with pytest.raises(RuntimeError, match="locked"): grid.compute(a)


def test_cli_defaults_and_spawn(tmp_path):
    data = build(tmp_path/"input")
    a = arguments(data, tmp_path/"out")
    argv = []
    for name in ("bcsd_root", "model", "scenario", "patch", "patch_manifest", "land_plan", "tech", "years", "output_root"):
        argv += ["--"+name.replace('_','-'), str(getattr(a, name))]
    parsed = grid.parser().parse_args(argv)
    assert parsed.processes == 8 and parsed.merge_final is False
    with pytest.raises(SystemExit): grid.parser().parse_args(argv+["--processes", "0"])
    subprocess.run([sys.executable, str(Path(grid.__file__)), *argv, "--time-chunk", "2000"], check=True, capture_output=True)
    subprocess.run([sys.executable, str(Path(grid.__file__)), *argv, "--time-chunk", "2000", "--merge-final"], check=True, capture_output=True)
    m = io.read_json(tmp_path/"out/M/s/P/solar/manifest.json")
    assert m["status"] == "COMPLETED" and all(r["reused"] for r in m["blocks"])


def test_relocated_manifest_and_explicit_units(tmp_path):
    data = build(tmp_path/"input")
    a = arguments(data, tmp_path/"out")
    for variable in io.VARS[a.tech]:
        manifest = data["root"]/"manifests"/f"M__s__{variable}__P.json"
        m = io.read_json(manifest)
        m["final_tasks"][0]["block_files"] = [p.replace(str(data["root"]), "/old/production_v1") for p in m["final_tasks"][0]["block_files"]]
        manifest.write_text(json.dumps(m))
        for p in (data["root"]/"blocks/M/s"/variable/"P").glob('*.nc'):
            with nc.Dataset(p, "a") as ds: ds[variable].delncattr("units")
    with pytest.raises(ValueError, match="units"): grid.compute(a)
    a.input_units = ["uas=m/s", "vas=m/s", "tas=K", "rsds=W/m2"]
    _, ds = outputs(grid.compute(a))
    np.testing.assert_allclose(ds.solar_cf, reference(data, a.tech), rtol=1e-6, atol=1e-6)
    a.input_units = ["tas=K", "tas=K"]
    with pytest.raises(ValueError, match="duplicate"): grid.compute(a)


def test_unit_conflict_and_sidecar_time_identity(tmp_path):
    data = build(tmp_path/"input")
    a = arguments(data, tmp_path/"out", input_units=["tas=degC"])
    with pytest.raises(ValueError, match="units conflict"): grid.compute(a)
    a.input_units = []
    p = next((data["root"]/"blocks/M/s/uas/P").glob('uas_0.nc'))
    side = p.with_name("block_"+p.stem+".json")
    m = io.read_json(side); m.update(time_start=0, time_stop=1); side.write_text(json.dumps(m))
    with pytest.raises(ValueError, match="time indices"): grid.compute(a)


def test_inputs_changed_during_computation(tmp_path, monkeypatch):
    data = build(tmp_path/"input")
    a = arguments(data, tmp_path/"out", tech="wind")
    original = grid.physics.compute_wind_cf_chunk
    changed = False
    def change(*args):
        nonlocal changed
        if not changed:
            changed = True
            p = data["root"]/"patches.json"
            p.write_text(p.read_text()+"\n")
        return original(*args)
    monkeypatch.setattr(grid.physics, "compute_wind_cf_chunk", change)
    with pytest.raises(ValueError, match="input changed"): grid.compute(a)
    assert not list((tmp_path/"out").rglob('*.nc'))


def test_same_axis_nan_does_not_contaminate_neighbor(tmp_path):
    data = build(tmp_path/"input")
    a = arguments(data, tmp_path/"out", tech="wind")
    a.input_units = {}
    plan = io.load_plan(a)
    reader = io.BlockReader(plan["blocks"]["uas"], "uas")
    try:
        times = plan["blocks"]["uas"][0]["numbers"]
        got = reader.read(times[2:5], np.array([0]))[:, 0]
        assert np.isfinite(got[[0, 2]]).all() and np.isnan(got[1])
    finally:
        reader.close()


@pytest.mark.parametrize("tech", ["wind", "solar"])
def test_same_source_final_station_reference(tmp_path, tech):
    import pandas as pd
    import patchify_station_cf as station
    data = build(tmp_path/"input")
    a = arguments(data, tmp_path/"out", tech=tech)
    _, got = outputs(grid.compute(a))
    for variable in io.VARS[tech]:
        times, values = data["source"][variable]
        # Compare finite grid results; the old station path can turn missing inputs into zero.
        full = np.full((len(times), *data["mask"].shape), np.nan, dtype="f4")
        full[:, data["yy"], data["xx"]] = values
        ds = xr.Dataset({variable: (("time", "lat", "lon"), full)},
                        coords={"time": times, "lat": data["lat"], "lon": data["lon"]})
        ds.time.attrs.update(units=data["units"], calendar=data["calendar"])
        ds[variable].attrs["units"] = {"tas": "K", "rsds": "W m-2"}.get(variable, "m s-1")
        path = data["root"]/"outputs/M/s"/variable/"P.nc"
        path.parent.mkdir(parents=True, exist_ok=True)
        ds.to_netcdf(path)
        Path(str(path)+'.json').write_text(json.dumps({"patch_id":"P","variable":variable}))
    stations = tmp_path/'stations.csv'
    pd.DataFrame({'lat':data['lat'][data['yy']], 'lon':data['lon'][data['xx']], 'type':tech}).to_csv(stations,index=False)
    args = station.parser().parse_args(['--bcsd-root',str(data['root']),'--model','M','--scenario','s',
                                       '--patch','P','--tech',tech,'--stations-csv',str(stations),
                                       '--output-root',str(tmp_path/'station')])
    with xr.open_dataset(station.compute(args)) as ref:
        grid_values = got[tech+'_cf'].values[:,data['yy'],data['xx']]
        expected = ref[tech+'_cf'].values
        valid = np.isfinite(grid_values)
        np.testing.assert_allclose(grid_values[valid],expected[valid],atol=1e-6,rtol=1e-6)


def test_negative_time_offset_and_input_change(tmp_path):
    data = build(tmp_path/'input')
    for variable in ('uas','vas','tas'):
        times, values = data['source'][variable]
        data['source'][variable] = (times+1.5, values)
        for p in (data['root']/'blocks/M/s'/variable/'P').glob('*.nc'):
            with nc.Dataset(p,'a') as ds: ds['time'][:] = ds['time'][:]+1.5
    times, values = data['source']['rsds']
    data['source']['rsds'] = (times-1.5,values)
    for p in (data['root']/'blocks/M/s/rsds/P').glob('*.nc'):
        with nc.Dataset(p,'a') as ds: ds['time'][:] = ds['time'][:]-1.5
    a=arguments(data,tmp_path/'out')
    _,got=outputs(grid.compute(a))
    np.testing.assert_allclose(got.solar_cf,reference(data,'solar'),rtol=1e-6,atol=1e-6)
    p=data['root']/'patches.json'
    p.write_text(p.read_text()+'\n')
    with pytest.raises(ValueError,match='different inputs/configuration'): grid.compute(a)


def publish_final_inputs(data):
    """Publish equivalent grid inputs with the production final sidecar contract."""
    for variable, (times, values) in data['source'].items():
        full = np.full((len(times), *data['mask'].shape), io.FILL, dtype='f4')
        full[:, data['yy'], data['xx']] = values
        path = data['root']/'outputs/M/s'/variable/(variable+'_M_s_P.nc')
        path.parent.mkdir(parents=True, exist_ok=True)
        with nc.Dataset(path, 'w') as ds:
            for name, vals in [('time', times), ('lat', data['lat']), ('lon', data['lon'])]:
                ds.createDimension(name, len(vals)); ds.createVariable(name, 'f8', (name,))[:] = vals
            ds['time'].setncatts({'units': data['units'], 'calendar': data['calendar']})
            v = ds.createVariable(variable, 'f4', ('time','lat','lon'), fill_value=io.FILL,
                                  zlib=True, chunksizes=(120, 2, 2))
            v[:] = full
            v.units = {'tas':'K','rsds':'W/m2'}.get(variable,'m/s')
            ds.patch_id = 'P'
        side = {'kind':'global-final-patch','patch_id':'P','variable':variable,'output':str(path),
                'time_count':len(times),'shape':list(full.shape),'land_point_count':int(data['mask'].sum())}
        Path(str(path)+'.json').write_text(json.dumps(side))
        mp = data['root']/'manifests'/f'M__s__{variable}__P.json'
        m = io.read_json(mp); m['final_tasks'][0]['output'] = str(path); mp.write_text(json.dumps(m))


@pytest.mark.parametrize('calendar,reverse', [('noleap',False),('proleptic_gregorian',True)])
@pytest.mark.parametrize('tech', ['wind','solar'])
@pytest.mark.parametrize('missing_kind', ['nc', 'sidecar'])
def test_final_fallback_matches_blocks(tmp_path, calendar, reverse, tech, missing_kind):
    data = build(tmp_path/'input', calendar=calendar, reverse=reverse)
    publish_final_inputs(data)
    a = arguments(data,tmp_path/'blocks',tech=tech)
    block_manifest, block_ds = outputs(grid.compute(a))
    assert block_manifest['provenance']['input_mode']=='blocks'
    # Any missing required input switches the entire unit.
    variable = io.VARS[tech][-1]
    name = f'{variable}_0.nc' if missing_kind == 'nc' else f'block_{variable}_0.json'
    (data['root']/f'blocks/M/s/{variable}/P'/name).unlink()
    a.output_root=str(tmp_path/'final')
    m, got = outputs(grid.compute(a))
    assert m['provenance']['input_mode']=='final'
    np.testing.assert_array_equal(got[tech+'_cf'],block_ds[tech+'_cf'])
    with nc.Dataset(m['blocks'][0]['path']) as ds: assert ds.input_mode=='final'
    assert all('/blocks/' not in e['path'] for e in m['provenance']['inputs'])
    a.input_mode='blocks'
    with pytest.raises(FileNotFoundError,match='missing BCSD block'):grid.compute(a)


@pytest.mark.parametrize('tech', ['wind','solar'])
def test_final_spawn_merge_and_streaming(tmp_path, tech):
    data=build(tmp_path/'input',years=8)
    publish_final_inputs(data)
    a=arguments(data,tmp_path/'serial',tech=tech,input_mode='final',time_chunk=2000,merge_final=True)
    _, serial=outputs(grid.compute(a))
    a.output_root=str(tmp_path/'parallel');a.processes=8
    m,parallel=outputs(grid.compute(a))
    np.testing.assert_array_equal(parallel[tech+'_cf'],serial[tech+'_cf'])
    with xr.open_dataset(m['final']['path']) as ds:np.testing.assert_array_equal(ds[tech+'_cf'],serial[tech+'_cf'])
    stamps=[Path(b['path']).stat().st_mtime_ns for b in m['blocks']]
    again=io.read_json(grid.compute(a))
    assert all(b['reused'] for b in again['blocks'])
    assert stamps==[Path(b['path']).stat().st_mtime_ns for b in again['blocks']]
    plan=io.load_plan(a);reader=io.input_reader(plan,io.VARS[tech][0])
    try:
        reader.read(plan['blocks'][io.VARS[tech][0]][0]['numbers'][0:2],np.array([0]))
        assert reader.read_bytes==2*4  # Two times, one cell; not the complete file.
    finally:reader.close()


@pytest.mark.parametrize('fault',['time','calendar','grid','units','sidecar','missing','coverage'])
def test_final_rejects_invalid_inputs(tmp_path,fault):
    data=build(tmp_path/'input');publish_final_inputs(data)
    a=arguments(data,tmp_path/'out',tech='wind',input_mode='final')
    p=data['root']/'outputs/M/s/uas/uas_M_s_P.nc'
    if fault=='sidecar':
        side=Path(str(p)+'.json');s=io.read_json(side);s['patch_id']='wrong';side.write_text(json.dumps(s))
    elif fault=='missing':p.unlink()
    else:
        with nc.Dataset(p,'a') as ds:
            if fault=='time':ds['time'][10]+=1
            elif fault=='calendar':ds['time'].calendar='360_day'
            elif fault=='grid':ds['lat'][:]=ds['lat'][:]+1
            elif fault=='units':ds['uas'].units='K'
            elif fault=='coverage':ds['time'][:]=ds['time'][:]+6
    with pytest.raises((ValueError,FileNotFoundError,OSError)):grid.compute(a)


def test_final_relocated_and_explicit_units(tmp_path):
    data=build(tmp_path/'input');publish_final_inputs(data)
    for variable in io.VARS['solar']:
        mp=data['root']/'manifests'/f'M__s__{variable}__P.json';m=io.read_json(mp)
        path=Path(m['final_tasks'][0]['output'])
        m['final_tasks'][0]['output']=str(path).replace(str(data['root']),'/old/production_v1');mp.write_text(json.dumps(m))
        side=Path(str(path)+'.json');s=io.read_json(side);s['output']=m['final_tasks'][0]['output'];side.write_text(json.dumps(s))
        with nc.Dataset(path,'a') as ds:ds[variable].delncattr('units')
    a=arguments(data,tmp_path/'out',input_mode='final',input_units=['uas=m/s','vas=m/s','tas=K','rsds=W/m2'])
    _,got=outputs(grid.compute(a))
    np.testing.assert_allclose(got.solar_cf,reference(data,'solar'),atol=1e-6,rtol=1e-6)


def test_corrupt_existing_blocks_do_not_fallback(tmp_path):
    data=build(tmp_path/'input');publish_final_inputs(data)
    p=data['root']/'blocks/M/s/uas/P/block_uas_0.json';p.write_text('{broken')
    with pytest.raises(ValueError):grid.compute(arguments(data,tmp_path/'out'))
