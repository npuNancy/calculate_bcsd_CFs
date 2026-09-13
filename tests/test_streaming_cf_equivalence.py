"""End-to-end equivalence test for the time-blocked streaming CF writer.

Builds a tiny synthetic BCSD patch (netCDF files per variable), a stations CSV
and a patch manifest, runs patchify_station_cf.compute() (which now streams in
time blocks), reads the output back and compares against a reference computed
with the original full-materialization formula (gather + physics applied in
one shot).
"""
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from patchify_station_cf import (_gather, _match, _wind_points, _solar_points,
                                 compute, parser)  # noqa: E402


NT, NY, NX = 48, 14, 16


def build(root: Path, interp_off_axis=False):
    rng = np.random.default_rng(3)
    times = pd.date_range("2015-01-01", periods=NT, freq="6h")
    lat = np.linspace(20.0, 23.0, NY)
    lon = np.linspace(100.0, 104.0, NX)
    bcsd = root / "bcsd" / "outputs" / "M" / "s"
    for v, scale in (("uas", 4.0), ("vas", 2.0), ("rsds", 250.0), ("tas", 292.0)):
        d = bcsd / v; d.mkdir(parents=True)
        data = (rng.normal(size=(NT, NY, NX)) * 0.1 + scale).astype(np.float32)
        if v == "rsds":
            data = np.abs(data)
        units = {"rsds": "W m-2", "tas": "K"}.get(v, "m s-1")
        # rsds carries the reference axis offset by 90 minutes, like the real
        # BCSD finals, so tas/uas/vas exercise the block-padded interpolation.
        vt = times if v == "rsds" else times + pd.Timedelta(minutes=90)
        ds = xr.Dataset({f"{v}_bcsd": (("time", "lat", "lon"), data, {"units": units})},
                        coords={"time": vt, "lat": lat, "lon": lon})
        fn = d / f"{v}_M_s_P1.nc"
        ds.to_netcdf(fn)
        fn.with_suffix(".nc.json").write_text(json.dumps(
            {"patch_id": "P1", "variable": v, "model": "M", "scenario": "s"}))
    stations = pd.DataFrame({
        "lon": [101.1, 102.7, 103.9, 100.4],
        "lat": [20.6, 21.8, 22.5, 23.0],
        "capacity_gw": [0.5, 0.7, 0.2, 1.0],
    })
    stations["type"] = ["wind"] * 4
    stations["year"] = [2030] * 4
    stations.to_csv(root / "stations.csv", index=False)
    manifest = {"patches": {"P1": {"core_bbox_360": [100.0, 104.5, 20.0, 23.5]}}}
    (root / "patch_manifest.json").write_text(json.dumps(manifest))
    return lat, lon, times


def reference(lat, lon, times, root, tech):
    import patchify_station_cf as m
    stations = m.load_stations(str(root / "stations.csv"), tech, "s")
    files = {v: root / "bcsd" / "outputs" / "M" / "s" / v / f"{v}_M_s_P1.nc" for v in m.VARS[tech]}
    opened = {v: xr.open_dataset(p) for v, p in files.items()}
    try:
        ref = opened[m.VARS[tech][0]]
        ref_times = ref.time.values
        i0, i1, w, dist = _match(ref.lat.values, ref.lon.values, stations, "nearest")
        arrays = {}
        for v, ds in opened.items():
            da = ds[f"{v}_bcsd"]
            if not np.array_equal(da.time.values, ref_times):
                da = da.interp({"time": ref_times}, method="linear")
            arrays[v] = _gather(np.asarray(da.transpose("time", "lat", "lon").values), i0, i1, w)
        if tech == "wind":
            return _wind_points(arrays["uas"], arrays["vas"])
        doy = np.asarray([pd.Timestamp(t).dayofyear for t in times], np.float32)
        hour = np.asarray([pd.Timestamp(t).hour for t in times], np.float32)
        out = np.empty_like(arrays["rsds"], dtype=np.float32)
        for k in range(arrays["rsds"].shape[1]):
            out[:, k] = m.compute_solar_cf_chunk(
                arrays["rsds"][:, k:k+1, None] / 1000.0,
                arrays["tas"][:, k:k+1, None] - 273.15,
                arrays["uas"][:, k:k+1, None], arrays["vas"][:, k:k+1, None],
                doy, hour, np.array([stations.lat.values[k]]), np.array([stations.lon.values[k]])).reshape(-1)
        return out
    finally:
        for ds in opened.values(): ds.close()


def test_streaming_matches_reference():
    for tech in ("wind", "solar"):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            lat, lon, times = build(root)
            stations_csv = root / "stations.csv"
            pd.DataFrame({"lon": [101.1, 102.7, 103.9, 100.4], "lat": [20.6, 21.8, 22.5, 23.0],
                          "capacity_gw": [0.5, 0.7, 0.2, 1.0], "type": [tech]*4,
                          "year": [2030]*4}).to_csv(stations_csv, index=False)
            out_root = root / "out"
            args = parser().parse_args([
                "--bcsd-root", str(root / "bcsd"), "--model", "M", "--scenario", "s",
                "--patch", "P1", "--patch-manifest", str(root / "patch_manifest.json"),
                "--stations-csv", str(stations_csv), "--tech", tech,
                "--years", "2015-2015", "--output-root", str(out_root), "--overwrite"])
            out = compute(args)
            got = xr.open_dataset(out)
            ref = reference(lat, lon, times, root, tech)
            val = got[f"{tech}_cf"].values
            assert val.shape == ref.shape, (val.shape, ref.shape)
            np.testing.assert_allclose(val, np.clip(ref, 0, 1), rtol=1e-6, atol=1e-6)
            assert got.attrs["model"] == "M" and got.attrs["patch_id"] == "P1"
            assert list(got.dims) == ["time", "station"] or set(got.dims) == {"time", "station"}
            got.close()
            side = json.loads(Path(str(out) + ".json").read_text())
            assert side["station_count"] == 4
    print("streaming CF writer equivalence passed")


if __name__ == "__main__":
    test_streaming_matches_reference()
