"""Regression test: station gather with grid cropping matches full-grid gather.

Synthesizes a tiny BCSD-like patch dataset, runs the cropping logic used in
patchify_station_cf.compute (unique-row/col isel + index remap) against the
original full-grid path for both nearest and bilinear, and asserts bit-equal
station series.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from patchify_station_cf import _match, _gather  # noqa: E402


def synth_grid(nt=24, ny=12, nx=15, seed=7):
    rng = np.random.default_rng(seed)
    times = pd.date_range("2015-01-01", periods=nt, freq="D")
    lat = np.linspace(-50.0, 50.0, ny)
    lon = np.linspace(10.0, 40.0, nx)
    data = rng.normal(size=(nt, ny, nx)).astype(np.float32)
    da = xr.DataArray(data, coords={"time": times, "lat": lat, "lon": lon},
                      dims=("time", "lat", "lon"))
    return da


def stations_frame(pts):
    return pd.DataFrame({"lon": [p[0] for p in pts], "lat": [p[1] for p in pts],
                         "capacity_mw": [1.0] * len(pts)})


def run_case(da, stations, method):
    i0, i1, weights, distance = _match(da.lat.values, da.lon.values, stations, method)
    full = _gather(np.asarray(da.transpose("time", "lat", "lon").values), i0, i1, weights)

    used_lat = np.unique(i0); used_lon = np.unique(i1)
    lat_pos = np.full(i0.max() + 1, -1, dtype=np.int64); lat_pos[used_lat] = np.arange(len(used_lat))
    lon_pos = np.full(i1.max() + 1, -1, dtype=np.int64); lon_pos[used_lon] = np.arange(len(used_lon))
    i0c, i1c = lat_pos[i0], lon_pos[i1]
    small = da.isel(lat=used_lat, lon=used_lon)
    cropped = _gather(np.asarray(small.transpose("time", "lat", "lon").values), i0c, i1c, weights)

    assert np.array_equal(full, cropped), f"gather mismatch for method={method}"
    assert used_lat.size <= da.sizes["lat"] and used_lon.size <= da.sizes["lon"]


def test_nearest_and_bilinear_match():
    da = synth_grid()
    pts = [(12.3, -11.7), (33.9, 41.2), (25.0, 0.05), (10.01, 49.5), (40.0, 33.3)]
    stations = stations_frame(pts)
    run_case(da, stations, "nearest")
    run_case(da, stations, "bilinear")


def test_many_stations_single_row():
    da = synth_grid()
    # many stations squeezed between two grid rows -> used rows tiny
    pts = [(lon, 3.0) for lon in np.linspace(11.0, 39.0, 30)]
    stations = stations_frame(pts)
    run_case(da, stations, "nearest")
    run_case(da, stations, "bilinear")


def test_distance_mask_semantics():
    da = synth_grid()
    pts = [(180.0, 80.0)]  # far outside -> nearest still returns an index
    stations = stations_frame(pts)
    for method in ("nearest", "bilinear"):
        i0, i1, weights, distance = _match(da.lat.values, da.lon.values, stations, method)
        used_lat = np.unique(i0); used_lon = np.unique(i1)
        lat_pos = np.full(i0.max() + 1, -1, dtype=np.int64); lat_pos[used_lat] = np.arange(len(used_lat))
        lon_pos = np.full(i1.max() + 1, -1, dtype=np.int64); lon_pos[used_lon] = np.arange(len(used_lon))
        i0c, i1c = lat_pos[i0], lon_pos[i1]
        full = _gather(np.asarray(da.transpose("time", "lat", "lon").values), i0, i1, weights)
        small = da.isel(lat=used_lat, lon=used_lon)
        cropped = _gather(np.asarray(small.transpose("time", "lat", "lon").values), i0c, i1c, weights)
        assert np.array_equal(full, cropped)


if __name__ == "__main__":
    test_nearest_and_bilinear_match()
    test_many_stations_single_row()
    test_distance_mask_semantics()
    print("all gather-crop regression tests passed")
