#!/usr/bin/env python3
"""Compute station-only capacity factors from global_bcsd final patch files.

This entry point intentionally does not call the historical region pipeline and
never writes a grid CF.  It gathers the required meteorological variables to
stations first, then applies the existing nonlinear CF kernels.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import xarray as xr


VARS = {"wind": ("uas", "vas"), "solar": ("rsds", "tas", "uas", "vas")}


def _safe(value: str, name: str) -> str:
    if not value or any(c not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-" for c in value):
        raise ValueError(f"{name} contains unsafe characters: {value!r}")
    return value


def find_final(root: str | Path, model: str, scenario: str, variable: str, patch: str) -> Path:
    """Resolve exactly one global final file; region globs are deliberately unsupported."""
    root = Path(root).expanduser().resolve()
    _safe(model, "model"); _safe(scenario, "scenario"); _safe(variable, "variable"); _safe(patch, "patch")
    direct = [
        root / "outputs" / model / scenario / variable / f"{patch}.nc",
        root / "outputs" / model / scenario / variable / f"{variable}_{model}_{scenario}_{patch}.nc",
        root / model / scenario / variable / f"{patch}.nc",
    ]
    matches = [p for p in direct if p.is_file()]
    if not matches:
        base = root / "outputs" / model / scenario / variable
        matches = sorted(base.glob(f"*{patch}*.nc")) if base.is_dir() else []
    if len(matches) != 1:
        raise FileNotFoundError(f"expected exactly one final file for {model}/{scenario}/{variable}/{patch}, got {matches}")
    sidecar = Path(str(matches[0]) + ".json")
    if not sidecar.is_file():
        raise FileNotFoundError(f"missing BCSD sidecar: {sidecar}")
    meta = json.loads(sidecar.read_text(encoding="utf-8"))
    for key, expected in (("patch_id", patch), ("variable", variable)):
        if str(meta.get(key, "")) != expected:
            raise ValueError(f"sidecar {sidecar} has {key}={meta.get(key)!r}, expected {expected!r}")
    return matches[0]


def _var(ds: xr.Dataset, preferred: str) -> xr.DataArray:
    if preferred in ds.data_vars:
        return ds[preferred]
    candidate = f"{preferred}_bcsd"
    if candidate in ds.data_vars:
        return ds[candidate]
    if len(ds.data_vars) == 1:
        return ds[next(iter(ds.data_vars))]
    raise KeyError(f"cannot resolve {preferred}; variables={list(ds.data_vars)}")


def _coord(ds: xr.Dataset, names: Iterable[str]) -> str:
    for name in names:
        if name in ds.coords or name in ds.dims:
            return name
    raise KeyError(f"missing coordinate; tried {tuple(names)}")


def load_stations(path: str | Path, tech: str, scenario: str = "ssp126") -> pd.DataFrame:
    df = pd.read_csv(path)
    def col(*names: str) -> str:
        for name in names:
            if name in df:
                return name
        raise ValueError(f"station CSV missing one of {names}")
    lon, lat = col("lon", "longitude", "x"), col("lat", "latitude", "y")
    type_col = next((x for x in ("tech", "type", "technology") if x in df), None)
    if type_col:
        vals = df[type_col].astype(str).str.lower()
        df = df.loc[vals.isin({tech.lower(), "wind" if tech == "wind" else "solar"})].copy()
    out = pd.DataFrame({"lon": pd.to_numeric(df[lon]), "lat": pd.to_numeric(df[lat])})
    if out.empty:
        return out
    out["lon"] = ((out["lon"] + 180.0) % 360.0) - 180.0
    cap_col = next((x for x in ("capacity_mw", "capacity_gw", "capacity") if x in df), None)
    if cap_col:
        out["capacity_mw"] = pd.to_numeric(df[cap_col], errors="coerce").fillna(0).to_numpy(float) * (1000.0 if cap_col == "capacity_gw" else 1.0)
    else:
        out["capacity_mw"] = 1.0
    out = out.dropna().drop_duplicates(["lon", "lat"]).reset_index(drop=True)
    out["station_id"] = [hashlib.sha1(f"{scenario}|{tech}|{float(x):.4f}|{float(y):.4f}".encode()).hexdigest()[:20] for x, y in zip(out.lon, out.lat)]
    return out


def _match(lat: np.ndarray, lon: np.ndarray, stations: pd.DataFrame, method: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if method not in {"nearest", "bilinear"}:
        raise ValueError("spatial method must be nearest or bilinear")
    lat = np.asarray(lat, dtype=float); lon = ((np.asarray(lon, dtype=float) + 180) % 360) - 180
    slat = stations.lat.to_numpy(float); slon = stations.lon.to_numpy(float)
    ilat = np.abs(lat[:, None] - slat[None, :]).argmin(axis=0)
    dlon = np.abs(((lon[:, None] - slon[None, :] + 180) % 360) - 180)
    ilon = dlon.argmin(axis=0)
    dist = np.hypot(lat[ilat] - slat, dlon[ilon, np.arange(len(stations))])
    if method == "nearest":
        return ilat[:, None], ilon[:, None], np.ones((len(stations), 1), np.float32), dist
    # Four-point regular-grid interpolation weights.  Out-of-domain points are
    # clamped and retain the nearest-grid distance check.
    ys = np.argsort(lat); xs = np.argsort(lon); la = lat[ys]; lo = lon[xs]
    i0 = np.empty((len(stations), 4), np.int64); i1 = np.empty_like(i0); w = np.empty_like(i0, dtype=np.float32)
    for k, (y, x) in enumerate(zip(slat, slon)):
        py = int(np.searchsorted(la, y, side="right")); py = max(1, min(py, len(la)-1))
        px = int(np.searchsorted(lo, x, side="right")) % len(lo)
        a, b = py-1, py; c, d = (px-1) % len(lo), px % len(lo)
        wy = 0.0 if la[b] == la[a] else float(np.clip((y-la[a])/(la[b]-la[a]), 0, 1))
        dx = (lo[d]-lo[c]) % 360.0; wx = 0.0 if dx == 0 else float(np.clip(((x-lo[c]) % 360.0)/dx, 0, 1))
        i0[k] = [ys[a], ys[a], ys[b], ys[b]]; i1[k] = [xs[c], xs[d], xs[c], xs[d]]
        w[k] = [(1-wy)*(1-wx), (1-wy)*wx, wy*(1-wx), wy*wx]
    return i0, i1, w, dist


def _gather(values: np.ndarray, i0: np.ndarray, i1: np.ndarray, weights: np.ndarray) -> np.ndarray:
    picked = values[:, i0, i1].astype(np.float32)
    finite = np.isfinite(picked)
    ww = weights[None, :, :]
    denom = np.where(finite, ww, 0).sum(axis=2)
    out = np.full(picked.shape[:2], np.nan, np.float32)
    np.divide(np.where(finite, picked * ww, 0).sum(axis=2), denom, out=out, where=denom > 0)
    return out


def _load_kernel(name: str):
    path = Path(__file__).with_name("S01E01_Simulate_Solar_CF_BCSD.py" if name == "solar" else "S01E02_Simulate_Wind_CF_BCSD.py")
    spec = importlib.util.spec_from_file_location(f"_cf_{name}_kernel", path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def _solar_points(rsds, tas, uas, vas, times, lat, lon):
    mod = _load_kernel("solar")
    out = np.empty_like(rsds, dtype=np.float32)
    def _doy(t):
        value = getattr(t, "dayofyr", None)
        if value is None:
            value = getattr(t, "dayofyear", None)
        return value if value is not None else pd.Timestamp(t).dayofyear
    doy = np.asarray([_doy(t) for t in times], np.float32)
    hour = np.asarray([t.hour + getattr(t, "minute", 0)/60 + getattr(t, "second", 0)/3600 for t in times], np.float32)
    # The legacy kernel is grid-shaped; one-point calls preserve its formula
    # exactly while keeping the persisted output station-only.
    for k in range(rsds.shape[1]):
        out[:, k] = mod.compute_solar_cf_chunk(rsds[:, k:k+1, None], tas[:, k:k+1, None], uas[:, k:k+1, None], vas[:, k:k+1, None], doy, hour, np.array([lat[k]]), np.array([lon[k]])).reshape(-1)
    return out


def _wind_points(uas, vas):
    mod = _load_kernel("wind")
    ws, pw, rated = mod._get_power_curve_arrays()
    # _get_power_curve_arrays returns (ws, power, rated); existing kernel needs ratio.
    ratio = mod.power_law_ratio(100.0, 10.0, 1/7)
    return mod.compute_wind_cf_chunk(uas, vas, ws, pw, rated, ratio).astype(np.float32)


def compute(args: argparse.Namespace) -> Path:
    stations = load_stations(args.stations_csv, args.tech, args.scenario)
    if args.patch_manifest:
        manifest = json.loads(Path(args.patch_manifest).read_text(encoding="utf-8"))
        try:
            w, e, s, n = map(float, manifest["patches"][args.patch]["core_bbox_360"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid patch manifest entry for {args.patch}") from exc
        w, e = ((w + 180) % 360) - 180, ((e + 180) % 360) - 180
        lon = stations["lon"].to_numpy(float); lat = stations["lat"].to_numpy(float)
        in_lon = (lon >= w) & (lon < e) if w < e else ((lon >= w) | (lon < e))
        stations = stations.loc[in_lon & (lat >= s) & ((lat < n) | ((lat == 90) if n == 90 else False))].reset_index(drop=True)
    files = {v: find_final(args.bcsd_root, args.model, args.scenario, v, args.patch) for v in VARS[args.tech]}
    if stations.empty:
        out = Path(args.output_root) / args.model / args.scenario / args.patch / f"{args.tech}.nc"
        out.parent.mkdir(parents=True, exist_ok=True)
        empty = xr.Dataset(coords={"time": np.array([], dtype="datetime64[ns]"), "station": np.array([], dtype=np.int32)})
        empty.attrs.update(source="global_bcsd_patch", model=args.model, scenario=args.scenario, patch_id=args.patch, tech=args.tech, station_only="true", status="SKIPPED_NO_STATIONS")
        empty.to_netcdf(out)
        Path(str(out)+".json").write_text(json.dumps({"status":"SKIPPED_NO_STATIONS","model":args.model,"scenario":args.scenario,"patch_id":args.patch,"tech":args.tech,"output":str(out)}, indent=2)+"\n", encoding="utf-8")
        return out
    opened = {v: xr.open_dataset(p) for v, p in files.items()}
    try:
        ref = opened[VARS[args.tech][0]]; tn = _coord(ref, ("time", "valid_time")); ln = _coord(ref, ("lat", "latitude")); on = _coord(ref, ("lon", "longitude"))
        y0, y1 = (int(x) for x in args.years.split("-", 1)) if "-" in args.years else (int(args.years), int(args.years))
        raw_times = ref[tn].values
        def _year(t):
            value = getattr(t, "year", None)
            return value if value is not None else pd.Timestamp(t).year
        year_values = np.asarray([int(_year(t)) for t in raw_times])
        selected = np.flatnonzero((year_values >= y0) & (year_values <= y1))
        if selected.size == 0: raise ValueError(f"no time points in --years {args.years}")
        ref = ref.isel({tn: selected}); times = ref[tn].values
        i0, i1, weights, distance = _match(ref[ln].values, ref[on].values, stations, args.spatial_method)
        arrays = {}
        for v, ds in opened.items():
            da = _var(ds, v)
            dt = _coord(ds, ("time", "valid_time"))
            if dt != tn or not np.array_equal(ds[dt].values, times):
                if v == "rsds":
                    raise ValueError(f"{v} must provide the reference time axis")
                da = da.interp({dt: ref[tn]}, method="linear")
            arrays[v] = _gather(np.asarray(da.transpose(dt, ln, on).values), i0, i1, weights)
            arrays[v][:, distance > args.max_distance_deg] = np.nan
            units = str(da.attrs.get("units", "")).lower().replace(" ", "")
            if v == "rsds":
                # Existing solar kernel consumes kW m-2; global_bcsd final is
                # normally W m-2, while preserving an explicit kW input.
                if "kw" not in units:
                    arrays[v] = arrays[v] / 1000.0
            elif v == "tas":
                if units in {"k", "kelvin"} or (not units and np.nanmedian(arrays[v]) > 100):
                    arrays[v] = arrays[v] - 273.15
        if args.tech == "wind":
            cf = _wind_points(arrays["uas"], arrays["vas"]); name = "wind_cf"
        else:
            cf = _solar_points(arrays["rsds"], arrays["tas"], arrays["uas"], arrays["vas"], times, stations.lat.values, stations.lon.values); name = "solar_cf"
        dsout = xr.Dataset({name: (("time", "station"), np.clip(cf, 0, 1).astype(np.float32))}, coords={"time": times, "station": np.arange(len(stations), dtype=np.int32), "station_id": ("station", stations.station_id.values), "lon": ("station", stations.lon.values), "lat": ("station", stations.lat.values), "capacity_mw": ("station", stations.capacity_mw.values.astype(np.float64))})
        dsout.attrs.update(source="global_bcsd_patch", model=args.model, scenario=args.scenario, patch_id=args.patch, tech=args.tech, spatial_method=args.spatial_method, bcsd_files=json.dumps({k: str(v) for k,v in files.items()}), station_only="true")
        out = Path(args.output_root) / args.model / args.scenario / args.patch / f"{args.tech}.nc"; out.parent.mkdir(parents=True, exist_ok=True)
        if out.exists() and not args.overwrite: return out
        tmp = out.with_suffix(out.suffix + f".partial.{os.getpid()}"); dsout.to_netcdf(tmp, engine="netcdf4", encoding={name: {"zlib": True, "complevel": args.compress_level, "dtype": "f4"}}); os.replace(tmp, out)
        side = Path(str(out) + ".json"); side.write_text(json.dumps({"model": args.model, "scenario": args.scenario, "patch_id": args.patch, "tech": args.tech, "station_count": len(stations), "output": str(out)}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        return out
    finally:
        for ds in opened.values(): ds.close()


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="global_bcsd patch -> station-only CF")
    p.add_argument("--bcsd-root", "--bcsd_root", required=True); p.add_argument("--model", required=True); p.add_argument("--scenario", required=True); p.add_argument("--patch", required=True); p.add_argument("--patch-manifest", default=None); p.add_argument("--stations-csv", "--stations_csv", required=True); p.add_argument("--tech", choices=("wind", "solar"), required=True); p.add_argument("--years", default="2015-2060"); p.add_argument("--output-root", "--output_root", required=True); p.add_argument("--spatial-method", choices=("nearest", "bilinear"), default="nearest"); p.add_argument("--max-distance-deg", type=float, default=0.15); p.add_argument("--compress-level", type=int, default=4); p.add_argument("--overwrite", action="store_true"); return p


if __name__ == "__main__":
    compute(parser().parse_args())
