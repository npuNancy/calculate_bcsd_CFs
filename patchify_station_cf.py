#!/usr/bin/env python3
"""Compute station-only capacity factors from global_bcsd final patch files.

This entry point is the sole production entry and
never writes a grid CF.  It gathers the required meteorological variables to
stations first, then applies the existing nonlinear CF kernels.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import timedelta
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import xarray as xr

from cf_physics import compute_solar_cf_chunk, compute_wind_cf_chunk, _get_power_curve_arrays, power_law_ratio


VARS = {"wind": ("uas", "vas"), "solar": ("rsds", "tas", "uas", "vas")}


class BlockInputsMissing(FileNotFoundError):
    """A block input is absent; auto mode may use the validated final files."""


def _safe(value: str, name: str) -> str:
    if not value or any(c not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-" for c in value):
        raise ValueError(f"{name} contains unsafe characters: {value!r}")
    return value


def find_final(root: str | Path, model: str, scenario: str, variable: str, patch: str) -> Path:
    """Resolve exactly one global_bcsd patch file; broad globs are unsupported."""
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


def _atomic_json(path: str | Path, payload: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".partial.{os.getpid()}")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _year(value) -> int:
    year = getattr(value, "year", None)
    return int(year if year is not None else pd.Timestamp(value).year)


def _time_numbers(values) -> np.ndarray:
    raw = np.asarray(values)
    if np.issubdtype(raw.dtype, np.datetime64):
        return raw.astype("datetime64[ns]").astype(np.int64).astype(np.float64) / 3.6e12
    if raw.size == 0:
        return np.array([], dtype=np.float64)
    from netCDF4 import date2num
    calendar = str(getattr(raw.flat[0], "calendar", "standard"))
    return np.asarray(date2num(list(raw), "hours since 0001-01-01 00:00:00", calendar=calendar), dtype=np.float64)


def _validate_time_axis(values, expected_start: int, expected_end: int, path: Path) -> None:
    if len(values) == 0:
        raise ValueError(f"empty block time axis: {path}")
    numbers = _time_numbers(values)
    if np.any(np.diff(numbers) <= 0) or not np.allclose(np.diff(numbers), 3.0, rtol=0.0, atol=1e-8):
        raise ValueError(f"block is not a continuous 3-hour time axis: {path}")
    if _year(values[0]) != expected_start or _year(values[-1]) != expected_end:
        raise ValueError(
            f"block years mismatch for {path}: {_year(values[0])}-{_year(values[-1])}, "
            f"expected {expected_start}-{expected_end}"
        )


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


def _select_patch_stations(stations: pd.DataFrame, patch_manifest: str | Path | None, patch: str) -> pd.DataFrame:
    if not patch_manifest:
        return stations
    manifest = json.loads(Path(patch_manifest).read_text(encoding="utf-8"))
    try:
        w, e, s, n = map(float, manifest["patches"][patch]["core_bbox_360"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"invalid patch manifest entry for {patch}") from exc
    w, e = ((w + 180) % 360) - 180, ((e + 180) % 360) - 180
    lon = stations["lon"].to_numpy(float); lat = stations["lat"].to_numpy(float)
    in_lon = (lon >= w) & (lon < e) if w < e else ((lon >= w) | (lon < e))
    return stations.loc[in_lon & (lat >= s) & ((lat < n) | ((lat == 90) if n == 90 else False))].reset_index(drop=True)


def _block_manifest_path(root: str | Path, model: str, scenario: str, variable: str, patch: str) -> Path:
    root = Path(root).expanduser().resolve()
    return root / "manifests" / f"{model}__{scenario}__{variable}__{patch}.json"


def _block_sidecar(path: Path, variable: str) -> Path:
    prefix = f"{variable}_"
    if not path.stem.startswith(prefix):
        raise ValueError(f"unexpected block filename for {variable}: {path}")
    return path.with_name(f"block_{variable}_{path.stem[len(prefix):]}.json")


def _block_identity(path: Path, variable: str, patch: str, sidecar: Path) -> dict:
    if not path.is_file():
        raise BlockInputsMissing(path)
    if not sidecar.is_file():
        raise BlockInputsMissing(sidecar)
    try:
        meta = json.loads(sidecar.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid block sidecar: {sidecar}") from exc
    if meta.get("kind") != "global-step6-block":
        raise ValueError(f"unexpected block sidecar kind: {sidecar}")
    if str(meta.get("variable")) != variable or patch not in [str(x) for x in meta.get("patch_ids", [])]:
        raise ValueError(f"block sidecar identity mismatch: {sidecar}")
    outputs = meta.get("outputs", {})
    if str(outputs.get(patch, "")) != str(path.resolve()):
        raise ValueError(f"block sidecar output mismatch: {sidecar}")
    stat = path.stat()
    return {"path": str(path.resolve()), "sidecar": str(sidecar.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _read_point_metadata(path: Path, variable: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with xr.open_dataset(path, decode_times=True, chunks=None) as ds:
        if variable not in ds or ds[variable].dims != ("time", "point"):
            raise ValueError(f"block must contain {variable}(time, point): {path}")
        required = ("global_point", "y_index", "x_index", "flat_index", "lat", "lon")
        if any(name not in ds for name in required):
            raise ValueError(f"point metadata incomplete: {path}")
        times = np.asarray(ds["time"].values)
        _validate_time_axis(times, _year(times[0]), _year(times[-1]), path)
        return (
            times,
            np.asarray(ds["global_point"].values, dtype=np.int64),
            np.asarray(ds["y_index"].values, dtype=np.int64),
            np.asarray(ds["x_index"].values, dtype=np.int64),
            np.asarray(ds["flat_index"].values, dtype=np.int64),
        )


def _load_block_plan(args: argparse.Namespace, stations: pd.DataFrame) -> dict:
    if not args.patch_manifest:
        raise ValueError("block mode requires --patch-manifest")
    required_vars = VARS[args.tech]
    manifests: dict[str, dict] = {}
    all_contracts = None
    var_blocks: dict[str, list[dict]] = {}
    for variable in required_vars:
        manifest_path = _block_manifest_path(args.bcsd_root, args.model, args.scenario, variable, args.patch)
        if not manifest_path.is_file():
            raise BlockInputsMissing(manifest_path)
        try:
            meta = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid BCSD block manifest: {manifest_path}") from exc
        if meta.get("kind") != "global-bcsd-three-stage":
            raise ValueError(f"unexpected BCSD manifest kind: {manifest_path}")
        contracts = meta.get("time_block_contract", {}).get("blocks")
        final_tasks = [x for x in meta.get("final_tasks", []) if x.get("patch_id") == args.patch and x.get("variable") == variable]
        if not isinstance(contracts, list) or len(final_tasks) != 1:
            raise ValueError(f"manifest lacks block contract/final task: {manifest_path}")
        block_paths = [Path(x).expanduser().resolve() for x in final_tasks[0].get("block_files", [])]
        if len(block_paths) != len(contracts):
            raise ValueError(f"manifest block count mismatch: {manifest_path}")
        if all_contracts is None:
            all_contracts = contracts
        elif contracts != all_contracts:
            raise ValueError(f"variables have different year-block contracts: {manifest_path}")
        blocks = []
        for index, (contract, path) in enumerate(zip(contracts, block_paths)):
            start_year, end_year = int(contract["start_year"]), int(contract["end_year"])
            side = _block_sidecar(path, variable)
            identity = _block_identity(path, variable, args.patch, side)
            times, global_point, y_index, x_index, flat_index = _read_point_metadata(path, variable)
            _validate_time_axis(times, start_year, end_year, path)
            blocks.append({"index": index, "start_year": start_year, "end_year": end_year, "path": str(path), "identity": identity,
                           "times": times, "global_point": global_point, "y_index": y_index, "x_index": x_index, "flat_index": flat_index})
        manifests[variable] = {"path": str(manifest_path.resolve()), "sha256": _sha256(manifest_path), "meta": meta}
        var_blocks[variable] = blocks
    if all_contracts is None:
        raise ValueError("empty BCSD block contract")
    y0, y1 = (int(x) for x in args.years.split("-", 1)) if "-" in args.years else (int(args.years), int(args.years))
    selected_indices = [i for i, c in enumerate(all_contracts) if int(c["end_year"]) >= y0 and int(c["start_year"]) <= y1]
    if not selected_indices:
        raise ValueError(f"requested years do not intersect BCSD blocks: {args.years}")
    return {"contracts": all_contracts, "indices": selected_indices, "blocks": var_blocks, "manifests": manifests,
            "land_plan": str(Path(manifests[required_vars[0]]["meta"]["land_plan"]).resolve()) if manifests[required_vars[0]]["meta"].get("land_plan") else str(Path(manifests[required_vars[0]]["meta"]["block_tasks"][0]["land_plan_path"]).resolve())}


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


def _gather_points(values: np.ndarray, positions: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Gather (time, point) data using nearest/bilinear point positions."""
    safe = np.where(positions >= 0, positions, 0)
    picked = values[:, safe].astype(np.float32)
    finite = np.isfinite(picked) & (positions[None, :, :] >= 0)
    ww = weights[None, :, :]
    denom = np.where(finite, ww, 0).sum(axis=2)
    out = np.full(picked.shape[:2], np.nan, np.float32)
    np.divide(np.where(finite, picked * ww, 0).sum(axis=2), denom, out=out, where=denom > 0)
    return out


def _point_map(land_plan: str | Path, patch_manifest: str | Path, patch: str, stations: pd.DataFrame,
               method: str, reference_block: Path) -> dict:
    """Build regular-grid-compatible point positions from a Step6 point block."""
    manifest = json.loads(Path(patch_manifest).read_text(encoding="utf-8"))
    try:
        w, e, s, n = map(float, manifest["patches"][patch]["core_bbox_360"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"invalid patch manifest entry for {patch}") from exc
    with xr.open_dataset(land_plan, decode_times=False, chunks=None) as plan:
        if "lat" not in plan or "lon" not in plan:
            raise ValueError(f"land plan lacks lat/lon coordinates: {land_plan}")
        lat = np.asarray(plan["lat"].values, dtype=float)
        lon = np.asarray(plan["lon"].values, dtype=float) % 360.0
        if "point_y" not in plan or "point_x" not in plan:
            raise ValueError(f"land plan lacks point_y/point_x: {land_plan}")
        point_y = np.asarray(plan["point_y"].values, dtype=np.int64)
        point_x = np.asarray(plan["point_x"].values, dtype=np.int64)
    w360, e360 = float(w) % 360.0, float(e) % 360.0
    xkeep = (lon >= w360) & (lon < e360) if w360 < e360 else ((lon >= w360) | (lon < e360))
    ykeep = (lat >= s) & ((lat < n) | ((lat == 90) if n == 90 else False))
    yidx, xidx = np.flatnonzero(ykeep), np.flatnonzero(xkeep)
    if not yidx.size or not xidx.size:
        raise ValueError(f"patch has no target grid coordinates: {patch}")
    local_i0, local_i1, weights, distance = _match(lat[yidx], lon[xidx], stations, method)
    global_y, global_x = yidx[local_i0], xidx[local_i1]
    with xr.open_dataset(reference_block, decode_times=False, chunks=None) as ds:
        block_y = np.asarray(ds["y_index"].values, dtype=np.int64)
        block_x = np.asarray(ds["x_index"].values, dtype=np.int64)
        block_flat = np.asarray(ds["flat_index"].values, dtype=np.int64)
        if not np.array_equal(block_flat, point_y * len(lon) + point_x):
            # The equality is checked only when block and land-plan point order
            # can be compared; a different flat-index convention is rejected.
            expected = block_y * len(lon) + block_x
            if not np.array_equal(block_flat, expected):
                raise ValueError(f"block flat_index is incompatible with land plan: {reference_block}")
    lookup = {(int(y), int(x)): i for i, (y, x) in enumerate(zip(block_y, block_x))}
    positions = np.full(global_y.shape, -1, dtype=np.int64)
    for k in range(positions.shape[0]):
        for j in range(positions.shape[1]):
            positions[k, j] = lookup.get((int(global_y[k, j]), int(global_x[k, j])), -1)
    return {"positions": positions, "weights": weights.astype(np.float32), "distance": distance.astype(float),
            "reference_y": block_y, "reference_x": block_x, "reference_flat": block_flat,
            "point_count": int(len(block_flat))}


def _validate_block_point_order(ds: xr.Dataset, mapping: dict, path: Path) -> None:
    for name, expected in (("y_index", mapping["reference_y"]), ("x_index", mapping["reference_x"]), ("flat_index", mapping["reference_flat"])):
        if name not in ds or not np.array_equal(np.asarray(ds[name].values, dtype=np.int64), expected):
            raise ValueError(f"point order mismatch in block: {path}")


def _read_block_chunk(entries: list[tuple[xr.Dataset, np.ndarray]], variable: str, positions: np.ndarray,
                      weights: np.ndarray, target_numbers: np.ndarray) -> np.ndarray:
    """Read only needed point columns and interpolate one time chunk."""
    valid = positions[positions >= 0]
    unique = np.unique(valid)
    if unique.size == 0:
        return np.full((len(target_numbers), positions.shape[0]), np.nan, np.float32)
    rank = np.full(int(max(unique.max(), 0)) + 1, -1, dtype=np.int64)
    rank[unique] = np.arange(unique.size, dtype=np.int64)
    local = np.where(positions >= 0, rank[np.where(positions >= 0, positions, 0)], -1)
    source_times: list[np.ndarray] = []
    source_values: list[np.ndarray] = []
    lo, hi = float(target_numbers[0]) - 6.0, float(target_numbers[-1]) + 6.0
    for ds, numbers in entries:
        selected = np.flatnonzero((numbers >= lo) & (numbers <= hi))
        if selected.size == 0:
            continue
        values = np.asarray(ds[variable].isel(time=slice(int(selected[0]), int(selected[-1]) + 1), point=unique).values, dtype=np.float32)
        source_times.append(numbers[selected[0]:selected[-1] + 1])
        source_values.append(values)
    if not source_times:
        raise ValueError(f"no source time surrounds target block for {variable}")
    numbers = np.concatenate(source_times)
    values = np.concatenate(source_values, axis=0)
    order = np.argsort(numbers, kind="stable")
    numbers, values = numbers[order], values[order]
    keep = np.concatenate(([True], np.diff(numbers) > 0))
    numbers, values = numbers[keep], values[keep]
    interpolated = xr.DataArray(values, dims=("time_num", "point"), coords={"time_num": numbers}).interp(time_num=target_numbers).values
    return _gather_points(np.asarray(interpolated, dtype=np.float32), local, weights)


def _block_worker(payload: dict) -> dict:
    """Process one BCSD year block in a spawned child process."""
    args = argparse.Namespace(**payload["args"])
    stations = pd.DataFrame(payload["stations"])
    mapping = {k: np.asarray(v) for k, v in payload["mapping"].items()}
    block = payload["block"]
    required = VARS[args.tech]
    paths_by_var = {v: {k: (Path(x) if x else None) for k, x in payload["paths"][v].items()} for v in required}
    opened: dict[str, list[tuple[xr.Dataset, np.ndarray]]] = {}
    output = Path(payload["output"])
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not args.overwrite:
        side = Path(str(output) + ".json")
        if side.is_file() and json.loads(side.read_text(encoding="utf-8")).get("status") == "COMPLETED":
            return {"output": str(output), "sidecar": str(side), "skipped": True}
        raise ValueError(f"existing block output lacks a valid sidecar: {output}")
    try:
        for variable in required:
            entries = []
            for path in (paths_by_var[variable][k] for k in ("previous", "current", "next")):
                if path is None:
                    continue
                ds = xr.open_dataset(path, decode_times=True, chunks=None)
                _validate_block_point_order(ds, mapping, path)
                times = np.asarray(ds["time"].values)
                entries.append((ds, _time_numbers(times)))
            opened[variable] = entries
        ref_var = required[0]
        ref_ds = opened[ref_var][1 if paths_by_var[ref_var]["previous"] is not None else 0][0]
        ref_times = np.asarray(ref_ds["time"].values)
        ref_numbers = _time_numbers(ref_times)
        years = np.asarray([_year(t) for t in ref_times])
        y0, y1 = (int(x) for x in args.years.split("-", 1)) if "-" in args.years else (int(args.years), int(args.years))
        selected = np.flatnonzero((years >= max(y0, int(block["start_year"]))) & (years <= min(y1, int(block["end_year"]))))
        if selected.size == 0:
            raise ValueError(f"block has no requested years: {block}")
        times = ref_times[selected]
        target_numbers = ref_numbers[selected]
        nstations = len(stations)
        weights = np.asarray(mapping["weights"], dtype=np.float32)
        bad = np.asarray(mapping["distance"], dtype=float) > args.max_distance_deg
        lats, lons = stations.lat.to_numpy(float), stations.lon.to_numpy(float)
        name = "wind_cf" if args.tech == "wind" else "solar_cf"
        args._bcsd_files = {v: str(paths_by_var[v]["current"]) for v in required}
        args._input_mode = "blocks"
        tmp = output.with_suffix(output.suffix + f".partial.{os.getpid()}")
        nt = len(times)
        tb = max(32, min(512, int(2.5e8 // max(1, nstations * max(1, weights.shape[1]) * 4 * 6))))
        with _StreamingStationWriter(tmp, times, stations, name, args) as writer:
            for start in range(0, nt, tb):
                stop = min(nt, start + tb)
                values: dict[str, np.ndarray] = {}
                for variable in required:
                    values[variable] = _read_block_chunk(opened[variable], variable, mapping["positions"], weights, target_numbers[start:stop])
                    values[variable][:, bad] = np.nan
                    units = str(opened[variable][0][0][variable].attrs.get("units", "")).lower().replace(" ", "")
                    if variable == "rsds" and "kw" not in units:
                        values[variable] /= 1000.0
                    elif variable == "tas" and (units in {"k", "kelvin"} or (not units and np.nanmedian(values[variable]) > 100)):
                        values[variable] -= 273.15
                if args.tech == "wind":
                    cf = _wind_points(values["uas"], values["vas"])
                else:
                    doy, hour = _doy_hour(times[start:stop])
                    cf = _solar_points_block(values["rsds"], values["tas"], values["uas"], values["vas"], doy, hour, lats, lons)
                writer.write(np.clip(cf, 0, 1).astype(np.float32), start)
        os.replace(tmp, output)
        side = Path(str(output) + ".json")
        _atomic_json(side, {"status": "COMPLETED", "input_mode": "blocks", "block_index": int(block["index"]),
                            "start_year": int(block["start_year"]), "end_year": int(block["end_year"]),
                            "model": args.model, "scenario": args.scenario, "patch_id": args.patch, "tech": args.tech,
                            "station_count": len(stations), "output": str(output.resolve())})
        return {"output": str(output), "sidecar": str(side), "block_index": int(block["index"]), "skipped": False}
    finally:
        for entries in opened.values():
            for ds, _ in entries:
                ds.close()


def _compute_blocks(args: argparse.Namespace, stations: pd.DataFrame) -> Path:
    plan = _load_block_plan(args, stations)
    selected = plan["indices"]
    ref_path = Path(plan["blocks"][VARS[args.tech][0]][selected[0]]["path"])
    mapping = _point_map(plan["land_plan"], args.patch_manifest, args.patch, stations,
                         args.spatial_method, ref_path)
    parts_root = Path(args.parts_root or (Path(args.output_root) / ".cf_parts")) / args.model / args.scenario / args.patch / args.tech
    parts_root.mkdir(parents=True, exist_ok=True)
    payloads = []
    for index in selected:
        paths = {}
        for variable in VARS[args.tech]:
            blocks = plan["blocks"][variable]
            paths[variable] = {
                "previous": blocks[index - 1]["path"] if index > 0 else None,
                "current": blocks[index]["path"],
                "next": blocks[index + 1]["path"] if index + 1 < len(blocks) else None,
            }
        payloads.append({"args": {k: v for k, v in vars(args).items() if not k.startswith("_")},
                         "stations": stations.to_dict(orient="list"),
                         "mapping": {k: np.asarray(v).tolist() for k, v in mapping.items() if k in {"positions", "weights", "distance", "reference_y", "reference_x", "reference_flat"}},
                         "block": plan["blocks"][VARS[args.tech][0]][index],
                         "paths": paths,
                         "output": str(parts_root / f"block_{index:02d}.nc")})
    results = []
    process_count = max(1, int(args.processes))
    if process_count == 1 or len(payloads) == 1:
        results = [_block_worker(p) for p in payloads]
    else:
        ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=min(process_count, len(payloads)), mp_context=ctx) as pool:
            futures = [pool.submit(_block_worker, p) for p in payloads]
            for future in as_completed(futures):
                results.append(future.result())
    results.sort(key=lambda x: json.loads(Path(x["sidecar"]).read_text(encoding="utf-8")).get("block_index", 0))
    out = Path(args.output_root) / args.model / args.scenario / args.patch / f"{args.tech}.nc"
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and not args.overwrite:
        return out
    first = xr.open_dataset(results[0]["output"], decode_times=True)
    try:
        times = first["time"].values
        station_ids = np.asarray(first["station_id"].values).astype(str)
        for result in results[1:]:
            with xr.open_dataset(result["output"], decode_times=True) as ds:
                if not np.array_equal(np.asarray(ds["station_id"].values).astype(str), station_ids):
                    raise ValueError("station identity differs between CF block outputs")
                times = np.concatenate((times, ds["time"].values))
        numbers = _time_numbers(times)
        if np.any(np.diff(numbers) <= 0):
            raise ValueError("CF block outputs have overlapping or non-monotonic times")
        tmp = out.with_suffix(out.suffix + f".partial.{os.getpid()}")
        args._bcsd_files = {v: plan["blocks"][v][selected[0]]["path"] for v in VARS[args.tech]}
        args._input_mode = "blocks"
        with _StreamingStationWriter(tmp, times, stations, "wind_cf" if args.tech == "wind" else "solar_cf", args) as writer:
            offset = 0
            for result in results:
                with xr.open_dataset(result["output"], decode_times=True) as ds:
                    values = ds[writer.var_name]
                    for start in range(0, values.sizes["time"], 512):
                        chunk = np.asarray(values.isel(time=slice(start, start + 512)).values, dtype=np.float32)
                        writer.write(chunk, offset + start)
                    offset += values.sizes["time"]
        os.replace(tmp, out)
    finally:
        first.close()
    side = Path(str(out) + ".json")
    _atomic_json(side, {"status": "COMPLETED", "input_mode": "blocks", "model": args.model,
                        "scenario": args.scenario, "patch_id": args.patch, "tech": args.tech,
                        "station_count": len(stations), "processes": process_count,
                        "block_indices": selected, "parts": [x["output"] for x in results],
                        "manifest_sha256": {v: plan["manifests"][v]["sha256"] for v in VARS[args.tech]},
                        "output": str(out.resolve())})
    return out


def _solar_points(rsds, tas, uas, vas, times, lat, lon):
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
        out[:, k] = compute_solar_cf_chunk(rsds[:, k:k+1, None], tas[:, k:k+1, None], uas[:, k:k+1, None], vas[:, k:k+1, None], doy, hour, np.array([lat[k]]), np.array([lon[k]])).reshape(-1)
    return out


def _wind_points(uas, vas):
    ws, pw, rated = _get_power_curve_arrays()
    return compute_wind_cf_chunk(uas, vas, ws, pw, rated, power_law_ratio()).astype(np.float32)


def _doy_hour(times) -> tuple[np.ndarray, np.ndarray]:
    def _stamp(t):
        return t if hasattr(t, "hour") else pd.Timestamp(t)
    def _doy(t):
        value = getattr(t, "dayofyr", None)
        if value is None:
            value = getattr(t, "dayofyear", None)
        return value if value is not None else pd.Timestamp(t).dayofyear
    doy = np.asarray([_doy(t) for t in times], np.float32)
    hour = np.asarray([_stamp(t).hour + _stamp(t).minute / 60 + _stamp(t).second / 3600
                       for t in times], np.float32)
    return doy, hour


def _solar_points_block(rsds, tas, uas, vas, doy, hour, lats, lons):
    # The solar kernels index stations as a leading batch axis crossed with
    # lon (cos_sza would broadcast (T,B,B)), so keep the legacy per-station
    # call shape (T,1,None); values are identical to the original loop.
    out = np.empty(rsds.shape, dtype=np.float32)
    for k in range(rsds.shape[1]):
        got = compute_solar_cf_chunk(
            rsds[:, k:k+1, None], tas[:, k:k+1, None], uas[:, k:k+1, None], vas[:, k:k+1, None],
            doy, hour, np.array([lats[k]]), np.array([lons[k]]))
        out[:, k] = np.asarray(got).reshape(-1)
    return out


class _StreamingStationWriter:
    """netCDF4 writer that appends station CF time blocks without holding
    the full (time, station) array in memory."""

    def __init__(self, path, times, stations, var_name, args):
        import netCDF4
        self.path = Path(path)
        self.var_name = var_name
        self.ds = netCDF4.Dataset(self.path, "w", format="NETCDF4")
        self.ds.createDimension("time", len(times))
        self.ds.createDimension("station", len(stations))
        t = self.ds.createVariable("time", "f8", ("time",))
        from xarray.coding.times import encode_cf_datetime
        encoded, units, calendar = encode_cf_datetime(
            times.astype("datetime64[ns]") if np.issubdtype(np.asarray(times).dtype, np.datetime64) else list(times),
            "hours since 1970-01-01")
        t.units = units; t.calendar = calendar; t.standard_name = "time"
        t[:] = encoded
        s = self.ds.createVariable("station", "i4", ("station",)); s[:] = np.arange(len(stations), dtype=np.int32)
        s.long_name = "station index"
        for name, values, dtype in (
            ("station_id", np.asarray(stations.station_id.tolist(), dtype=object), str),
            ("lon", stations.lon.to_numpy(np.float64), "f8"),
            ("lat", stations.lat.to_numpy(np.float64), "f8"),
            ("capacity_mw", stations.capacity_mw.to_numpy(np.float64).astype(np.float64), "f8"),
        ):
            v = self.ds.createVariable(name, dtype if dtype is str else dtype, ("station",))
            if dtype is str:
                v[:] = np.asarray(values, dtype=object)
            else:
                v[:] = values
        cf = self.ds.createVariable(var_name, "f4", ("time", "station"),
                                    zlib=True, complevel=args.compress_level)
        cf.long_name = f"{args.tech} capacity factor"
        self.cf = cf
        self.ds.setncatts({
            "source": "global_bcsd_patch", "model": args.model, "scenario": args.scenario,
            "patch_id": args.patch, "tech": args.tech, "spatial_method": args.spatial_method,
            "bcsd_files": json.dumps({k: str(v) for k, v in args._bcsd_files.items()}),
            "station_only": "true",
            "input_mode": getattr(args, "_input_mode", "final"),
        })

    def write(self, block, t_start):
        self.cf[t_start:t_start + block.shape[0], :] = block

    def __enter__(self): return self

    def __exit__(self, *exc):
        self.ds.close()
        return False


def compute(args: argparse.Namespace) -> Path:
    stations = load_stations(args.stations_csv, args.tech, args.scenario)
    stations = _select_patch_stations(stations, args.patch_manifest, args.patch)
    mode = getattr(args, "input_mode", "final")
    if mode == "block":
        mode = "blocks"
    if mode not in {"final", "blocks", "auto"}:
        raise ValueError(f"invalid input mode: {mode}")
    if mode in {"blocks", "auto"} and not stations.empty:
        try:
            return _compute_blocks(args, stations)
        except BlockInputsMissing:
            if mode == "blocks":
                raise
            args._input_mode = "final_fallback"
    else:
        args._input_mode = "final"
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
        # Station gather only needs the grid rows/columns the interpolation
        # touches; crop to them before touching data. Dense patches (tens of
        # thousands of stations spread over the whole grid) still need far
        # more memory than a node offers if materialized at once, so the
        # remaining pipeline streams in time blocks written straight to disk.
        used_lat = np.unique(i0); used_lon = np.unique(i1)
        lat_pos = np.full(i0.max() + 1, -1, dtype=np.int64); lat_pos[used_lat] = np.arange(len(used_lat))
        lon_pos = np.full(i1.max() + 1, -1, dtype=np.int64); lon_pos[used_lon] = np.arange(len(used_lon))
        i0 = lat_pos[i0]; i1 = lon_pos[i1]
        # Keep variables lazy here: interpolating off-axis time to the
        # reference axis is eager and would materialize the full cropped grid
        # (dense patches = the whole patch). Interpolation moves into the time
        # block loop below, padded by one source step so block edges keep
        # their linear neighbours.
        cropped = {}
        off_axis = {}
        for v, ds in opened.items():
            da = _var(ds, v)
            dt = _coord(ds, ("time", "valid_time"))
            if dt != tn or not np.array_equal(ds[dt].values, times):
                if v == "rsds":
                    raise ValueError(f"{v} must provide the reference time axis")
                off_axis[v] = dt
            cropped[v] = (da.isel({ln: used_lat, on: used_lon}), dt)
        units = {v: str(da.attrs.get("units", "")).lower().replace(" ", "") for v, (da, _dt) in cropped.items()}
        bad = distance > args.max_distance_deg
        doy, hour = _doy_hour(times)
        lats = stations.lat.to_numpy(float); lons = stations.lon.to_numpy(float)
        if args.tech == "wind":
            name = "wind_cf"
        else:
            name = "solar_cf"
        nt, ns = len(times), len(stations)
        tb = max(64, min(1024, int(1.5e9 // max(1, ns * 4 * 4))))
        out = Path(args.output_root) / args.model / args.scenario / args.patch / f"{args.tech}.nc"
        out.parent.mkdir(parents=True, exist_ok=True)
        if out.exists() and not args.overwrite: return out
        tmp = out.with_suffix(out.suffix + f".partial.{os.getpid()}")
        args._bcsd_files = files
        with _StreamingStationWriter(tmp, times, stations, name, args) as writer:
            for t_start in range(0, nt, tb):
                t_stop = min(nt, t_start + tb)
                blk = {}
                for v, (da, dt) in cropped.items():
                    if v in off_axis:
                        src_times = da[dt].values
                        pad_lo = 1 if t_start > 0 else 0
                        pad_hi = 1 if t_stop < len(src_times) else 0
                        sub = da.isel({dt: slice(t_start - pad_lo, t_stop + pad_hi)})
                        sub = sub.interp({dt: times[t_start:t_stop]}, method="linear")
                    else:
                        sub = da.isel({dt: slice(t_start, t_stop)})
                    arr = _gather(np.asarray(sub.transpose(dt, ln, on).values), i0, i1, weights)
                    arr[:, bad] = np.nan
                    u = units[v]
                    if v == "rsds" and "kw" not in u:
                        arr = arr / 1000.0
                    elif v == "tas" and (u in {"k", "kelvin"} or (not u and np.nanmedian(arr) > 100)):
                        arr = arr - 273.15
                    blk[v] = arr
                if args.tech == "wind":
                    cf = _wind_points(blk["uas"], blk["vas"])
                else:
                    cf = _solar_points_block(blk["rsds"], blk["tas"], blk["uas"], blk["vas"], doy[t_start:t_stop], hour[t_start:t_stop], lats, lons)
                writer.write(np.clip(cf, 0, 1).astype(np.float32), t_start)
        os.replace(tmp, out)
        side = Path(str(out) + ".json"); side.write_text(json.dumps({"status": "COMPLETED", "input_mode": getattr(args, "_input_mode", "final"), "model": args.model, "scenario": args.scenario, "patch_id": args.patch, "tech": args.tech, "station_count": len(stations), "output": str(out)}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        return out
    finally:
        for ds in opened.values(): ds.close()


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="global_bcsd patch -> station-only CF")
    p.add_argument("--bcsd-root", "--bcsd_root", required=True); p.add_argument("--model", required=True); p.add_argument("--scenario", required=True); p.add_argument("--patch", required=True); p.add_argument("--patch-manifest", default=None); p.add_argument("--stations-csv", "--stations_csv", required=True); p.add_argument("--tech", choices=("wind", "solar"), required=True); p.add_argument("--years", default="2015-2060"); p.add_argument("--output-root", "--output_root", required=True); p.add_argument("--spatial-method", choices=("nearest", "bilinear"), default="nearest"); p.add_argument("--max-distance-deg", type=float, default=0.15); p.add_argument("--compress-level", type=int, default=4); p.add_argument("--overwrite", action="store_true"); p.add_argument("--input-mode", choices=("final", "block", "blocks", "auto"), default="final"); p.add_argument("--processes", type=int, default=8); p.add_argument("--parts-root", default=None); return p


if __name__ == "__main__":
    compute(parser().parse_args())
