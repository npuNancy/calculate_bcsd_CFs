#!/usr/bin/env python3
"""Compute native-grid CF from BCSD blocks or final files in year segments."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import ExitStack
import importlib.metadata
import json
import multiprocessing
import os
from pathlib import Path
import resource
import subprocess
import time

# Spawned workers must not each start another BLAS/OpenMP thread pool.
for _key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import netCDF4 as nc
import numpy as np

import cf_physics as physics
import grid_cf_io as io


def positive(value):
    value = int(value)
    if value <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    for flag in ("bcsd-root", "model", "scenario", "patch", "patch-manifest", "land-plan", "output-root"):
        p.add_argument("--"+flag, required=True)
    p.add_argument("--tech", choices=tuple(io.VARS), required=True)
    p.add_argument("--years", default="2015-2060", help="Complete year range in the BCSD manifest")
    p.add_argument("--processes", type=positive, default=8)
    p.add_argument("--input-mode", choices=("auto", "blocks", "final"), default="auto",
                   help="Prefer complete blocks; otherwise read final files for the whole unit")
    p.add_argument("--time-chunk", type=positive, default=240)
    p.add_argument("--tile-shape", nargs=2, type=positive, default=[64, 64], metavar=("NY", "NX"))
    p.add_argument("--compress-level", type=int, choices=range(10), default=2)
    p.add_argument("--parts-root")
    p.add_argument("--input-units", nargs="+", default=[], metavar="VAR=UNIT",
                   help="Explicit units for blocks missing metadata, e.g. tas=K rsds=W/m2 uas=m/s vas=m/s")
    p.add_argument("--merge-final", action="store_true", default=False)
    p.add_argument("--overwrite", action="store_true")
    return p


def implementation_identity(tech):
    root = Path(__file__).resolve().parent
    files = {name: io.file_identity(root/name, True)["sha256"]
             for name in ("patchify_grid_cf.py", "grid_cf_io.py", "grid_cf_sources.py", "cf_physics.py")}
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, stderr=subprocess.DEVNULL).decode().strip()
    except (OSError, subprocess.CalledProcessError):
        sha = "unavailable"
    result = {"code_sha": sha, "source_sha256": files,
              "libraries": {name: importlib.metadata.version(name) for name in ("numpy", "netCDF4", "cftime", "windpowerlib")},
              "physical_parameters": {k: getattr(physics, k) for k in (
                  "TURBINE_TYPE", "HUB_HEIGHT", "REF_HEIGHT", "POWER_LAW_ALPHA", "CUT_OUT", "SYS_COEF",
                  "SOLAR_CONSTANT_KW", "MIN_COS_ZENITH", "T_STC", "GAMMA_TEMP", "C1", "C2", "C3", "C4")}}
    if tech == "wind":
        ws, power, rated = physics._get_power_curve_arrays()
        result["power_curve_sha256"] = io.array_digest(ws, power, np.asarray(rated))
    return result


def block_axis(plan, tech, index):
    entry = plan["blocks"][io.VARS[tech][0]][index]
    return {k: entry[k] for k in ("raw", "units", "calendar", "numbers")}


def block_identity(identity, index):
    return io.digest({"unit": identity, "block_index": index})


def _time_fields(axis, start, stop):
    dates = nc.num2date(axis["raw"][start:stop], axis["units"], axis["calendar"], only_use_cftime_datetimes=True)
    doy = np.asarray([d.dayofyr for d in dates], dtype=np.float32)
    hour = np.asarray([d.hour+d.minute/60+d.second/3600 for d in dates], dtype=np.float32)
    return doy, hour


def _compute_block(payload):
    args, plan, index, path, identity, provenance = payload
    path = Path(path)
    axis = block_axis(plan, args.tech, index)
    with io.output_lock(str(path)+".lock"):
        if io.completed(path, plan, axis, args.tech, identity) and not args.overwrite:
            return {"index": index, "path": str(path), "reused": True}
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = io.partial_path(path)
        elapsed = {"read_align": 0.0, "physics": 0.0, "write": 0.0}
        begin = time.monotonic()
        try:
            with ExitStack() as stack:
                readers = {v: io.input_reader(plan, v) for v in io.VARS[args.tech]}
                for reader in readers.values():
                    stack.callback(reader.close)
                ds = stack.enter_context(io.create_output(tmp, plan, axis, args, identity))
                ds.provenance = json.dumps(provenance, sort_keys=True)
                curve = physics._get_power_curve_arrays() if args.tech == "wind" else None
                grid_tiles = []
                for bounds in io.tiles(plan, args.tile_shape):
                    y0, y1, x0, x1 = bounds
                    pos = np.flatnonzero((plan["y"] >= y0) & (plan["y"] < y1) & (plan["x"] >= x0) & (plan["x"] < x1))
                    if pos.size:
                        grid_tiles.append((bounds, pos, plan["y"][pos]-y0, plan["x"][pos]-x0))
                for start in range(0, len(axis["raw"]), args.time_chunk):
                    stop = min(start+args.time_chunk, len(axis["raw"]))
                    doy, hour = _time_fields(axis, start, stop) if args.tech == "solar" else (None, None)
                    for (y0, y1, x0, x1), pos, yy, xx in grid_tiles:
                        mark = time.monotonic()
                        values = {}
                        for variable, reader in readers.items():
                            arr = np.full((stop-start, y1-y0, x1-x0), np.nan, dtype=np.float32)
                            arr[:, yy, xx] = reader.read(axis["numbers"][start:stop], pos)
                            values[variable] = arr
                        valid = np.logical_and.reduce([np.isfinite(v) for v in values.values()])
                        elapsed["read_align"] += time.monotonic()-mark
                        mark = time.monotonic()
                        # Invalid cells are restored after the legacy kernels' NaN-to-zero rules.
                        if args.tech == "wind":
                            out = physics.compute_wind_cf_chunk(values["uas"], values["vas"], *curve, physics.power_law_ratio())
                        else:
                            out = physics.compute_solar_cf_chunk(values["rsds"], values["tas"], values["uas"], values["vas"],
                                                                doy, hour, plan["lat"][y0:y1],
                                                                ((plan["lon"][x0:x1].astype(float)+180) % 360)-180)
                        out = np.where(valid & np.isfinite(out), out, io.FILL).astype(np.float32)
                        elapsed["physics"] += time.monotonic()-mark
                        mark = time.monotonic()
                        ds[args.tech+"_cf"][start:stop, y0:y1, x0:x1] = out
                        elapsed["write"] += time.monotonic()-mark
                read_bytes = sum(r.read_bytes for r in readers.values())
            io.unchanged(plan["inputs"])
            elapsed["total"] = time.monotonic()-begin
            metadata = {"block_index": index, **plan["contracts"][index], "time_count": len(axis["raw"]),
                        "provenance": provenance, "timing_seconds": elapsed, "read_array_bytes": read_bytes,
                        "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
            io.publish(tmp, path, plan, axis, args.tech, identity, metadata)
            return {"index": index, "path": str(path), "reused": False, **elapsed}
        finally:
            tmp.unlink(missing_ok=True)


def merge_axis(plan, tech):
    axes = [block_axis(plan, tech, i) for i in range(len(plan["contracts"]))]
    first = axes[0]
    # Different units across blocks are converted to the first reference encoding.
    raws = []
    for a in axes:
        raw = a["raw"]
        if a["units"] != first["units"] or a["calendar"] != first["calendar"]:
            raw = nc.date2num(nc.num2date(raw, a["units"], a["calendar"]), first["units"], first["calendar"])
        raws.append(raw)
    return {"raw": np.concatenate(raws), "units": first["units"], "calendar": first["calendar"],
            "numbers": np.concatenate([a["numbers"] for a in axes])}


def _merge(args, plan, paths, output, identity, provenance):
    axis = merge_axis(plan, args.tech)
    final_identity = io.digest({"unit": identity, "stage": "final"})
    with io.output_lock(str(output)+".lock"):
        if io.completed(output, plan, axis, args.tech, final_identity) and not args.overwrite:
            return {"status": "COMPLETED", "path": str(output), "reused": True}
        tmp = io.partial_path(output)
        begin = time.monotonic()
        try:
            with io.create_output(tmp, plan, axis, args, final_identity) as dest:
                dest.provenance = json.dumps(provenance, sort_keys=True)
                offset = 0
                for index, path in enumerate(paths):
                    current = block_axis(plan, args.tech, index)
                    if not io.completed(path, plan, current, args.tech, block_identity(identity, index)):
                        raise ValueError(f"incomplete CF block during merge: {path}")
                    with nc.Dataset(path) as src:
                        src.set_auto_mask(False)
                        dest.set_auto_mask(False)
                        for start in range(0, len(current["raw"]), args.time_chunk):
                            stop = min(start+args.time_chunk, len(current["raw"]))
                            for y0, y1, x0, x1 in io.tiles(plan, args.tile_shape):
                                if not plan["mask"][y0:y1, x0:x1].any():
                                    continue
                                data = src[args.tech+"_cf"][start:stop, y0:y1, x0:x1]
                                section = (slice(offset+start, offset+stop), slice(y0, y1), slice(x0, x1))
                                dest[args.tech+"_cf"][section] = data
                                if not np.array_equal(dest[args.tech+"_cf"][section], data):
                                    raise ValueError("merged CF differs from block")
                    offset += len(current["raw"])
            io.unchanged(plan["inputs"])
            io.publish(tmp, output, plan, axis, args.tech, final_identity,
                       {"provenance": provenance, "blocks": [io.file_identity(p) for p in paths],
                        "time_count": len(axis["raw"]), "merge_seconds": time.monotonic()-begin})
        finally:
            tmp.unlink(missing_ok=True)
    return {"status": "COMPLETED", "path": str(output), "reused": False}


def compute(args):
    for name in ("model", "scenario", "patch"):
        if not getattr(args, name) or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-" for c in getattr(args, name)) or getattr(args, name) in {".", ".."}:
            raise ValueError(f"unsafe {name}")
    if args.processes <= 0 or args.time_chunk <= 0 or min(args.tile_shape) <= 0 or not 0 <= args.compress_level <= 9:
        raise ValueError("invalid processes/chunk/tile/compression")
    if not isinstance(args.input_units, dict):
        unit_map = {}
        for spec in args.input_units:
            key, sep, value = spec.partition("=")
            if not sep or key not in io.VARS[args.tech] or key in unit_map:
                raise ValueError(f"invalid/duplicate input units: {spec}")
            io.conversion(key, value)
            unit_map[key] = value
        args.input_units = unit_map
    unit = Path(args.output_root).expanduser().resolve() / args.model / args.scenario / args.patch
    unit.mkdir(parents=True, exist_ok=True)
    with io.output_lock(unit / f".{args.tech}.lock"):
        plan = io.load_plan(args)
        provenance = {"implementation": implementation_identity(args.tech), "inputs": plan["inputs"],
                      "grid_fingerprint": plan["grid_fingerprint"], "years": args.years,
                      "model": args.model, "scenario": args.scenario, "patch": args.patch, "tech": args.tech,
                      "schema": 2, "input_mode": plan["input_mode"],
                      "tile_shape": list(args.tile_shape), "time_chunk": args.time_chunk,
                      "compress_level": args.compress_level, "input_units": args.input_units,
                      "time_axes": {v: [{"sha256": io.array_digest(e["raw"]), "units": e["units"],
                                         "calendar": e["calendar"], "unit_scale": e["scale"],
                                         "unit_offset": e["offset"]} for e in entries]
                                    for v, entries in plan["blocks"].items()}}
        identity = io.digest(provenance)
        # Reject conflicting runs even when identity-specific block directories differ.
        run_manifest = unit / args.tech / "manifest.json"
        if run_manifest.exists():
            try:
                old = io.read_json(run_manifest)
            except (ValueError, OSError):
                old = {}
            if old.get("identity") not in (None, identity) and not args.overwrite:
                raise ValueError("existing unit has different inputs/configuration; use a new output root or --overwrite")
        parts = (Path(args.parts_root).expanduser().resolve() / args.model / args.scenario / args.patch / args.tech
                 if args.parts_root else unit / args.tech / "blocks") / identity
        parts.mkdir(parents=True, exist_ok=True)
        paths = [parts / f"cf_{c['start_year']}-{c['end_year']}.nc" for c in plan["contracts"]]
        state = {"identity": identity, "status": "RUNNING", "merge_final": args.merge_final,
                 "provenance": provenance, "blocks": [], "final": {"status": "not_requested" if not args.merge_final else "pending"}}
        def save():
            io.atomic_json(run_manifest, state)
            io.atomic_json(parts / "manifest.json", state)
        payloads = []
        for i, path in enumerate(paths):
            reusable = io.completed(path, plan, block_axis(plan, args.tech, i), args.tech, block_identity(identity, i)) and not args.overwrite
            state["blocks"].append({"index": i, **plan["contracts"][i], "path": str(path),
                                    "status": "COMPLETED" if reusable else "pending", "reused": reusable})
            if not reusable:
                payloads.append((args, plan, i, str(path), block_identity(identity, i), provenance))
        save()
        try:
            count = min(args.processes, len(payloads), 8)
            if count <= 1:
                for payload in payloads:
                    result = _compute_block(payload)
                    state["blocks"][result["index"]].update(result, status="COMPLETED")
                    save()
            else:
                with ProcessPoolExecutor(max_workers=count, mp_context=multiprocessing.get_context("spawn")) as pool:
                    pending = [pool.submit(_compute_block, payload) for payload in payloads]
                    for future in as_completed(pending):
                        result = future.result()
                        state["blocks"][result["index"]].update(result, status="COMPLETED")
                        save()
            for i, path in enumerate(paths):
                if not io.completed(path, plan, block_axis(plan, args.tech, i), args.tech, block_identity(identity, i)):
                    raise ValueError(f"block validation failed: {path}")
            io.unchanged(plan["inputs"])
            if args.merge_final:
                state["status"] = "merge_pending"; save()
                state["final"] = _merge(args, plan, paths, unit / (args.tech+".nc"), identity, provenance)
            state["status"] = "COMPLETED"; save()
        except BaseException as exc:
            state["status"] = "FAILED"
            state["error"] = str(exc)
            save()
            raise
        return run_manifest


if __name__ == "__main__":
    print(compute(parser().parse_args()))
