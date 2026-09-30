"""Validated BCSD point-block inputs and streaming native-grid NetCDF output."""
from __future__ import annotations

from contextlib import contextmanager, ExitStack
import fcntl
import hashlib
import json
import os
from pathlib import Path
import uuid

import netCDF4 as nc
import numpy as np

VARS = {"wind": ("uas", "vas"), "solar": ("rsds", "tas", "uas", "vas")}
FILL = np.float32(9.96921e36)
TIME_UNITS = "hours since 1970-01-01 00:00:00"
CALENDARS = {"standard": "gregorian", "gregorian": "gregorian",
             "proleptic_gregorian": "proleptic_gregorian", "365_day": "noleap", "noleap": "noleap"}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def array_digest(*arrays):
    h = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(a)
        h.update(str((a.shape, a.dtype.str)).encode())
        h.update(a.tobytes())
    return h.hexdigest()


def read_json(path):
    value = json.loads(Path(path).read_text())
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def file_identity(path, content=False):
    path = Path(path).resolve()
    stat = path.stat()
    result = {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if content:
        result["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def unchanged(identities):
    for old in identities:
        if file_identity(old["path"], "sha256" in old) != old:
            raise ValueError(f"input changed during run: {old['path']}")


def partial_path(path):
    return Path(str(path) + f".partial.{os.getpid()}.{uuid.uuid4().hex}")


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = partial_path(path)
    try:
        tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


@contextmanager
def output_lock(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"output already locked: {path}") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def numbers(raw, units, calendar):
    dates = nc.num2date(raw, units, calendar, only_use_cftime_datetimes=True)
    return np.asarray(nc.date2num(dates, TIME_UNITS, CALENDARS[calendar]), dtype=np.float64)


def conversion(variable, units):
    u = units.lower().replace(" ", "").replace("**", "^").replace("°", "deg")
    if variable in {"uas", "vas"} and u in {"ms-1", "ms^-1", "m/s"}:
        return 1.0, 0.0
    if variable == "rsds":
        if u in {"wm-2", "wm^-2", "w/m2", "w/m^2"}:
            return 0.001, 0.0
        if u in {"kwm-2", "kwm^-2", "kw/m2", "kw/m^2"}:
            return 1.0, 0.0
    if variable == "tas":
        if u in {"k", "kelvin"}:
            return 1.0, -273.15
        if u in {"c", "degc", "celsius", "degree_celsius", "degrees_celsius"}:
            return 1.0, 0.0
    raise ValueError(f"unsupported/missing units for {variable}: {units!r}")


def _path(value, root):
    p = Path(value).expanduser()
    return (p if p.is_absolute() else root / p).resolve()


def production_block_path(value, root, model, scenario, variable, patch):
    # Production manifests can retain their original absolute root after relocation.
    parts = Path(value).parts
    suffix = ("blocks", model, scenario, variable, patch)
    if len(parts) < 6 or tuple(parts[-6:-1]) != suffix:
        raise ValueError(f"unexpected block path identity: {value}")
    return (root.joinpath(*suffix, parts[-1])).resolve()


def _integers(ds, key):
    v = ds[key]
    if v.dimensions != ("point",) or v.dtype.kind not in "iu":
        raise ValueError(f"invalid point indices: {key}")
    return np.asarray(v[:], dtype=np.int64)


def _axis(ds, key):
    v = ds[key]
    a = np.asarray(v[:])
    if v.dimensions != (key,) or a.size < 1 or not np.isfinite(a).all():
        raise ValueError(f"invalid coordinate: {key}")
    if a.size > 1 and not (np.all(np.diff(a) > 0) or np.all(np.diff(a) < 0)):
        raise ValueError(f"non-monotonic coordinate: {key}")
    return a


def load_plan(args):
    """Read metadata only; validate all variables/years before starting workers."""
    root = Path(args.bcsd_root).expanduser().resolve()
    patch_path = Path(args.patch_manifest).expanduser().resolve()
    land_path = Path(args.land_plan).expanduser().resolve()
    inputs = [file_identity(patch_path, True), file_identity(land_path)]
    patch_meta = read_json(patch_path)
    w, e, south, north = map(float, patch_meta["patches"][args.patch]["core_bbox_360"])
    with nc.Dataset(land_path) as ds:
        lat, lon = _axis(ds, "lat"), _axis(ds, "lon")
        py, px = _integers(ds, "point_y"), _integers(ds, "point_x")
        if np.any(py < 0) or np.any(py >= len(lat)) or np.any(px < 0) or np.any(px >= len(lon)):
            raise ValueError("land plan point indices out of range")
        mask = np.asarray(ds["land_mask"][:])
        if ds["land_mask"].dimensions != ("lat", "lon") or not np.isin(mask, [0, 1]).all():
            raise ValueError("invalid land_mask")
    ysel = np.flatnonzero((lat >= south) & ((lat <= north) if north == 90 else (lat < north)))
    xsel = np.flatnonzero(((lon.astype(float) - w) % 360) < (e - w))
    if not (ysel.size and xsel.size) or not (0 < e-w <= 360 and -90 <= south < north <= 90):
        raise ValueError("empty/invalid patch core")
    expected_mask = mask[np.ix_(ysel, xsel)].astype(np.int8)
    ymap = np.full(len(lat), -1, dtype=np.int64); ymap[ysel] = np.arange(len(ysel))
    xmap = np.full(len(lon), -1, dtype=np.int64); xmap[xsel] = np.arange(len(xsel))
    blocks, contracts, reference_points, ref_calendar = {}, None, None, None
    for variable in VARS[args.tech]:
        manifest = root / "manifests" / f"{args.model}__{args.scenario}__{variable}__{args.patch}.json"
        meta = read_json(manifest)
        inputs.append(file_identity(manifest, True))
        if meta.get("kind") != "global-bcsd-three-stage":
            raise ValueError(f"invalid BCSD manifest: {manifest}")
        for key, value in (("model", args.model), ("scenario", args.scenario), ("variable", variable)):
            if key in meta and meta[key] != value:
                raise ValueError(f"manifest {key} mismatch: {manifest}")
        current = meta["time_block_contract"]["blocks"]
        current = [{"start_year": int(c["start_year"]), "end_year": int(c["end_year"])} for c in current]
        if not current or len(current) > 8 or any(c["start_year"] > c["end_year"] for c in current):
            raise ValueError("invalid year-block contract")
        if any(b["start_year"] != a["end_year"] + 1 for a, b in zip(current, current[1:])):
            raise ValueError("year-block contract has gaps/overlaps")
        if contracts is None:
            contracts = current
        elif contracts != current:
            raise ValueError("variables have different year-block contracts")
        tasks = [t for t in meta["final_tasks"] if t["patch_id"] == args.patch and t["variable"] == variable]
        if len(tasks) != 1 or len(tasks[0]["block_files"]) != len(contracts):
            raise ValueError("manifest final task/block count mismatch")
        plans = [meta["land_plan"]] if meta.get("land_plan") else [t["land_plan_path"] for t in meta["block_tasks"]]
        if not plans or any(_path(p, root) != land_path for p in plans):
            raise ValueError("land plan differs from manifest")
        entries = []
        for index, (contract, value) in enumerate(zip(contracts, tasks[0]["block_files"])):
            path = production_block_path(value, root, args.model, args.scenario, variable, args.patch)
            if not path.is_relative_to((root / "blocks").resolve()):
                raise ValueError(f"block outside production blocks directory: {path}")
            side = path.with_name("block_" + path.stem + ".json")
            sm = read_json(side)
            if (sm.get("kind") != "global-step6-block" or sm.get("variable") != variable
                    or args.patch not in sm.get("patch_ids", [])
                    or _path(sm.get("outputs", {}).get(args.patch, ""), root) != path):
                raise ValueError(f"block sidecar identity mismatch: {side}")
            inputs.extend([file_identity(path), file_identity(side, True)])
            with nc.Dataset(path) as ds:
                if getattr(ds, "patch_id", None) != args.patch or ds[variable].dimensions != ("time", "point"):
                    raise ValueError(f"block schema/patch mismatch: {path}")
                yi, xi, flat, gp = [_integers(ds, k) for k in ("y_index", "x_index", "flat_index", "global_point")]
                if (np.any(yi < 0) or np.any(yi >= len(lat)) or np.any(xi < 0) or np.any(xi >= len(lon))
                        or np.any(gp < 0) or np.any(gp >= len(py))):
                    raise ValueError(f"point indices out of range: {path}")
                if (len(np.unique(flat)) != len(flat) or not np.array_equal(flat, yi * len(lon) + xi)
                        or not np.array_equal(py[gp], yi) or not np.array_equal(px[gp], xi)
                        or not np.array_equal(ds["lat"][:], lat[yi]) or not np.array_equal(ds["lon"][:], lon[xi])):
                    raise ValueError(f"point identity mismatch: {path}")
                points = (yi, xi, flat, gp)
                if reference_points is None:
                    reference_points = points
                elif any(not np.array_equal(a, b) for a, b in zip(reference_points, points)):
                    raise ValueError(f"point order differs: {path}")
                t = ds["time"]
                raw = np.asarray(t[:])
                units, calendar = getattr(t, "units", ""), getattr(t, "calendar", "standard")
                if calendar not in CALENDARS or t.dimensions != ("time",) or len(raw) < 2 or not np.isfinite(raw).all():
                    raise ValueError(f"invalid calendar/time axis: {path}")
                cal = CALENDARS[calendar]
                if ref_calendar is None:
                    ref_calendar = cal
                elif ref_calendar != cal:
                    raise ValueError("incompatible calendars across input blocks")
                numeric = numbers(raw, units, calendar)
                if not np.allclose(np.diff(numeric), 3.0, rtol=0, atol=1e-8):
                    raise ValueError(f"time gap/duplicate or non-3-hour axis: {path}")
                dates = nc.num2date(raw[[0, -1]], units, calendar)
                if (dates[0].year, dates[1].year) != (contract["start_year"], contract["end_year"]):
                    raise ValueError(f"block years mismatch: {path}")
                import cftime
                begin = nc.date2num(cftime.datetime(contract["start_year"], 1, 1, calendar=cal), TIME_UNITS, cal)
                end = nc.date2num(cftime.datetime(contract["end_year"] + 1, 1, 1, calendar=cal), TIME_UNITS, cal)
                if not (0 <= numeric[0]-begin <= 3 and 0 < end-numeric[-1] <= 3):
                    raise ValueError(f"incomplete year coverage: {path}")
                declared_units = getattr(ds[variable], "units", "")
                supplied_units = getattr(args, "input_units", {}).get(variable)
                if declared_units and supplied_units and conversion(variable, declared_units) != conversion(variable, supplied_units):
                    raise ValueError(f"explicit units conflict with file metadata: {path}")
                scale, offset = conversion(variable, declared_units or supplied_units or "")
                if "time_start" in sm or "time_stop" in sm:
                    if (int(sm["time_stop"])-int(sm["time_start"]) != len(raw)
                            or getattr(ds, "time_start_index", None) != sm["time_start"]
                            or getattr(ds, "time_stop_index", None) != sm["time_stop"]):
                        raise ValueError(f"sidecar time indices mismatch: {side}")
                land_input = sm.get("inputs", {}).get("land_plan")
                if land_input:
                    expected = file_identity(land_path)
                    if (_path(land_input["path"], root) != land_path or
                            any(land_input[k] != expected[k] for k in ("size", "mtime_ns"))):
                        raise ValueError(f"sidecar land plan identity mismatch: {side}")
                if entries and not np.isclose(numeric[0]-entries[-1]["numbers"][-1], 3, rtol=0, atol=1e-8):
                    raise ValueError(f"time gap/overlap across blocks: {path}")
                if entries and (scale, offset) != (entries[0]["scale"], entries[0]["offset"]):
                    raise ValueError(f"units change across year blocks: {path}")
                entries.append({"path": str(path), "raw": raw, "units": units, "calendar": calendar,
                                "numbers": numeric, "scale": scale, "offset": offset})
        blocks[variable] = entries
    year_range = f"{contracts[0]['start_year']}-{contracts[-1]['end_year']}"
    if args.years != year_range:
        raise ValueError(f"--years must match complete manifest contract: {year_range}")
    yi, xi, _, _ = reference_points
    ly, lx = ymap[yi], xmap[xi]
    if np.any(ly < 0) or np.any(lx < 0):
        raise ValueError("block point outside patch core")
    domain = np.zeros(expected_mask.shape, dtype=np.int8); domain[ly, lx] = 1
    if not np.array_equal(domain, expected_mask):
        raise ValueError("block points differ from land plan domain")
    return {"contracts": contracts, "blocks": blocks, "inputs": inputs,
            "lat": lat[ysel], "lon": lon[xsel], "mask": domain, "y": ly, "x": lx,
            "grid_fingerprint": array_digest(lat[ysel], lon[xsel], domain)}


def tiles(plan, shape):
    for y0 in range(0, len(plan["lat"]), shape[0]):
        for x0 in range(0, len(plan["lon"]), shape[1]):
            yield (y0, min(y0+shape[0], len(plan["lat"])), x0, min(x0+shape[1], len(plan["lon"])))


class BlockReader:
    """Read contiguous point runs, using actual time indices across year boundaries."""
    def __init__(self, entries, variable):
        self.entries, self.variable = entries, variable
        self.stack = ExitStack()
        self.opened = {}
        self.offsets = np.r_[0, np.cumsum([len(e["numbers"]) for e in entries])]
        self.times = np.concatenate([e["numbers"] for e in entries])
        self.read_bytes = 0

    def close(self):
        self.stack.close()

    def read(self, targets, positions):
        right = np.searchsorted(self.times, targets)
        right = np.clip(right, 0, len(self.times)-1)
        exact = np.isclose(self.times[right], targets, rtol=0, atol=1e-8)
        left = np.where(exact, right, np.maximum(right-1, 0))
        valid = (targets >= self.times[0]-1e-8) & (targets <= self.times[-1]+1e-8)
        lo, hi = int(left.min()), int(right.max())+1
        source = np.empty((hi-lo, len(positions)), dtype=np.float32)
        # positions are sorted; adjacent disk columns are read in a single slice.
        edges = np.r_[0, np.flatnonzero(np.diff(positions) != 1)+1, len(positions)]
        for i, e in enumerate(self.entries):
            a, b = max(lo, int(self.offsets[i])), min(hi, int(self.offsets[i+1]))
            if a >= b:
                continue
            if i not in self.opened:
                self.opened[i] = self.stack.enter_context(nc.Dataset(e["path"]))
                self.opened[i][self.variable].set_var_chunk_cache(64*1024**2, 1009, 0.75)
            v = self.opened[i][self.variable]
            for r0, r1 in zip(edges[:-1], edges[1:]):
                raw = v[a-self.offsets[i]:b-self.offsets[i], int(positions[r0]):int(positions[r1-1])+1]
                values = np.asarray(np.ma.filled(raw, np.nan), dtype=np.float32)
                values[~np.isfinite(values)] = np.nan
                source[a-lo:b-lo, r0:r1] = values
                self.read_bytes += values.nbytes
        lval, rval = source[left-lo], source[right-lo]
        denom = self.times[right]-self.times[left]
        weight = np.divide(targets-self.times[left], denom, out=np.zeros_like(targets), where=denom != 0)
        result = np.where(exact[:, None], rval, lval + (rval-lval)*weight[:, None]).astype(np.float32)
        result[~valid] = np.nan
        if self.entries[0]["scale"] == 0.001:
            result /= 1000.0
        if self.entries[0]["offset"]:
            result += self.entries[0]["offset"]
        return result


def create_output(path, plan, axis, args, identity):
    ds = nc.Dataset(path, "w", format="NETCDF4")
    try:
        for name, size in (("time", len(axis["raw"])), ("lat", len(plan["lat"])), ("lon", len(plan["lon"]))):
            ds.createDimension(name, size)
        for name in ("lat", "lon"):
            v = ds.createVariable(name, plan[name].dtype, (name,))
            v[:] = plan[name]
            v.standard_name = "latitude" if name == "lat" else "longitude"
            v.units = "degrees_north" if name == "lat" else "degrees_east"
        t = ds.createVariable("time", axis["raw"].dtype, ("time",))
        t[:] = axis["raw"]; t.units = axis["units"]; t.calendar = axis["calendar"]; t.standard_name = "time"
        m = ds.createVariable("domain_mask", "i1", ("lat", "lon"), zlib=True, complevel=1)
        m[:] = plan["mask"]
        m.flag_values = np.array([0, 1], dtype="i1"); m.flag_meanings = "outside_domain inside_domain"
        v = ds.createVariable(args.tech+"_cf", "f4", ("time", "lat", "lon"), fill_value=FILL,
                              zlib=args.compress_level > 0, complevel=args.compress_level, shuffle=True,
                              chunksizes=(min(args.time_chunk, len(axis["raw"])),
                                          min(args.tile_shape[0], len(plan["lat"])), min(args.tile_shape[1], len(plan["lon"]))))
        v.units = "1"; v.long_name = args.tech + " capacity factor"; v.valid_range = np.array([0, 1], dtype="f4")
        ds.setncatts({"schema_version": 1, "model": args.model, "scenario": args.scenario, "patch_id": args.patch,
                     "tech": args.tech, "identity": identity, "grid_fingerprint": plan["grid_fingerprint"],
                     "input_mode": "blocks", "missing_policy": "outside domain or any required input nonfinite"})
        return ds
    except BaseException:
        ds.close()
        raise


def validate_output(path, plan, axis, tech, identity):
    with nc.Dataset(path) as ds:
        if ds.getncattr("identity") != identity:
            raise ValueError("output identity mismatch")
        for name, expected in (("lat", plan["lat"]), ("lon", plan["lon"]), ("time", axis["raw"]), ("domain_mask", plan["mask"])):
            if not np.array_equal(ds[name][:], expected) or ds[name].dtype != expected.dtype:
                raise ValueError(f"output coordinate/mask mismatch: {name}")
        v = ds[tech+"_cf"]
        if (v.dimensions != ("time", "lat", "lon") or v.dtype != np.dtype("f4")
                or getattr(v, "_FillValue", None) != FILL or v.shape != (len(axis["raw"]), len(plan["lat"]), len(plan["lon"]))):
            raise ValueError("output CF schema mismatch")
        if ds["time"].units != axis["units"] or ds["time"].calendar != axis["calendar"]:
            raise ValueError("output time encoding mismatch")


def completed(path, plan, axis, tech, identity):
    try:
        side = read_json(str(path)+".json")
        if side["status"] != "COMPLETED" or side["identity"] != identity or side["file"] != file_identity(path):
            return False
        validate_output(path, plan, axis, tech, identity)
        return True
    except (OSError, ValueError, KeyError, AttributeError, IndexError, RuntimeError, TypeError):
        return False


def publish(tmp, path, plan, axis, tech, identity, metadata):
    validate_output(tmp, plan, axis, tech, identity)
    os.replace(tmp, path)
    atomic_json(str(path)+".json", {**metadata, "status": "COMPLETED", "identity": identity, "file": file_identity(path)})
