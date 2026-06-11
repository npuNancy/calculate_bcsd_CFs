"""ERA5-Land 全球逐月数据通用工具函数。"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import xarray as xr
from netCDF4 import Dataset

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def parse_months(months: str) -> list[int]:
    if not months.strip():
        return list(range(1, 13))
    out = [int(m) for m in months.split(",") if m.strip()]
    bad = [m for m in out if m < 1 or m > 12]
    if bad:
        raise ValueError(f"月份必须在 1-12 之间，收到：{bad}")
    return out


def iter_year_months(start_year: int, end_year: int, months: list[int] | None = None) -> list[tuple[int, int]]:
    ms = months if months is not None else list(range(1, 13))
    return [(y, m) for y in range(start_year, end_year + 1) for m in ms]


def era5land_file(data_dir: str | Path, var: str, year: int, month: int) -> Path:
    return Path(data_dir) / var / f"{var}_{year:04d}_{month:02d}.nc"


def find_coord_name(ds: xr.Dataset, candidates: Iterable[str]) -> str:
    for c in candidates:
        if c in ds.coords or c in ds.dims:
            return c
    raise KeyError(f"找不到坐标名，候选：{list(candidates)}；数据维度：{list(ds.dims)}")


def prepare_dataarray(
    ds: xr.Dataset,
    preferred_var: str,
    time_name: str,
    lat_name: str,
    lon_name: str,
) -> xr.DataArray:
    var = preferred_var if preferred_var in ds.data_vars else (
        list(ds.data_vars)[0] if len(ds.data_vars) == 1 else None
    )
    if var is None:
        raise KeyError(f"无法确定变量 {preferred_var}，可用变量：{list(ds.data_vars)}")
    da = ds[var]
    for dim in list(da.dims):
        if dim not in {time_name, lat_name, lon_name}:
            if da.sizes[dim] == 1:
                da = da.isel({dim: 0}, drop=True)
            else:
                raise ValueError(f"变量 {var} 存在非单例额外维度 {dim}={da.sizes[dim]}")
    return da.transpose(time_name, lat_name, lon_name)


def open_era5land(path: Path, var: str) -> tuple[xr.Dataset, str, str, str]:
    ds = xr.open_dataset(path)
    time_name = find_coord_name(ds, ["valid_time", "time"])
    lat_name = find_coord_name(ds, ["latitude", "lat"])
    lon_name = find_coord_name(ds, ["longitude", "lon"])
    return ds, time_name, lat_name, lon_name


def get_lat_lon(ds: xr.Dataset, lat_name: str, lon_name: str) -> tuple[np.ndarray, np.ndarray]:
    return ds[lat_name].values.astype(np.float32), ds[lon_name].values.astype(np.float32)


def resolve_grid(
    data_dir: str, start_year: int, end_year: int, months: list[int], probe_var: str = "u10",
) -> tuple[np.ndarray, np.ndarray]:
    for year, month in iter_year_months(start_year, end_year, months):
        p = era5land_file(data_dir, probe_var, year, month)
        if p.exists():
            ds, _, latname, lonname = open_era5land(p, probe_var)
            lats, lons = get_lat_lon(ds, latname, lonname)
            ds.close()
            return lats, lons
    raise RuntimeError("找不到任何 ERA5-Land 文件，请检查 data_dir 和年份范围")


def save_result(
    out_file: Path,
    data: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
    var_name: str,
    long_name: str,
    units: str,
    attrs: dict,
    compress_level: int = 4,
) -> Path:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_file = out_file.with_name(f".{out_file.name}.tmp")
    if tmp_file.exists():
        tmp_file.unlink()

    nc = Dataset(tmp_file, "w", format="NETCDF4")
    nc.createDimension("lat", len(lats))
    nc.createDimension("lon", len(lons))

    latvar = nc.createVariable("lat", "f4", ("lat",))
    lonvar = nc.createVariable("lon", "f4", ("lon",))
    latvar[:] = lats
    lonvar[:] = lons
    latvar.units = "degrees_north"
    latvar.long_name = "latitude"
    lonvar.units = "degrees_east"
    lonvar.long_name = "longitude"

    v = nc.createVariable(
        var_name, "f4", ("lat", "lon"),
        zlib=True, complevel=compress_level,
        fill_value=np.float32(np.nan),
    )
    v[:] = data.astype(np.float32)
    v.long_name = long_name
    v.units = units

    for k, val in attrs.items():
        try:
            setattr(nc, k, val)
        except Exception:
            pass

    nc.close()
    os.replace(tmp_file, out_file)
    logger.info(f"已保存：{out_file}")
    return out_file
