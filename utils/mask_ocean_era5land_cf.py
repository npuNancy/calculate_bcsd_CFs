"""
将 ERA5-Land 容量因子文件中的海洋格点（值为 0）改为 NaN
========================================================

适用对象
--------
仅处理 ERA5-Land 容量因子输出（S02E01 / S02E02 产生的规则网格文件），
即坐标为一维 lat/lon、变量为 solar_cf 或 wind_cf 的 NetCDF。

背景
----
早期版本的 CF 计算脚本在单位转换阶段用 np.nan_to_num 把海洋的 NaN
填成了 0；ERA5-Land 海洋约占 66%，导致全局/区域平均被严重稀释。
本脚本作为「后处理」工具，用 global-land-mask 判断海陆，把所有海洋
格点重新设为 NaN，而不重新计算容量因子。

处理逻辑
--------
  1. 读取输入文件的一维 lat/lon，用 global-land-mask 构建陆地掩膜；
  2. 逐时间块读取 CF 变量，把海洋格点（掩膜为 False）置为 NaN；
  3. 逐块写入输出文件，避免一次性载入整个大数组。

注意：脚本只改海洋格点，陆地格点的数值原样保留。

示例
----
# 单文件
python utils/mask_ocean_era5land_cf.py \
  --input  output/CFs_of_solar_ERA5Land/solar_cf_2015_06.nc \
  --output output/CFs_of_solar_ERA5Land_masked/solar_cf_2015_06.nc

# 原地覆盖（先写临时文件再替换，安全）
python utils/mask_ocean_era5land_cf.py \
  --input  output/CFs_of_wind_ERA5Land/wind_cf_2015_06.nc \
  --inplace
"""

from __future__ import annotations

import argparse
import gc
import logging
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import xarray as xr
from netCDF4 import Dataset
from tqdm import tqdm
from global_land_mask import globe

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ERA5-Land CF 文件可能出现的坐标 / 变量命名
TIME_CANDIDATES = ["time", "valid_time"]
LAT_CANDIDATES = ["lat", "latitude"]
LON_CANDIDATES = ["lon", "longitude"]
CF_CANDIDATES = ["solar_cf", "wind_cf"]


def find_name(ds: xr.Dataset, candidates: Iterable[str], *, kind: str) -> str:
    """在坐标/维度或数据变量里找出第一个匹配的名字。"""
    for c in candidates:
        if c in ds.coords or c in ds.dims or c in ds.data_vars:
            return c
    raise KeyError(f"找不到{kind}，候选：{list(candidates)}；可用：{list(ds.variables)}")


def get_cf_var_name(ds: xr.Dataset) -> str:
    """获取容量因子变量名：优先 solar_cf/wind_cf，否则取唯一数据变量。"""
    for c in CF_CANDIDATES:
        if c in ds.data_vars:
            return c
    if len(ds.data_vars) == 1:
        return list(ds.data_vars)[0]
    raise KeyError(f"无法确定 CF 变量，候选：{CF_CANDIDATES}，可用：{list(ds.data_vars)}")


def build_land_mask(lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    """根据一维 lat/lon 生成陆地掩膜 (n_lat, n_lon)，True=陆地。

    global_land_mask 要求经度 ∈ [-180, 180]、纬度 ∈ [-90, 90]，
    因此对 ERA5-Land [0, 360) 经度做归一化，并夹紧纬度避免浮点越界。
    """
    lats = np.asarray(lats, dtype=np.float64)
    lons = np.asarray(lons, dtype=np.float64)

    lons_conv = ((lons + 180.0) % 360.0) - 180.0  # [0,360) -> [-180,180)
    lats_clip = np.clip(lats, -90.0, 90.0)

    lon2d, lat2d = np.meshgrid(lons_conv, lats_clip)  # (n_lat, n_lon)
    return globe.is_land(lat2d, lon2d).astype(bool)


def mask_ocean_file(
    input_file: Path,
    output_file: Path,
    *,
    chunk_time: int = 24,
    compress_level: int = 4,
    overwrite: bool = False,
) -> Path:
    """把单个 ERA5-Land CF 文件的海洋格点（0）改为 NaN，写入新文件。"""
    if output_file.exists() and not overwrite:
        logger.info(f"✓ 输出已存在，跳过：{output_file}")
        return output_file

    # decode_times=False 保留原始数值型时间，避免写出时单位/日历问题
    ds = xr.open_dataset(input_file, decode_times=False)

    time_name = find_name(ds, TIME_CANDIDATES, kind="时间坐标")
    lat_name = find_name(ds, LAT_CANDIDATES, kind="纬度坐标")
    lon_name = find_name(ds, LON_CANDIDATES, kind="经度坐标")
    cf_name = get_cf_var_name(ds)

    if ds[lat_name].ndim != 1 or ds[lon_name].ndim != 1:
        ds.close()
        raise ValueError(
            f"本脚本仅支持一维 lat/lon 的 ERA5-Land 规则网格，"
            f"当前 lat.ndim={ds[lat_name].ndim}, lon.ndim={ds[lon_name].ndim}"
        )

    da = ds[cf_name].transpose(time_name, lat_name, lon_name)
    n_time = da.sizes[time_name]
    lats = ds[lat_name].values.astype(np.float32)
    lons = ds[lon_name].values.astype(np.float32)
    n_lat, n_lon = len(lats), len(lons)
    logger.info(f"输入：{input_file}")
    logger.info(f"变量：{cf_name}；维度：time={n_time}, lat={n_lat}, lon={n_lon}")

    land_mask = build_land_mask(lats, lons)  # (n_lat, n_lon) bool
    logger.info(f"陆地格点占比：{100.0 * land_mask.mean():.1f}%")

    # ── 准备输出文件（复制坐标与属性，CF 变量重建） ──────────────────
    output_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_file = output_file.with_name(f".{output_file.name}.tmp")
    if tmp_file.exists():
        tmp_file.unlink()

    raw_time = ds[time_name].values[:n_time]
    time_attrs = dict(ds[time_name].attrs)
    lat_attrs = dict(ds[lat_name].attrs)
    lon_attrs = dict(ds[lon_name].attrs)
    cf_attrs = dict(ds[cf_name].attrs)
    global_attrs = dict(ds.attrs)

    nc = Dataset(tmp_file, "w", format="NETCDF4")
    nc.createDimension("time", n_time)
    nc.createDimension("lat", n_lat)
    nc.createDimension("lon", n_lon)

    tvar = nc.createVariable("time", raw_time.dtype, ("time",))
    tvar[:] = raw_time
    for k, v in time_attrs.items():
        try:
            setattr(tvar, k, v)
        except Exception:
            pass

    latvar = nc.createVariable("lat", "f4", ("lat",))
    lonvar = nc.createVariable("lon", "f4", ("lon",))
    latvar[:] = lats
    lonvar[:] = lons
    for k, v in lat_attrs.items():
        try:
            setattr(latvar, k, v)
        except Exception:
            pass
    for k, v in lon_attrs.items():
        try:
            setattr(lonvar, k, v)
        except Exception:
            pass

    chunksizes = (min(chunk_time, n_time), n_lat, n_lon)
    cf_out = nc.createVariable(
        cf_name,
        "f4",
        ("time", "lat", "lon"),
        zlib=True,
        complevel=compress_level,
        chunksizes=chunksizes,
        fill_value=np.float32(np.nan),
    )
    for k, v in cf_attrs.items():
        if k == "_FillValue":
            continue
        try:
            setattr(cf_out, k, v)
        except Exception:
            pass

    for k, v in global_attrs.items():
        try:
            setattr(nc, k, v)
        except Exception:
            pass
    nc.ocean_handling = "ocean grid cells set to NaN via global-land-mask (post-processing)"

    # ── 逐时间块改写：海洋格点 -> NaN ───────────────────────────────
    ocean = ~land_mask  # (n_lat, n_lon)
    try:
        for start in tqdm(range(0, n_time, chunk_time), desc=input_file.name, unit="块"):
            end = min(start + chunk_time, n_time)
            block = da.isel({time_name: slice(start, end)}).values.astype(np.float32)
            block[:, ocean] = np.nan
            cf_out[start:end, :, :] = block
            del block
            gc.collect()
    finally:
        nc.close()
        ds.close()

    os.replace(tmp_file, output_file)
    logger.info(f"✓ 已保存：{output_file}")
    return output_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description="把 ERA5-Land 容量因子文件中的海洋格点(0)改为 NaN（global-land-mask 判断海陆）"
    )
    parser.add_argument("--input", required=True, help="输入 NetCDF 文件（ERA5-Land 的 solar_cf / wind_cf）")
    parser.add_argument("--output", help="输出 NetCDF 文件；与 --inplace 二选一")
    parser.add_argument(
        "--inplace",
        action="store_true",
        help="原地覆盖输入文件（内部先写临时文件再替换）；与 --output 二选一",
    )
    parser.add_argument("--chunk_time", type=int, default=24, help="每块处理的时间步数，默认 24")
    parser.add_argument("--compress_level", type=int, default=4, help="NetCDF 压缩级别 0-9，默认 4")
    parser.add_argument("--overwrite", action="store_true", help="若输出文件已存在则覆盖")
    args = parser.parse_args()

    input_file = Path(args.input)
    if not input_file.exists():
        raise FileNotFoundError(f"找不到输入文件：{input_file}")

    if args.inplace and args.output:
        raise SystemExit("--inplace 与 --output 不能同时指定。")
    if not args.inplace and not args.output:
        raise SystemExit("必须指定 --output 或 --inplace 其一。")

    output_file = input_file if args.inplace else Path(args.output)
    # 原地模式总是覆盖
    overwrite = args.overwrite or args.inplace

    mask_ocean_file(
        input_file=input_file,
        output_file=output_file,
        chunk_time=args.chunk_time,
        compress_level=args.compress_level,
        overwrite=overwrite,
    )

    logger.info("✓ 完成！")


if __name__ == "__main__":
    main()
