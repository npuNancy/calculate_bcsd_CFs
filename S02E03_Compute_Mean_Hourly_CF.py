"""
多年平均逐小时容量因子计算（ERA5-Land 风光共用）
==================================================

从逐月保存的容量因子文件计算多年平均逐小时容量因子。
对每个月份（1-12），逐年对齐小时，逐格点求平均。
2 月仅保留前 28 天（672 小时），确保闰年与非闰年对齐。
最终拼接 12 个月得到 8760 小时（非闰年）的多年平均逐小时容量因子。

步骤
----
  1. 逐月计算多年平均 → 保存至 monthly/ 子目录
  2. 可选：合并 12 个月为一个年度文件

输入
----
  {cf_dir}/{solar|wind}_cf_YYYY_MM.nc
  变量名：solar_cf 或 wind_cf

输出路径
--------
  {save_dir}/mean_CFs_of_{solar|wind}_ERA5Land/{start_year}_{end_year}/
    monthly/mean_{solar|wind}_cf_MM.nc
    mean_{solar|wind}_cf_{start_year}_{end_year}.nc

示例
----
# 太阳能
python S02E03_Compute_Mean_Hourly_CF.py \
  --cf_dir output/CFs_of_solar_ERA5Land \
  --type solar \
  --start_year 2015 --end_year 2025 \
  --save_dir output

# 风能，不合并
python S02E03_Compute_Mean_Hourly_CF.py \
  --cf_dir output/CFs_of_wind_ERA5Land \
  --type wind \
  --start_year 2015 --end_year 2025 \
  --save_dir output \
  --no_merge
"""

from __future__ import annotations

import argparse
import calendar
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

# 非闰年每月天数
DAYS_PER_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def find_coord_name(ds: xr.Dataset, candidates: Iterable[str]) -> str:
    for c in candidates:
        if c in ds.coords or c in ds.dims:
            return c
    raise KeyError(f"找不到坐标名，候选：{list(candidates)}；数据维度：{list(ds.dims)}")


def get_var_name(ds: xr.Dataset, preferred: str) -> str:
    if preferred in ds.data_vars:
        return preferred
    if len(ds.data_vars) == 1:
        return list(ds.data_vars)[0]
    raise KeyError(f"无法确定变量 {preferred}，可用变量：{list(ds.data_vars)}")


def build_land_mask(lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    """根据输出网格 lat/lon 生成陆地掩膜 (n_lat, n_lon)，True=陆地。

    输入数据（output/CFs_of_*_ERA5Land）的海洋格点曾被错误填成 0，
    本脚本在求多年平均后用 global-land-mask 把海洋格点重新设为 NaN。
    global_land_mask 要求经度 ∈ [-180, 180]、纬度 ∈ [-90, 90]，
    因此对 ERA5-Land [0, 360) 经度做归一化，并夹紧纬度避免浮点越界。
    """
    lats = np.asarray(lats, dtype=np.float64)
    lons = np.asarray(lons, dtype=np.float64)

    lons_conv = ((lons + 180.0) % 360.0) - 180.0  # [0,360) -> [-180,180)
    lats_clip = np.clip(lats, -90.0, 90.0)

    lon2d, lat2d = np.meshgrid(lons_conv, lats_clip)  # (n_lat, n_lon)
    return globe.is_land(lat2d, lon2d).astype(bool)


def cf_file_path(cf_dir: str | Path, cf_type: str, year: int, month: int) -> Path:
    """构建容量因子文件路径。"""
    return Path(cf_dir) / f"{cf_type}_cf_{year:04d}_{month:02d}.nc"


def reference_year(start_year: int, end_year: int) -> int:
    """获取参考年份（优先非闰年），用于输出时间轴。"""
    for y in range(start_year, end_year + 1):
        if not calendar.isleap(y):
            return y
    return start_year


def output_base_dir(save_dir: str, cf_type: str, start_year: int, end_year: int) -> Path:
    """构建输出根目录。"""
    return Path(save_dir) / f"mean_CFs_of_{cf_type}_ERA5Land" / f"{start_year}_{end_year}"


def prepare_dataarray(ds: xr.Dataset, var_name: str, time_name: str, lat_name: str, lon_name: str) -> xr.DataArray:
    """获取变量，压缩额外单例维度。"""
    var = get_var_name(ds, var_name)
    da = ds[var]
    for dim in list(da.dims):
        if dim not in {time_name, lat_name, lon_name}:
            if da.sizes[dim] == 1:
                da = da.isel({dim: 0}, drop=True)
            else:
                raise ValueError(f"变量 {var} 存在非单例额外维度 {dim}={da.sizes[dim]}")
    return da


# ─────────────────────────────────────────────────────────────────
# 1. 逐月多年平均
# ─────────────────────────────────────────────────────────────────


def compute_monthly_mean(
    cf_dir: str,
    cf_type: str,
    start_year: int,
    end_year: int,
    month: int,
    output_file: Path,
    chunk_time: int = 24,
    overwrite: bool = False,
    compress_level: int = 4,
) -> Path:
    """计算单个月的多年平均逐小时容量因子。"""
    var_name = f"{cf_type}_cf"
    target_hours = DAYS_PER_MONTH[month - 1] * 24

    # 收集可用文件
    files = []
    for year in range(start_year, end_year + 1):
        f = cf_file_path(cf_dir, cf_type, year, month)
        if f.exists():
            files.append((year, f))
        else:
            logger.warning(f"文件不存在，跳过：{f}")

    if not files:
        raise FileNotFoundError(f"没有找到 {cf_type} {month:02d} 月的任何文件（范围 {start_year}-{end_year}）")

    logger.info(f"月份 {month:02d}：共 {len(files)} 个年份，目标 {target_hours} 小时")

    # 从第一个文件获取网格信息与时间坐标（decode_times=False 保留原始数值型时间）
    ref_ds = xr.open_dataset(files[0][1], decode_times=False)
    time_name = find_coord_name(ref_ds, ["time", "valid_time"])
    lat_name = find_coord_name(ref_ds, ["lat", "latitude"])
    lon_name = find_coord_name(ref_ds, ["lon", "longitude"])

    lats = ref_ds[lat_name].values.astype(np.float32)
    lons = ref_ds[lon_name].values.astype(np.float32)
    n_lat, n_lon = len(lats), len(lons)

    # 输入数据海洋格点曾被填成 0，这里构建陆地掩膜，写盘前把海洋置为 NaN。
    land_mask = build_land_mask(lats, lons)  # (n_lat, n_lon) bool
    ocean = ~land_mask
    logger.info(f"陆地格点占比：{100.0 * land_mask.mean():.1f}%")

    ref_time = ref_ds[time_name].values[:target_hours]
    time_attrs = dict(ref_ds[time_name].attrs)
    lat_attrs = dict(ref_ds[lat_name].attrs) if lat_name in ref_ds else {}
    lon_attrs = dict(ref_ds[lon_name].attrs) if lon_name in ref_ds else {}
    ref_ds.close()

    actual_hours = len(ref_time)

    # 输出文件
    output_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_file = output_file.with_name(f".{output_file.name}.tmp")
    if output_file.exists() and not overwrite:
        logger.info(f"✓ 已存在，跳过：{output_file}")
        return output_file
    if tmp_file.exists():
        tmp_file.unlink()

    nc = Dataset(tmp_file, "w", format="NETCDF4")
    nc.createDimension("time", actual_hours)
    nc.createDimension("lat", n_lat)
    nc.createDimension("lon", n_lon)

    tvar = nc.createVariable("time", ref_time.dtype, ("time",))
    tvar[:] = ref_time
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

    chunksizes = (min(chunk_time, actual_hours), n_lat, n_lon)
    cf_var = nc.createVariable(
        var_name,
        "f4",
        ("time", "lat", "lon"),
        zlib=True,
        complevel=compress_level,
        chunksizes=chunksizes,
        fill_value=np.float32(np.nan),
    )
    cf_var.long_name = f"Multi-year mean hourly {cf_type} capacity factor"
    cf_var.units = "1"

    nc.source = "ERA5-Land"
    nc.month = month
    nc.start_year = start_year
    nc.end_year = end_year
    nc.n_years_used = len(files)
    nc.years_used = ", ".join(str(y) for y, _ in files)

    try:
        for t_start in tqdm(range(0, actual_hours, chunk_time), desc=f"{cf_type}-{month:02d}", unit="块"):
            t_end = min(t_start + chunk_time, actual_hours)
            chunk_len = t_end - t_start

            sum_arr = np.zeros((chunk_len, n_lat, n_lon), dtype=np.float64)
            count_arr = np.zeros((chunk_len, n_lat, n_lon), dtype=np.int32)

            for year, f in files:
                try:
                    ds = xr.open_dataset(f, decode_times=False)
                    tname = find_coord_name(ds, ["time", "valid_time"])
                    lname = find_coord_name(ds, ["lat", "latitude"])
                    oname = find_coord_name(ds, ["lon", "longitude"])
                    da = prepare_dataarray(ds, var_name, tname, lname, oname)

                    file_n_hours = da.sizes[tname]
                    if t_start >= file_n_hours:
                        ds.close()
                        continue

                    file_end = min(t_end, file_n_hours)
                    read_len = file_end - t_start
                    chunk = da.isel({tname: slice(t_start, file_end)}).values.astype(np.float32)
                    valid = np.isfinite(chunk)

                    if read_len < chunk_len:
                        sum_arr[:read_len][valid] += chunk[valid].astype(np.float64)
                        count_arr[:read_len][valid] += 1
                    else:
                        sum_arr[valid] += chunk[valid].astype(np.float64)
                        count_arr[valid] += 1

                    ds.close()
                    del chunk, valid, da
                except Exception as e:
                    logger.warning(f"读取 {f} 失败：{e}")
                    continue

            mean_arr = np.where(count_arr > 0, sum_arr / count_arr, np.nan).astype(np.float32)
            # 海洋格点置为 NaN（修复输入数据把海洋填 0 的问题）
            mean_arr[:, ocean] = np.nan
            cf_var[t_start:t_end, :, :] = mean_arr

            del sum_arr, count_arr, mean_arr
            gc.collect()
    finally:
        nc.close()

    os.replace(tmp_file, output_file)
    logger.info(f"✓ 已保存月均值：{output_file}")
    return output_file


# ─────────────────────────────────────────────────────────────────
# 2. 合并 12 个月
# ─────────────────────────────────────────────────────────────────


def merge_monthly_means(
    monthly_dir: Path,
    cf_type: str,
    start_year: int,
    end_year: int,
    output_file: Path,
    overwrite: bool = False,
    compress_level: int = 4,
    chunk_time: int = 24,
) -> Path:
    """合并 12 个月的多年平均容量因子为一个年度文件。"""
    var_name = f"{cf_type}_cf"

    if output_file.exists() and not overwrite:
        logger.info(f"✓ 已存在，跳过合并：{output_file}")
        return output_file

    # 收集月文件信息
    month_info = []
    all_times = []
    n_lat = n_lon = None
    lat_vals = lon_vals = None
    lat_attrs: dict = {}
    lon_attrs: dict = {}
    time_attrs: dict = {}
    offset = 0

    for month in range(1, 13):
        f = monthly_dir / f"mean_{cf_type}_cf_{month:02d}.nc"
        if not f.exists():
            raise FileNotFoundError(f"月均值文件不存在：{f}")

        ds = xr.open_dataset(f, decode_times=False)
        n_hours = ds.sizes["time"]
        all_times.append(ds["time"].values[:])

        if n_lat is None:
            n_lat = ds.sizes["lat"]
            n_lon = ds.sizes["lon"]
            lat_vals = ds["lat"].values.astype(np.float32)
            lon_vals = ds["lon"].values.astype(np.float32)
            lat_attrs = dict(ds["lat"].attrs)
            lon_attrs = dict(ds["lon"].attrs)
            time_attrs = dict(ds["time"].attrs)
        else:
            if ds.sizes["lat"] != n_lat or ds.sizes["lon"] != n_lon:
                raise ValueError(f"月份 {month:02d} 网格尺寸不一致")

        month_info.append({"month": month, "file": f, "n_hours": n_hours, "offset": offset})
        offset += n_hours
        ds.close()

    total_hours = offset
    merged_time = np.concatenate(all_times)

    logger.info(f"合并 {total_hours} 小时（{total_hours / 24:.0f} 天）")

    # 创建输出文件
    tmp_file = output_file.with_name(f".{output_file.name}.tmp")
    if tmp_file.exists():
        tmp_file.unlink()

    nc = Dataset(tmp_file, "w", format="NETCDF4")
    nc.createDimension("time", total_hours)
    nc.createDimension("lat", n_lat)
    nc.createDimension("lon", n_lon)

    tvar = nc.createVariable("time", merged_time.dtype, ("time",))
    tvar[:] = merged_time
    for k, v in time_attrs.items():
        try:
            setattr(tvar, k, v)
        except Exception:
            pass

    latvar = nc.createVariable("lat", "f4", ("lat",))
    lonvar = nc.createVariable("lon", "f4", ("lon",))
    latvar[:] = lat_vals
    lonvar[:] = lon_vals
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

    chunksizes = (min(chunk_time, total_hours), n_lat, n_lon)
    cf_var = nc.createVariable(
        var_name,
        "f4",
        ("time", "lat", "lon"),
        zlib=True,
        complevel=compress_level,
        chunksizes=chunksizes,
        fill_value=np.float32(np.nan),
    )
    cf_var.long_name = f"Multi-year mean hourly {cf_type} capacity factor ({start_year}-{end_year})"
    cf_var.units = "1"

    nc.source = "ERA5-Land"
    nc.start_year = start_year
    nc.end_year = end_year
    nc.total_hours = total_hours

    try:
        for info in tqdm(month_info, desc=f"merge-{cf_type}", unit="月"):
            ds = xr.open_dataset(info["file"], decode_times=False)
            vname = get_var_name(ds, var_name)
            n_hours = info["n_hours"]
            off = info["offset"]

            for t_start in range(0, n_hours, chunk_time):
                t_end = min(t_start + chunk_time, n_hours)
                chunk = ds[vname].isel(time=slice(t_start, t_end)).values
                cf_var[off + t_start : off + t_end, :, :] = chunk
                del chunk

            ds.close()
            gc.collect()
    finally:
        nc.close()

    os.replace(tmp_file, output_file)
    logger.info(f"✓ 已保存合并文件：{output_file}")
    return output_file


# ─────────────────────────────────────────────────────────────────
# 3. 主入口
# ─────────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(description="多年平均逐小时容量因子计算（ERA5-Land 风光共用）")
    parser.add_argument("--cf_dir", required=True, help="容量因子数据目录，含 {solar|wind}_cf_YYYY_MM.nc")
    parser.add_argument(
        "--type",
        required=True,
        choices=["solar", "wind"],
        help="容量因子类型：solar 或 wind",
    )
    parser.add_argument("--start_year", type=int, default=2015, help="起始年份（含），默认 2015")
    parser.add_argument("--end_year", type=int, default=2025, help="结束年份（含），默认 2025")
    parser.add_argument("--save_dir", required=True, help="保存根目录")
    parser.add_argument("--no_merge", action="store_true", help="不合并逐月文件（默认合并）")
    parser.add_argument("--chunk_time", type=int, default=24, help="每块处理的小时数，默认 24")
    parser.add_argument("--compress_level", type=int, default=4, help="NetCDF 压缩级别 0-9，默认 4")
    parser.add_argument("--overwrite", action="store_true", help="覆盖已有文件")
    args = parser.parse_args()

    cf_type = args.type
    out_base = output_base_dir(args.save_dir, cf_type, args.start_year, args.end_year)
    monthly_dir = out_base / "monthly"

    logger.info(f"容量因子类型：{cf_type}")
    logger.info(f"年份范围：{args.start_year}-{args.end_year}")
    logger.info(f"输出目录：{out_base}")
    logger.info(f"是否合并：{not args.no_merge}")

    # Step 1: 逐月计算多年平均
    for month in range(1, 13):
        logger.info(f"\n{'=' * 60}")
        logger.info(f"处理月份 {month:02d}")
        logger.info(f"{'=' * 60}")
        compute_monthly_mean(
            cf_dir=args.cf_dir,
            cf_type=cf_type,
            start_year=args.start_year,
            end_year=args.end_year,
            month=month,
            output_file=monthly_dir / f"mean_{cf_type}_cf_{month:02d}.nc",
            chunk_time=args.chunk_time,
            overwrite=args.overwrite,
            compress_level=args.compress_level,
        )

    # Step 2: 合并
    if not args.no_merge:
        logger.info(f"\n{'=' * 60}")
        logger.info("合并 12 个月均值")
        logger.info(f"{'=' * 60}")
        merged_file = out_base / f"mean_{cf_type}_cf_{args.start_year}_{args.end_year}.nc"
        merge_monthly_means(
            monthly_dir=monthly_dir,
            cf_type=cf_type,
            start_year=args.start_year,
            end_year=args.end_year,
            output_file=merged_file,
            overwrite=args.overwrite,
            compress_level=args.compress_level,
            chunk_time=args.chunk_time,
        )

    logger.info("\n✓ 全部完成！")


if __name__ == "__main__":
    main()

"""
使用方法示例：

# 太阳能（2015-2025，合并）
python S02E03_Compute_Mean_Hourly_CF.py \
  --cf_dir output/CFs_of_solar_ERA5Land \
  --type solar \
  --save_dir output

# 风能（自定义年份，不合并）
python S02E03_Compute_Mean_Hourly_CF.py \
  --cf_dir output/CFs_of_wind_ERA5Land \
  --type wind \
  --start_year 2016 --end_year 2024 \
  --save_dir output \
  --no_merge
"""
