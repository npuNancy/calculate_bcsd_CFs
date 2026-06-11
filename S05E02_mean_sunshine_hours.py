"""
ERA5-Land 多年平均日均太阳辐照时数
====================================

从 ERA5-Land ssrd 全球逐小时月度数据计算多年平均日均太阳辐照时数。

计算方法：
  1. 将 ssrd 累计量转为逐小时增量：hourly = diff(ssrd)，负值置 0
  2. 逐小时增量 / 3600 → GHI (W m⁻²)
  3. GHI > 120 W m⁻² 计为 1 小时日照（WMO 推荐阈值）
  4. 按日求和得到每日日照时数
  5. 多年所有日取均值

示例
----
python S05E02_mean_sunshine_hours.py \
  --start_year 2015 --end_year 2025 \
  --output_dir output/era5land_climatology
"""

from __future__ import annotations

import argparse
import gc

import numpy as np
from tqdm import tqdm

from era5land_utils import (
    era5land_file,
    iter_year_months,
    logger,
    open_era5land,
    parse_months,
    prepare_dataarray,
    resolve_grid,
    save_result,
)

GHI_THRESHOLD = 120.0  # W m⁻²，WMO 推荐日照阈值


def compute(
    data_dir: str,
    start_year: int,
    end_year: int,
    lats: np.ndarray,
    lons: np.ndarray,
    chunk_lat: int,
    months: list[int],
) -> np.ndarray:
    n_lat, n_lon = len(lats), len(lons)
    sum_arr = np.zeros((n_lat, n_lon), dtype=np.float64)
    count_arr = np.zeros((n_lat, n_lon), dtype=np.int32)
    logger.info(f"日照判定阈值: GHI > {GHI_THRESHOLD} W/m2")

    for year, month in tqdm(iter_year_months(start_year, end_year, months), desc="sunshine_hours", unit="月"):
        ssrd_path = era5land_file(data_dir, "ssrd", year, month)
        if not ssrd_path.exists():
            logger.warning(f"ssrd 文件缺失，跳过 {year}-{month:02d}")
            continue

        ds, tname, latname, lonname = open_era5land(ssrd_path, "ssrd")
        da = prepare_dataarray(ds, "ssrd", tname, latname, lonname)

        times = ds[tname].values
        dates = np.array([np.datetime64(t, "D") for t in times])
        unique_days = np.unique(dates)

        try:
            for i0 in range(0, n_lat, chunk_lat):
                i1 = min(i0 + chunk_lat, n_lat)
                raw = da.isel({latname: slice(i0, i1)}).values.astype(np.float32)

                # 累计量 → 逐小时增量
                hourly_ssrd = np.diff(raw, axis=0, prepend=np.zeros_like(raw[:1]))
                hourly_ssrd = np.maximum(hourly_ssrd, 0.0)

                # J m⁻² → W m⁻²
                ghi = hourly_ssrd / 3600.0
                sunshine = (ghi > GHI_THRESHOLD).astype(np.float32)

                # 逐纬度行按日聚合，避免 chained boolean indexing 问题
                for day in unique_days:
                    day_mask = dates == day
                    if not np.any(day_mask):
                        continue
                    daily_sun = sunshine[day_mask, :, :]
                    daily_raw = raw[day_mask, :, :]
                    daily_sum = np.nansum(daily_sun, axis=0)
                    daily_valid = np.any(np.isfinite(daily_raw), axis=0)
                    valid_vals = daily_sum[daily_valid].astype(np.float64)
                    sl = np.s_[i0:i1]
                    tmp_sum = sum_arr[sl].copy()
                    tmp_cnt = count_arr[sl].copy()
                    tmp_sum[daily_valid] += valid_vals
                    tmp_cnt[daily_valid] += 1
                    sum_arr[sl] = tmp_sum
                    count_arr[sl] = tmp_cnt

                del raw, hourly_ssrd, ghi, sunshine
        finally:
            ds.close()
            del da
            gc.collect()

    with np.errstate(invalid="ignore"):
        result = np.where(count_arr > 0, sum_arr / count_arr, np.nan)
    logger.info(f"日照时数统计：有效格点 {np.sum(count_arr > 0)}, 全局均值 {np.nanmean(result):.2f} h/day")
    return result.astype(np.float32)


def main() -> None:
    parser = argparse.ArgumentParser(description="ERA5-Land 多年平均日均太阳辐照时数")
    parser.add_argument("--data_dir", default="data/ERA5_land/global")
    parser.add_argument("--start_year", type=int, default=2015)
    parser.add_argument("--end_year", type=int, default=2025)
    parser.add_argument("--months", default="", help="逗号分隔月份；默认全年")
    parser.add_argument("--output_dir", default="output")
    parser.add_argument("--chunk_lat", type=int, default=100)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--compress_level", type=int, default=4)
    args = parser.parse_args()

    month_list = parse_months(args.months)
    out_dir = args.output_dir
    grid_lats, grid_lons = resolve_grid(args.data_dir, args.start_year, args.end_year, month_list, "ssrd")

    logger.info(f"网格：lat={len(grid_lats)}, lon={len(grid_lons)}")
    logger.info(f"年份：{args.start_year}-{args.end_year}")

    from pathlib import Path
    out_file = Path(out_dir) / "solar_radiation_hours.nc"
    if out_file.exists() and not args.overwrite:
        logger.info(f"已存在，跳过：{out_file}")
        return

    result = compute(args.data_dir, args.start_year, args.end_year,
                     grid_lats, grid_lons, args.chunk_lat, month_list)
    save_result(out_file, result, grid_lats, grid_lons,
                "sunshine_hours", "Multi-year mean daily sunshine hours", "h day-1",
                {"source": "ERA5-Land", "start_year": args.start_year, "end_year": args.end_year,
                 "method": f"hourly_ssrd = diff(cumulative); GHI = hourly_ssrd/3600; sunshine when GHI > {GHI_THRESHOLD} W/m2; daily sum then mean"},
                args.compress_level)


if __name__ == "__main__":
    main()
