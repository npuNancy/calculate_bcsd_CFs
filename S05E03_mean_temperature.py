"""
ERA5-Land 多年平均地表空气温度
================================

从 ERA5-Land t2m 全球逐小时月度数据计算多年平均 2 m 气温。

计算方法：
  1. t2m (K) - 273.15 → °C
  2. 所有时刻取均值

示例
----
python S05E03_mean_temperature.py \
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

    for year, month in tqdm(iter_year_months(start_year, end_year, months), desc="mean_temperature", unit="月"):
        t2m_path = era5land_file(data_dir, "t2m", year, month)
        if not t2m_path.exists():
            logger.warning(f"t2m 文件缺失，跳过 {year}-{month:02d}")
            continue

        ds, tname, latname, lonname = open_era5land(t2m_path, "t2m")
        da = prepare_dataarray(ds, "t2m", tname, latname, lonname)

        try:
            for i0 in range(0, n_lat, chunk_lat):
                i1 = min(i0 + chunk_lat, n_lat)
                data = da.isel({latname: slice(i0, i1)}).values.astype(np.float32)
                t_c = data - 273.15
                valid = np.isfinite(t_c)
                t_c_clean = np.where(valid, t_c, 0.0).astype(np.float64)
                sum_arr[i0:i1] += np.sum(t_c_clean, axis=0)
                count_arr[i0:i1] += np.sum(valid, axis=0).astype(np.int32)
                del data, t_c, valid, t_c_clean
        finally:
            ds.close()
            del da
            gc.collect()

    with np.errstate(invalid="ignore"):
        result = np.where(count_arr > 0, sum_arr / count_arr, np.nan)
    logger.info(f"温度统计：有效格点 {np.sum(count_arr > 0)}, 全局均值 {np.nanmean(result):.2f} °C")
    return result.astype(np.float32)


def main() -> None:
    parser = argparse.ArgumentParser(description="ERA5-Land 多年平均地表空气温度")
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
    grid_lats, grid_lons = resolve_grid(args.data_dir, args.start_year, args.end_year, month_list, "t2m")

    logger.info(f"网格：lat={len(grid_lats)}, lon={len(grid_lons)}")
    logger.info(f"年份：{args.start_year}-{args.end_year}")

    from pathlib import Path
    out_file = Path(out_dir) / "surface_air_temperature.nc"
    if out_file.exists() and not args.overwrite:
        logger.info(f"已存在，跳过：{out_file}")
        return

    result = compute(args.data_dir, args.start_year, args.end_year,
                     grid_lats, grid_lons, args.chunk_lat, month_list)
    save_result(out_file, result, grid_lats, grid_lons,
                "t2m_mean", "Multi-year mean 2m air temperature", "degC",
                {"source": "ERA5-Land", "start_year": args.start_year, "end_year": args.end_year,
                 "method": "t2m (K) - 273.15 -> degC; mean over all time steps"},
                args.compress_level)


if __name__ == "__main__":
    main()
