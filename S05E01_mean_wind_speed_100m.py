"""
ERA5-Land 多年平均 100 m 风速
==============================

从 ERA5-Land u10/v10 全球逐小时月度数据计算多年平均 100 m 风速。

计算方法：
  1. 合成 10 m 风速：ws10 = sqrt(u10² + v10²)
  2. 幂律外推至 100 m：ws100 = ws10 * (100/10)^(1/7)
  3. 所有时刻取均值

示例
----
python S05E01_mean_wind_speed_100m.py \
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

POWER_LAW_ALPHA = 1.0 / 7.0
HUB_HEIGHT = 100.0
REF_HEIGHT = 10.0


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
    extrap_ratio = np.float32((HUB_HEIGHT / REF_HEIGHT) ** POWER_LAW_ALPHA)
    logger.info(f"风速外推系数: {float(extrap_ratio):.4f} (100m/10m)^(1/7)")

    for year, month in tqdm(iter_year_months(start_year, end_year, months), desc="wind_speed_100m", unit="月"):
        u10_path = era5land_file(data_dir, "u10", year, month)
        v10_path = era5land_file(data_dir, "v10", year, month)
        if not u10_path.exists() or not v10_path.exists():
            logger.warning(f"文件缺失，跳过 {year}-{month:02d}")
            continue

        ds_u, tname, latname, lonname = open_era5land(u10_path, "u10")
        ds_v, *_ = open_era5land(v10_path, "v10")
        da_u = prepare_dataarray(ds_u, "u10", tname, latname, lonname)
        da_v = prepare_dataarray(ds_v, "v10", tname, latname, lonname)

        try:
            for i0 in range(0, n_lat, chunk_lat):
                i1 = min(i0 + chunk_lat, n_lat)
                u = da_u.isel({latname: slice(i0, i1)}).values.astype(np.float32)
                v = da_v.isel({latname: slice(i0, i1)}).values.astype(np.float32)
                ws10 = np.sqrt(u ** 2 + v ** 2)
                ws100 = ws10 * extrap_ratio
                valid = np.isfinite(ws100)
                ws100_clean = np.where(valid, ws100, 0.0).astype(np.float64)
                sum_arr[i0:i1] += np.sum(ws100_clean, axis=0)
                count_arr[i0:i1] += np.sum(valid, axis=0).astype(np.int32)
                del u, v, ws10, ws100, valid, ws100_clean
        finally:
            ds_u.close()
            ds_v.close()
            del da_u, da_v
            gc.collect()

    with np.errstate(invalid="ignore"):
        result = np.where(count_arr > 0, sum_arr / count_arr, np.nan)
    logger.info(f"风速统计：有效格点 {np.sum(count_arr > 0)}, 全局均值 {np.nanmean(result):.2f} m/s")
    return result.astype(np.float32)


def main() -> None:
    parser = argparse.ArgumentParser(description="ERA5-Land 多年平均 100m 风速")
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
    grid_lats, grid_lons = resolve_grid(args.data_dir, args.start_year, args.end_year, month_list, "u10")

    logger.info(f"网格：lat={len(grid_lats)}, lon={len(grid_lons)}")
    logger.info(f"年份：{args.start_year}-{args.end_year}")

    from pathlib import Path
    out_file = Path(out_dir) / "wind_speed_100m.nc"
    if out_file.exists() and not args.overwrite:
        logger.info(f"已存在，跳过：{out_file}")
        return

    result = compute(args.data_dir, args.start_year, args.end_year,
                     grid_lats, grid_lons, args.chunk_lat, month_list)
    save_result(out_file, result, grid_lats, grid_lons,
                "wind_speed_100m", "Multi-year mean 100m wind speed", "m s-1",
                {"source": "ERA5-Land", "start_year": args.start_year, "end_year": args.end_year,
                 "method": "power law V100=V10*(100/10)^(1/7)", "alpha": POWER_LAW_ALPHA},
                args.compress_level)


if __name__ == "__main__":
    main()
