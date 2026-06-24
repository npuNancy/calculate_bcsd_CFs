#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""保存各个区域的网格 nc 文件。

只保存网格（经纬度坐标），不保存任何数据变量与时间维度。
结果输出到 data/grid_of_regions/<Region>_grid.nc。

区域来源：
  类型1 BCSD 各国：从 data/bcsd_outputs/NESM3/<Country>/NESM3 下任取一个 nc 读取网格。
  类型2 china：从 output/CFs_of_solar_china/NESM3 下的 CF nc 读取网格。
  类型3 NAM-12：已预先保存好，跳过。
"""

import glob
import os

import xarray as xr

# 项目根目录（utils 的上一级）
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

BCSD_DIR = os.path.join(ROOT, "data", "bcsd_outputs", "NESM3")
CHINA_NC = os.path.join(
    ROOT,
    "output",
    "CFs_of_solar_china",
    "NESM3",
    "solar_CF_china_NESM3_ssp126_2015-2060_allmonths.nc",
)
OUT_DIR = os.path.join(ROOT, "data", "grid_of_regions")

# 只保留网格相关坐标，丢弃时间维度
GRID_COORDS = ["lat", "lon"]


def extract_grid(src_nc):
    """从源 nc 中提取网格（经纬度坐标），返回一个不含数据变量与时间的 Dataset。"""
    with xr.open_dataset(src_nc) as ds:
        coords = {c: ds[c] for c in GRID_COORDS if c in ds.coords or c in ds.variables}
        grid = xr.Dataset(coords=coords)
        # 保留网格相关的全局属性中可能有用的信息
        grid.load()
    return grid


def save_grid(region, src_nc):
    out_path = os.path.join(OUT_DIR, f"{region}_grid.nc")
    grid = extract_grid(src_nc)
    grid.attrs["region"] = region
    grid.attrs["grid_source"] = os.path.relpath(src_nc, ROOT)
    grid.to_netcdf(out_path)
    sizes = dict(grid.sizes)
    print(f"  [OK] {region:20s} -> {os.path.relpath(out_path, ROOT)}  {sizes}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # 类型1：BCSD 各国
    print("=== 类型1：BCSD 各国 ===")
    countries = sorted(
        d
        for d in os.listdir(BCSD_DIR)
        if os.path.isdir(os.path.join(BCSD_DIR, d))
    )
    for country in countries:
        country_dir = os.path.join(BCSD_DIR, country, "NESM3")
        nc_files = sorted(glob.glob(os.path.join(country_dir, "*.nc")))
        if not nc_files:
            print(f"  [SKIP] {country}: 未找到 nc 文件 ({country_dir})")
            continue
        # 都读取同一个文件（排序后的第一个）以保持一致
        save_grid(country, nc_files[0])

    # 类型2：china
    print("=== 类型2：china ===")
    if os.path.exists(CHINA_NC):
        save_grid("china", CHINA_NC)
    else:
        print(f"  [SKIP] china: 未找到 {CHINA_NC}")

    # 类型3：NAM-12 已保存好，跳过
    print("=== 类型3：NAM-12 已存在，跳过 ===")

    print("完成。输出目录：", os.path.relpath(OUT_DIR, ROOT))


if __name__ == "__main__":
    main()
