"""
各国 & 全球「逐年平均容量因子」计算
=====================================

基于已生成的 BCSD 降尺度风/光容量因子（CF）NetCDF 文件，计算两类指标：

1. 各国逐年平均 CF：对每个 BCSD 降尺度国家 + china，按未来年份给出一个标量
   平均容量因子（先按年内时间求均值，再对陆地格点求空间均值）。
2. 全球逐年平均 CF：对「N 个降尺度国家 + china」共 N+1 个区域的逐年平均值，
   按各区域的陆地格点数加权求平均：
       global = Σ_r (mean_cf_r * n_land_cells_r) / Σ_r n_land_cells_r
   （等价于把所有区域的陆地格点汇到一起求平均，大国权重更大）。

风电、光伏分别独立计算。

输入文件结构
------------
- 各国：{cf_root}/CFs_of_{energy}/{model}/{region}/
          {energy}_CF_{region}_{model}_{scenario}_*_allmonths.nc
- china：{cf_root}/CFs_of_{energy}_china/{model}/
          {energy}_CF_china_{model}_{scenario}_*_allmonths.nc
  变量：solar_cf / wind_cf，维度 (time, lat, lon)。

海洋掩膜
--------
默认用 global_land_mask 把海洋格点置 NaN 后再求平均（与仓库 S01E0x /
utils/cal_mean.py 口径一致）；这样结果不依赖磁盘文件是否已掩膜海洋。
可用 --no-land-mask 关闭。

输出（CSV, long format）
------------------------
output/annual_mean_cf/{model}/
  ├── per_country_annual_cf_{model}.csv
  └── global_annual_cf_{model}.csv

运行示例
--------
# 默认：NESM3，三种情景，风+光
python S06_Annual_Mean_CF.py

# 只算 ssp245 的光伏
python S06_Annual_Mean_CF.py --model NESM3 --scenarios ssp245 --energy solar
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("annual_mean_cf")

VAR_NAME = {"solar": "solar_cf", "wind": "wind_cf"}

# per_country csv 的去重主键
COUNTRY_KEY = ["energy", "model", "scenario", "region", "year"]
GLOBAL_KEY = ["energy", "model", "scenario", "year"]


# --------------------------------------------------------------------------- #
# 工具函数
# --------------------------------------------------------------------------- #
def build_land_mask(lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    """根据 lat/lon 生成陆地掩膜 (n_lat, n_lon)，True=陆地。

    global_land_mask 要求经度 ∈ [-180, 180]、纬度 ∈ [-90, 90]，
    因此对 [0, 360) 经度做归一化，并夹紧纬度避免浮点越界。
    """
    from global_land_mask import globe  # 仅在需要掩膜时才导入

    lats = np.asarray(lats, dtype=np.float64)
    lons = np.asarray(lons, dtype=np.float64)

    lons_conv = ((lons + 180.0) % 360.0) - 180.0
    lats_clip = np.clip(lats, -90.0, 90.0)

    lon2d, lat2d = np.meshgrid(lons_conv, lats_clip)  # (n_lat, n_lon)
    return globe.is_land(lat2d, lon2d).astype(bool)


def parse_years(years: str | None) -> tuple[int, int] | None:
    if not years:
        return None
    y0, y1 = years.split("-")
    return int(y0), int(y1)


def _pick_file(cands: list[Path]) -> Path | None:
    """多匹配时优先 2015-2060，否则取第一个。"""
    if not cands:
        return None
    pref = [c for c in cands if "2015-2060" in c.name]
    return (pref or sorted(cands))[0]


def find_country_file(cf_root: Path, energy: str, model: str, region: str, scenario: str) -> Path | None:
    base = cf_root / f"CFs_of_{energy}" / model / region
    return _pick_file(sorted(base.glob(f"{energy}_CF_{region}_{model}_{scenario}_*_allmonths.nc")))


def find_china_file(cf_root: Path, energy: str, model: str, scenario: str) -> Path | None:
    base = cf_root / f"CFs_of_{energy}_china" / model
    return _pick_file(sorted(base.glob(f"{energy}_CF_china_{model}_{scenario}_*_allmonths.nc")))


def list_regions(cf_root: Path, energy: str, model: str) -> list[str]:
    """扫描 CFs_of_{energy}/{model}/ 下的国家子目录，再追加 china。"""
    base = cf_root / f"CFs_of_{energy}" / model
    if not base.is_dir():
        return ["china"]
    countries = sorted(p.name for p in base.iterdir() if p.is_dir())
    return countries + ["china"]


# --------------------------------------------------------------------------- #
# 核心计算
# --------------------------------------------------------------------------- #
def compute_region_annual(
    nc_path: Path,
    var_name: str,
    land_mask_enabled: bool,
    years_range: tuple[int, int] | None,
) -> list[dict]:
    """计算单个区域的逐年平均 CF。

    步骤：按年逐片读取 -> 年内时间平均（每格点）-> 海洋掩膜 -> 陆地格点空间平均。
    逐年切片读取以控制内存（china 单文件很大）。
    """
    ds = xr.open_dataset(nc_path)  # lazy
    try:
        da = ds[var_name]

        land_da = None
        if land_mask_enabled:
            land = build_land_mask(ds["lat"].values, ds["lon"].values)
            land_da = xr.DataArray(
                land, dims=("lat", "lon"), coords={"lat": ds["lat"], "lon": ds["lon"]}
            )

        years = np.unique(da["time"].dt.year.values)
        if years_range:
            y0, y1 = years_range
            years = years[(years >= y0) & (years <= y1)]

        records = []
        for y in years:
            sub = da.sel(time=str(int(y)))  # 部分字符串索引，选出该年所有时间步
            n_t = int(sub["time"].size)
            if n_t == 0:
                continue
            cell_mean = sub.mean("time", skipna=True)  # (lat, lon)，触发读取该年数据
            if land_da is not None:
                cell_mean = cell_mean.where(land_da)
            cell_mean = cell_mean.load()
            vals = cell_mean.values
            finite = np.isfinite(vals)
            n_cells = int(finite.sum())
            mean_cf = float(np.nanmean(vals)) if n_cells > 0 else np.nan
            records.append(
                {
                    "year": int(y),
                    "mean_cf": mean_cf,
                    "n_timesteps": n_t,
                    "n_land_cells": n_cells,
                }
            )
        return records
    finally:
        ds.close()


def merge_csv(out_path: Path, df_new: pd.DataFrame, keys: list[str], overwrite: bool) -> pd.DataFrame:
    """增量合并到 csv：按 keys 去重保留新值；overwrite 则整表重写。"""
    if out_path.exists() and not overwrite:
        df_old = pd.read_csv(out_path)
        df = pd.concat([df_old, df_new], ignore_index=True)
        df = df.drop_duplicates(subset=keys, keep="last")
    else:
        df = df_new
    df = df.sort_values(keys).reset_index(drop=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    return df


# --------------------------------------------------------------------------- #
# 主流程
# --------------------------------------------------------------------------- #
def run(args) -> None:
    cf_root = Path(args.cf_root)
    out_dir = Path(args.output_dir) / args.model
    years_range = parse_years(args.years)

    scenarios = [s.strip() for s in args.scenarios.split(",") if s.strip()]
    energies = ["solar", "wind"] if args.energy == "both" else [args.energy]

    logger.info("模式=%s  情景=%s  能源=%s  海洋掩膜=%s", args.model, scenarios, energies, args.land_mask)

    country_rows: list[dict] = []
    global_rows: list[dict] = []

    for energy in energies:
        var_name = VAR_NAME[energy]
        regions = list_regions(cf_root, energy, args.model)
        logger.info("[%s] 共 %d 个区域（含 china）", energy, len(regions))

        for scenario in scenarios:
            # 收集本 (energy, scenario) 下各区域的逐年值 (mean_cf, n_land_cells)，
            # 供全球「按陆地格点数加权」平均使用
            per_year_values: dict[int, list[tuple[float, int]]] = {}

            for region in regions:
                if region == "china":
                    nc_path = find_china_file(cf_root, energy, args.model, scenario)
                else:
                    nc_path = find_country_file(cf_root, energy, args.model, region, scenario)

                if nc_path is None or not nc_path.exists():
                    logger.warning("  跳过：未找到 %s / %s / %s 的文件", energy, scenario, region)
                    continue

                logger.info("  处理 %s | %s | %s", energy, scenario, region)
                recs = compute_region_annual(nc_path, var_name, args.land_mask, years_range)

                for r in recs:
                    country_rows.append(
                        {
                            "energy": energy,
                            "model": args.model,
                            "scenario": scenario,
                            "region": region,
                            "year": r["year"],
                            "mean_cf": r["mean_cf"],
                            "n_timesteps": r["n_timesteps"],
                            "n_land_cells": r["n_land_cells"],
                        }
                    )
                    if np.isfinite(r["mean_cf"]) and r["n_land_cells"] > 0:
                        per_year_values.setdefault(r["year"], []).append(
                            (r["mean_cf"], r["n_land_cells"])
                        )

            # 全球逐年平均 = 各区域逐年平均「按陆地格点数加权」的平均（含 china）
            #   global = Σ_r (mean_cf_r * n_cells_r) / Σ_r n_cells_r
            for year, pairs in sorted(per_year_values.items()):
                cfs = np.array([p[0] for p in pairs], dtype=float)
                weights = np.array([p[1] for p in pairs], dtype=float)
                global_rows.append(
                    {
                        "energy": energy,
                        "model": args.model,
                        "scenario": scenario,
                        "year": year,
                        "mean_cf": float(np.average(cfs, weights=weights)),
                        "n_regions": len(pairs),
                        "total_land_cells": int(weights.sum()),
                    }
                )

    if not country_rows:
        logger.error("没有计算出任何结果，请检查输入路径。")
        return

    df_country = pd.DataFrame(country_rows)
    df_global = pd.DataFrame(global_rows)

    country_csv = out_dir / f"per_country_annual_cf_{args.model}.csv"
    global_csv = out_dir / f"global_annual_cf_{args.model}.csv"

    df_country = merge_csv(country_csv, df_country, COUNTRY_KEY, args.overwrite)
    df_global = merge_csv(global_csv, df_global, GLOBAL_KEY, args.overwrite)

    logger.info("各国逐年平均 -> %s （%d 行）", country_csv, len(df_country))
    logger.info("全球逐年平均 -> %s （%d 行）", global_csv, len(df_global))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="各国 & 全球逐年平均容量因子计算")
    p.add_argument("--model", default="NESM3", help="气象模式，默认 NESM3")
    p.add_argument(
        "--scenarios",
        default="ssp126,ssp245,ssp585",
        help="未来情景，逗号分隔，可只给一个，默认 ssp126,ssp245,ssp585",
    )
    p.add_argument("--energy", choices=["solar", "wind", "both"], default="both", help="能源类型")
    p.add_argument("--cf-root", default="output", help="CF 输出根目录，默认 output")
    p.add_argument("--output-dir", default="output/annual_mean_cf", help="结果输出根目录")
    p.add_argument("--years", default="", help="年份过滤 Y0-Y1，空=全部")
    p.add_argument(
        "--land-mask",
        dest="land_mask",
        action="store_true",
        default=True,
        help="用 global_land_mask 把海洋置 NaN 后再平均（默认开）",
    )
    p.add_argument(
        "--no-land-mask",
        dest="land_mask",
        action="store_false",
        help="关闭海洋掩膜，对全部格点（含海洋）求平均",
    )
    p.add_argument("--overwrite", action="store_true", help="整表重写而非增量合并")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
