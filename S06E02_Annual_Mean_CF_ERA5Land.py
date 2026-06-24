"""
各国 & 全球「逐年平均容量因子」计算（ERA5-Land 历史再分析版）
================================================================

与 ``S06E01_Annual_Mean_CF.py`` 口径完全对齐，但输入换成基于 **ERA5-Land 再分析**
算得的历史风/光容量因子（CF），给出两类逐年平均指标：

1. 各国逐年平均 CF：对 ``data/grid_of_regions`` 下每个区域网格界定的国家，按年给出
   一个标量平均容量因子（先按年内时间求均值，再对陆地格点求空间均值）。
2. 全球逐年平均 CF：对所有区域的逐年平均值，按各区域的陆地格点数加权求平均：
       global = Σ_r (mean_cf_r * n_land_cells_r) / Σ_r n_land_cells_r

风电、光伏分别独立计算。无 model / scenario 维度（再分析只有一套数据），
用常量列 ``source=ERA5Land`` 标识。

输入文件结构
------------
- ERA5-Land CF（全球 0.1° 逐小时，按月分文件，单文件 ~19GB）：
    {cf_root}/CFs_of_{energy}_ERA5Land/{energy}_cf_{year}_{month:02d}.nc
  变量：solar_cf / wind_cf，维度 (time, lat, lon)；lat 90→-90、lon [0,360)。
  海洋格点为 NaN（ERA5-Land 仅陆地有值）。
- 区域网格（关键：保证与各国 BCSD 区域一致）：
    {grid_dir}/{region}_grid.nc      （仅含 lat/lon 坐标，无数据变量）
  支持两类网格：
    * 规则 1D lat/lon（如 china、France 等）：正交索引截取，子集 (rlat, rlon)。
    * 曲线/旋转极点 2D lat(rlat,rlon)/lon(rlat,rlon)（如 NAM-12）：展平后逐点最近邻索引，子集 (N,)。
  注意：NAM-12 是北美域，与 México 等国在空间上重叠；二者同时计入全球加权时会重复计数，
  如需避免可用 --exclude-regions 排除其一（详见 --exclude-regions 说明）。

为什么读取区域网格
------------------
为保证 ERA5-Land 各国区域与 S06E01（BCSD 降尺度）一致、可直接对比，必须用
``grid_of_regions`` 界定每个国家的空间范围，再据此从全球 ERA5-Land 场中截取。

I/O 策略（重要）
----------------
ERA5-Land 文件的 HDF5 分块为 (24, 1801, 3600)，即每个 chunk 跨整个全球网格。
任何空间子集都必须解压整块，因此 **逐区域 .sel 会把每个 19GB 文件重复解压 N 次**。
本脚本改为：每个月文件 **只做一次单程扫描**，沿时间维按 chunk 大小（24 步）分块读取，
每读入一块全球场即为所有区域做最近邻索引截取并累加 sum/count，从而一次解压服务所有区域。

海洋掩膜
--------
默认用 global_land_mask 把海洋格点置 NaN 后再求平均（与 S06E01 / utils/cal_mean.py 一致）。
ERA5-Land 海洋本就是 NaN，此步多为冗余，仅为口径一致；可用 --no-land-mask 关闭。

输出（CSV, long format）
------------------------
output/annual_mean_cf_ERA5Land/
  ├── per_country_annual_cf_ERA5Land.csv
  └── global_annual_cf_ERA5Land.csv

运行示例
--------
# 默认：风+光，2015-2025
python S06E02_Annual_Mean_CF_ERA5Land.py

# 只算光伏、2015-2020
python S06E02_Annual_Mean_CF_ERA5Land.py --energy solar --years 2015-2020
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
logger = logging.getLogger("annual_mean_cf_era5land")

VAR_NAME = {"solar": "solar_cf", "wind": "wind_cf"}
SOURCE = "ERA5Land"

# 沿时间维分块读取的步长，对齐 ERA5-Land HDF5 chunk (24, 1801, 3600)
TIME_BLOCK = 24

# csv 去重主键
COUNTRY_KEY = ["energy", "source", "region", "year"]
GLOBAL_KEY = ["energy", "source", "year"]


# --------------------------------------------------------------------------- #
# 工具函数
# --------------------------------------------------------------------------- #
def _to_globe_inputs(lats: np.ndarray, lons: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """global_land_mask 要求经度 ∈ [-180, 180]、纬度 ∈ [-90, 90]，做归一化。"""
    lats = np.asarray(lats, dtype=np.float64)
    lons = np.asarray(lons, dtype=np.float64)
    lons_conv = ((lons + 180.0) % 360.0) - 180.0
    lats_clip = np.clip(lats, -90.0, 90.0)
    return lats_clip, lons_conv


def build_land_mask(lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    """根据 1D lat/lon 生成规则网格陆地掩膜 (n_lat, n_lon)，True=陆地。"""
    from global_land_mask import globe  # 仅在需要掩膜时才导入

    lats_clip, lons_conv = _to_globe_inputs(lats, lons)
    lon2d, lat2d = np.meshgrid(lons_conv, lats_clip)  # (n_lat, n_lon)
    return globe.is_land(lat2d, lon2d).astype(bool)


def build_land_mask_points(lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    """逐点（曲线网格 / 散点）陆地掩膜，lats/lons 同形，返回同形 bool。"""
    from global_land_mask import globe

    lats_clip, lons_conv = _to_globe_inputs(lats, lons)
    return globe.is_land(lats_clip, lons_conv).astype(bool)


def parse_years(years: str | None) -> tuple[int, int] | None:
    if not years:
        return None
    y0, y1 = years.split("-")
    return int(y0), int(y1)


def nearest_indices(grid: np.ndarray, target: np.ndarray) -> np.ndarray:
    """对 target 中每个值，返回其在 grid 中最近邻的整数索引。

    grid 可乱序（ERA5-Land lat 为降序）；通过对 grid 排序后 searchsorted 实现。
    """
    grid = np.asarray(grid, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    order = np.argsort(grid)
    g_sorted = grid[order]
    pos = np.searchsorted(g_sorted, target)
    pos = np.clip(pos, 1, len(g_sorted) - 1)
    left = g_sorted[pos - 1]
    right = g_sorted[pos]
    choose = np.where(np.abs(target - left) <= np.abs(target - right), pos - 1, pos)
    return order[choose]


# --------------------------------------------------------------------------- #
# 区域网格
# --------------------------------------------------------------------------- #
class Region:
    """缓存一个区域的截取索引与陆地掩膜。

    支持两类网格：
    - 规则 1D 网格（pointwise=False）：lat_idx (rlat,) 与 lon_idx (rlon,) 做**正交索引**，
      子集形状 (rlat, rlon)。
    - 曲线/旋转极点 2D 网格（pointwise=True，如 NAM-12）：lat_idx 与 lon_idx 等长 (N,)，
      做**逐点索引**，子集形状 (N,)（N = 原 2D 网格展平后的格点数）。
    """

    def __init__(
        self,
        name: str,
        lat_idx: np.ndarray,
        lon_idx: np.ndarray,
        land: np.ndarray | None,
        pointwise: bool = False,
    ):
        self.name = name
        self.lat_idx = lat_idx          # 行索引：规则为 (rlat,)，逐点为 (N,)
        self.lon_idx = lon_idx          # 列索引：规则为 (rlon,)，逐点为 (N,)
        self.land = land                # 与 sum 同形的 bool，或 None
        self.pointwise = pointwise
        shape = (lat_idx.size,) if pointwise else (lat_idx.size, lon_idx.size)
        self.sum = np.zeros(shape, dtype=np.float64)   # 逐年累加：每格点 CF 之和
        self.cnt = np.zeros(shape, dtype=np.int64)     # 逐年累加：每格点有限值计数
        self.n_t = 0                                   # 逐年累加：时间步数

    def reset(self) -> None:
        self.sum.fill(0.0)
        self.cnt.fill(0)
        self.n_t = 0

    def accumulate(self, block: np.ndarray) -> None:
        """累加一个时间块 (bt, 1801, 3600) 中本区域的 CF。"""
        if self.pointwise:
            sub = block[:, self.lat_idx, self.lon_idx]        # (bt, N) 逐点
        else:
            sub = block[:, self.lat_idx][:, :, self.lon_idx]  # (bt, rlat, rlon) 正交
        finite = np.isfinite(sub)
        self.sum += np.where(finite, sub, 0.0).sum(axis=0)
        self.cnt += finite.sum(axis=0)
        self.n_t += sub.shape[0]

    def finalize(self) -> tuple[float, int]:
        """年内累加结束后求该区域逐年平均 CF 与有效陆地格点数。"""
        with np.errstate(invalid="ignore", divide="ignore"):
            cell_mean = np.where(self.cnt > 0, self.sum / self.cnt, np.nan)
        if self.land is not None:
            cell_mean = np.where(self.land, cell_mean, np.nan)
        finite = np.isfinite(cell_mean)
        n_cells = int(finite.sum())
        mean_cf = float(np.nanmean(cell_mean)) if n_cells > 0 else np.nan
        return mean_cf, n_cells


def load_regions(
    grid_dir: Path,
    global_lat: np.ndarray,
    global_lon: np.ndarray,
    land_mask_enabled: bool,
    only: set[str] | None,
    exclude: set[str],
) -> list[Region]:
    """扫描 grid_dir 下 *_grid.nc，构建可用区域。

    - 规则 1D lat/lon（如 china、France）：正交索引截取。
    - 曲线/旋转极点 2D lat(rlat,rlon)/lon(rlat,rlon)（如 NAM-12）：展平后逐点最近邻索引截取。
    """
    regions: list[Region] = []
    for path in sorted(grid_dir.glob("*_grid.nc")):
        name = path.name[: -len("_grid.nc")]
        if only is not None and name not in only:
            continue
        if name in exclude:
            logger.info("  跳过区域（--exclude-regions）：%s", name)
            continue

        g = xr.open_dataset(path)
        try:
            if "lat" not in g.coords or "lon" not in g.coords:
                logger.warning("  跳过区域（无 lat/lon 坐标）：%s", name)
                continue
            lat_ndim, lon_ndim = g["lat"].ndim, g["lon"].ndim
            glat = g["lat"].values.astype(np.float64)
            glon = g["lon"].values.astype(np.float64)
        finally:
            g.close()

        glon = glon % 360.0  # 归一化到 [0,360) 对齐 ERA5-Land

        if lat_ndim == 1 and lon_ndim == 1:
            # 规则网格：正交索引
            lat_idx = nearest_indices(global_lat, glat)
            lon_idx = nearest_indices(global_lon, glon)
            land = build_land_mask(glat, glon) if land_mask_enabled else None
            regions.append(Region(name, lat_idx, lon_idx, land, pointwise=False))
        elif lat_ndim == 2 and lon_ndim == 2 and glat.shape == glon.shape:
            # 曲线/旋转极点 2D 网格：展平逐点最近邻
            flat_lat = glat.ravel()
            flat_lon = glon.ravel()
            lat_idx = nearest_indices(global_lat, flat_lat)
            lon_idx = nearest_indices(global_lon, flat_lon)
            land = build_land_mask_points(flat_lat, flat_lon) if land_mask_enabled else None
            logger.info("  曲线网格区域 %s：%s 个格点（展平逐点截取）", name, flat_lat.size)
            regions.append(Region(name, lat_idx, lon_idx, land, pointwise=True))
        else:
            logger.warning(
                "  跳过区域（lat/lon 维度不受支持：lat=%dD lon=%dD）：%s", lat_ndim, lon_ndim, name
            )
            continue
    return regions


# --------------------------------------------------------------------------- #
# 文件发现
# --------------------------------------------------------------------------- #
def month_file(cf_root: Path, energy: str, year: int, month: int) -> Path:
    return cf_root / f"CFs_of_{energy}_ERA5Land" / f"{energy}_cf_{year}_{month:02d}.nc"


def read_global_grid(cf_root: Path, energy: str) -> tuple[np.ndarray, np.ndarray]:
    """从任一存在的月文件读取全球 lat/lon（各文件一致）。"""
    base = cf_root / f"CFs_of_{energy}_ERA5Land"
    cands = sorted(base.glob(f"{energy}_cf_*.nc"))
    if not cands:
        raise FileNotFoundError(f"未找到 ERA5-Land CF 文件：{base}")
    ds = xr.open_dataset(cands[0])
    try:
        return ds["lat"].values.astype(np.float64), ds["lon"].values.astype(np.float64)
    finally:
        ds.close()


# --------------------------------------------------------------------------- #
# 核心计算
# --------------------------------------------------------------------------- #
def accumulate_month(nc_path: Path, var_name: str, regions: list[Region]) -> int:
    """单程扫描一个月文件，沿时间维按 TIME_BLOCK 分块为所有区域累加。返回时间步数。"""
    ds = xr.open_dataset(nc_path)
    try:
        da = ds[var_name]
        n_time = int(da["time"].size)
        for t0 in range(0, n_time, TIME_BLOCK):
            t1 = min(t0 + TIME_BLOCK, n_time)
            block = da.isel(time=slice(t0, t1)).values  # (bt, nlat, nlon)，触发解压
            for region in regions:
                region.accumulate(block)
        return n_time
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
    grid_dir = Path(args.grid_dir)
    out_dir = Path(args.output_dir)
    years_range = parse_years(args.years)

    energies = ["solar", "wind"] if args.energy == "both" else [args.energy]
    only = {s.strip() for s in args.regions.split(",") if s.strip()} or None
    exclude = {s.strip() for s in args.exclude_regions.split(",") if s.strip()}

    logger.info("能源=%s  年份=%s  海洋掩膜=%s", energies, args.years or "全部", args.land_mask)

    country_rows: list[dict] = []
    global_rows: list[dict] = []

    for energy in energies:
        var_name = VAR_NAME[energy]
        global_lat, global_lon = read_global_grid(cf_root, energy)
        regions = load_regions(grid_dir, global_lat, global_lon, args.land_mask, only, exclude)
        if not regions:
            logger.warning("[%s] 没有可用区域，跳过。", energy)
            continue
        logger.info("[%s] 共 %d 个区域", energy, len(regions))

        # 决定要计算的年份：扫描该能源下存在的文件年份，再按 --years 过滤
        base = cf_root / f"CFs_of_{energy}_ERA5Land"
        avail_years = sorted({int(p.name.split("_")[2]) for p in base.glob(f"{energy}_cf_*.nc")})
        if years_range:
            y0, y1 = years_range
            avail_years = [y for y in avail_years if y0 <= y <= y1]
        logger.info("[%s] 计算年份：%s", energy, avail_years)

        for year in avail_years:
            for region in regions:
                region.reset()

            months_done = 0
            for month in range(1, 13):
                nc_path = month_file(cf_root, energy, year, month)
                if not nc_path.exists():
                    continue
                logger.info("  扫描 %s | %d-%02d", energy, year, month)
                accumulate_month(nc_path, var_name, regions)
                months_done += 1

            if months_done == 0:
                logger.warning("  %d 年无任何月文件，跳过。", year)
                continue

            # 收集本年各区域逐年值，供全球「按陆地格点数加权」平均
            per_year_pairs: list[tuple[float, int]] = []
            for region in regions:
                mean_cf, n_cells = region.finalize()
                country_rows.append(
                    {
                        "energy": energy,
                        "source": SOURCE,
                        "region": region.name,
                        "year": int(year),
                        "mean_cf": mean_cf,
                        "n_timesteps": region.n_t,
                        "n_land_cells": n_cells,
                    }
                )
                if np.isfinite(mean_cf) and n_cells > 0:
                    per_year_pairs.append((mean_cf, n_cells))

            if per_year_pairs:
                cfs = np.array([p[0] for p in per_year_pairs], dtype=float)
                weights = np.array([p[1] for p in per_year_pairs], dtype=float)
                global_rows.append(
                    {
                        "energy": energy,
                        "source": SOURCE,
                        "year": int(year),
                        "mean_cf": float(np.average(cfs, weights=weights)),
                        "n_regions": len(per_year_pairs),
                        "total_land_cells": int(weights.sum()),
                    }
                )
            logger.info("  %d 年完成（%d 月，%d 区域）", year, months_done, len(per_year_pairs))

    if not country_rows:
        logger.error("没有计算出任何结果，请检查输入路径。")
        return

    df_country = pd.DataFrame(country_rows)
    df_global = pd.DataFrame(global_rows)

    country_csv = out_dir / f"per_country_annual_cf_{SOURCE}.csv"
    global_csv = out_dir / f"global_annual_cf_{SOURCE}.csv"

    df_country = merge_csv(country_csv, df_country, COUNTRY_KEY, args.overwrite)
    df_global = merge_csv(global_csv, df_global, GLOBAL_KEY, args.overwrite)

    logger.info("各国逐年平均 -> %s （%d 行）", country_csv, len(df_country))
    logger.info("全球逐年平均 -> %s （%d 行）", global_csv, len(df_global))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="各国 & 全球逐年平均容量因子计算（ERA5-Land 版）")
    p.add_argument("--energy", choices=["solar", "wind", "both"], default="both", help="能源类型")
    p.add_argument("--cf-root", default="output", help="ERA5-Land CF 根目录，默认 output")
    p.add_argument("--grid-dir", default="data/grid_of_regions", help="区域网格目录")
    p.add_argument("--output-dir", default="output/annual_mean_cf_ERA5Land", help="结果输出根目录")
    p.add_argument("--years", default="2015-2025", help="年份过滤 Y0-Y1，空=全部，默认 2015-2025")
    p.add_argument("--regions", default="", help="仅算指定区域（逗号分隔），空=全部")
    p.add_argument(
        "--exclude-regions",
        default="",
        help="排除区域（逗号分隔）。NAM-12 与 México 等空间重叠，"
        "若都计入全球加权会重复计数，可据需排除其一",
    )
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
        help="关闭海洋掩膜，对全部格点求平均",
    )
    p.add_argument("--overwrite", action="store_true", help="整表重写而非增量合并")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
