#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
使用 Cartopy 绘制风电/光伏容量因子年度平均图。

适用目录结构示例：

    output/CFs_of_solar/MPI-ESM1-2-HR/
        Austria/solar_CF_Austria_MPI-ESM1-2-HR_ssp126_2015-2060_allmonths.nc
        Germany/solar_CF_Germany_MPI-ESM1-2-HR_ssp126_2015-2060_allmonths.nc
        ...

    output/CFs_of_wind/MPI-ESM1-2-HR/
        Austria/wind_CF_Austria_MPI-ESM1-2-HR_ssp126_2015-2060_allmonths.nc
        Germany/wind_CF_Germany_MPI-ESM1-2-HR_ssp126_2015-2060_allmonths.nc
        ...

核心做法：
  1. 对每个国家/地区的 nc 文件，筛选指定 year 的时间步；
  2. 计算该年平均容量因子；
  3. 构建全球 ERA5-Land 0.1° × 0.1° 空白网格；
  4. 将每个国家/地区的矩形网格填入全球网格；
  5. 对多个矩形重叠的格点，默认取平均，避免被后处理顺序影响；
  6. 用 Cartopy 绘图，图片文件名会根据技术、模型、情景、年份、区域、extent 等信息自动生成；
  7. 可选保存合并后的年度平均 nc 文件。

依赖：
    numpy xarray netCDF4 matplotlib cartopy

"""

from __future__ import annotations

import argparse
import logging
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

# 服务器/超算无显示器时也能直接保存图片。
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except Exception as exc:  # pragma: no cover
    raise ImportError(
        "本脚本需要 cartopy。建议使用 conda 安装：\n"
        "  conda install -c conda-forge cartopy xarray netCDF4 matplotlib\n"
        "或 pip 安装：\n"
        "  pip install -i https://mirrors.aliyun.com/pypi/simple cartopy xarray netCDF4 matplotlib"
    ) from exc


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ERA5-Land 0.1° global grid: latitude 90 -> -90, longitude 0 -> 359.9。
GLOBAL_RES = 0.1
GLOBAL_NLAT = 1801
GLOBAL_NLON = 3600
GLOBAL_LATS = np.round(90.0 - np.arange(GLOBAL_NLAT, dtype=np.float64) * GLOBAL_RES, 1).astype(np.float32)
GLOBAL_LONS_0360 = np.round(np.arange(GLOBAL_NLON, dtype=np.float64) * GLOBAL_RES, 1).astype(np.float32)


@dataclass
class RegionResult:
    region: str
    file: Path
    n_time: int
    n_valid: int
    mean_cf: float
    lat_min: float
    lat_max: float
    lon_min_0360: float
    lon_max_0360: float


def split_country_args(values: list[str] | None) -> list[str] | None:
    """同时支持 `--countries Austria Germany` 和 `--countries Austria,Germany`。"""
    if not values:
        return None
    out: list[str] = []
    for item in values:
        for part in item.split(","):
            part = part.strip()
            if part:
                out.append(part)
    return out or None


def is_valid_region_dir(path: Path) -> bool:
    """过滤临时目录，并确认目录内至少有 nc 文件。"""
    name = path.name
    if not path.is_dir():
        return False
    if not name or name.startswith("_") or name.endswith("_") or name.endswith("_repeated"):
        return False
    return any(path.glob("*.nc")) or any(path.glob("**/*.nc"))


def discover_regions(cf_root: Path, countries: list[str] | None) -> list[str]:
    """发现需要绘制的国家/地区目录。"""
    if countries:
        return countries

    regions = [p.name for p in sorted(cf_root.iterdir()) if is_valid_region_dir(p)]
    if regions:
        return regions

    # 兼容 cf_root 本身就是单个国家目录的情况。
    if any(cf_root.glob("*.nc")):
        return [cf_root.name]

    raise RuntimeError(f"没有在 {cf_root} 下发现包含 NetCDF 文件的国家/地区目录。")


def infer_model_name(cf_root: Path) -> str:
    """从 cf_root 最后一层路径推断模型名，例如 MPI-ESM1-2-HR。"""
    return cf_root.name


def normalize_technology(technology: str, cf_root: Path) -> str:
    """如果 technology=auto，则根据路径名粗略推断。"""
    tech = technology.lower().strip()
    if tech != "auto":
        if tech not in {"solar", "wind"}:
            raise ValueError("--technology 只能是 auto、solar 或 wind。")
        return tech

    root_l = str(cf_root).lower()
    if "solar" in root_l or "pv" in root_l:
        return "solar"
    if "wind" in root_l:
        return "wind"
    return "auto"


def candidate_patterns(region: str, scenario: str, technology: str) -> list[str]:
    """针对当前输出命名规则构造候选文件 pattern。"""
    if technology == "solar":
        return [
            f"solar_CF_{region}_*_{scenario}_*.nc",
            f"solar*{region}*{scenario}*.nc",
            f"*solar*{scenario}*.nc",
            f"*{scenario}*.nc",
        ]
    if technology == "wind":
        return [
            f"wind_CF_{region}_*_{scenario}_*.nc",
            f"wind*{region}*{scenario}*.nc",
            f"*wind*{scenario}*.nc",
            f"*{scenario}*.nc",
        ]
    return [f"*{region}*{scenario}*.nc", f"*{scenario}*.nc"]


def find_region_dir(cf_root: Path, region: str) -> Path:
    """返回国家/地区目录；兼容 cf_root 本身就是该区域目录。"""
    direct = cf_root / region
    if direct.exists():
        return direct
    if cf_root.name == region and any(cf_root.glob("*.nc")):
        return cf_root
    raise FileNotFoundError(f"找不到 region 目录：{direct}")


def find_coord_name(ds: xr.Dataset | xr.DataArray, candidates: Iterable[str]) -> str:
    """兼容 time/valid_time, lat/latitude, lon/longitude 等坐标名。"""
    for name in candidates:
        if name in ds.coords or name in ds.dims:
            return name
    raise KeyError(f"找不到坐标名，候选={list(candidates)}；dims={list(ds.dims)}")


def infer_cf_var(ds: xr.Dataset, technology: str) -> str:
    """识别容量因子变量名。"""
    preferred: list[str]
    if technology == "solar":
        preferred = ["solar_cf", "pv_cf", "cf"]
    elif technology == "wind":
        preferred = ["wind_cf", "onshore_wind_cf", "cf"]
    else:
        preferred = ["solar_cf", "wind_cf", "pv_cf", "onshore_wind_cf", "cf"]

    for name in preferred:
        if name in ds.data_vars:
            return name

    cf_like = [name for name in ds.data_vars if "cf" in name.lower()]
    if len(cf_like) == 1:
        return cf_like[0]
    if len(ds.data_vars) == 1:
        return list(ds.data_vars)[0]
    raise KeyError(f"无法自动识别容量因子变量，可用变量：{list(ds.data_vars)}")


def dataset_has_year(ds: xr.Dataset, year: int) -> tuple[bool, int]:
    """判断一个 nc 文件是否包含目标年份。"""
    time_name = find_coord_name(ds, ["time", "valid_time"])
    years = ds[time_name].dt.year.values
    n_time = int(np.count_nonzero(years == year))
    return n_time > 0, n_time


def choose_region_file(region_dir: Path, region: str, scenario: str, year: int, technology: str) -> Path:
    """在一个国家/地区目录下，选择包含目标年份和情景的容量因子文件。"""
    files: list[Path] = []
    for pat in candidate_patterns(region, scenario, technology):
        files = sorted(region_dir.glob(pat))
        if files:
            break

    if not files:
        raise FileNotFoundError(
            f"{region}: 找不到 scenario={scenario} 的 nc 文件。\n"
            f"搜索目录：{region_dir}\n"
            f"尝试 pattern：{candidate_patterns(region, scenario, technology)}"
        )

    valid: list[Path] = []
    for f in files:
        try:
            with xr.open_dataset(f, decode_times=True) as ds:
                ok, _ = dataset_has_year(ds, year)
                if ok:
                    # 顺便确认变量可识别。
                    infer_cf_var(ds, technology)
                    valid.append(f)
        except Exception as exc:
            logger.warning(f"跳过无法读取或不含目标年份的文件：{f}；原因：{exc}")

    if not valid:
        raise FileNotFoundError(f"{region}: 匹配到 {len(files)} 个 {scenario} 文件，但没有文件包含 year={year}。")

    if len(valid) > 1:
        logger.warning(f"{region}: 有多个文件包含 year={year}，将使用第一个：{valid[0]}")
    return valid[0]


def lon_to_0360(lon: np.ndarray) -> np.ndarray:
    """经度统一到 [0, 360)。"""
    out = np.mod(lon.astype(np.float64), 360.0)
    out = np.where(np.isclose(out, 360.0), 0.0, out)
    return out.astype(np.float64)


def coord_to_global_indices(lats: np.ndarray, lons: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """将 0.1° 网格中心坐标映射到全球数组索引。"""
    lat_idx = np.rint((90.0 - lats.astype(np.float64)) / GLOBAL_RES).astype(np.int64)
    lon_idx = np.rint(lon_to_0360(lons) / GLOBAL_RES).astype(np.int64) % GLOBAL_NLON

    if np.any(lat_idx < 0) or np.any(lat_idx >= GLOBAL_NLAT):
        raise ValueError(f"纬度超出全球 0.1° 网格范围：min={lats.min()}, max={lats.max()}")
    if np.any(lon_idx < 0) or np.any(lon_idx >= GLOBAL_NLON):
        raise ValueError(f"经度超出全球 0.1° 网格范围：min={lons.min()}, max={lons.max()}")

    # 防止输入不是 0.1° 网格中心。
    lat_back = GLOBAL_LATS[lat_idx].astype(np.float64)
    lon_back = GLOBAL_LONS_0360[lon_idx].astype(np.float64)
    lon_in = lon_to_0360(lons)
    lon_diff = np.abs(((lon_back - lon_in + 180.0) % 360.0) - 180.0)
    if np.nanmax(np.abs(lat_back - lats.astype(np.float64))) > 1e-3 or np.nanmax(lon_diff) > 1e-3:
        raise ValueError("输入经纬度看起来不是 ERA5-Land 0.1° 网格中心，无法安全填入全球 0.1° 网格。")

    return lat_idx, lon_idx


def read_year_mean_cf(nc_file: Path, year: int, technology: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, str]:
    """读取单个 nc 文件，并返回该年的平均容量因子二维数组。"""
    with xr.open_dataset(nc_file, decode_times=True) as ds:
        var_name = infer_cf_var(ds, technology)
        time_name = find_coord_name(ds, ["time", "valid_time"])
        lat_name = find_coord_name(ds, ["lat", "latitude"])
        lon_name = find_coord_name(ds, ["lon", "longitude"])

        da = ds[var_name]
        if time_name not in da.dims:
            raise ValueError(f"{nc_file}: 变量 {var_name} 不含时间维 {time_name}。")

        year_mask = ds[time_name].dt.year == year
        n_time = int(year_mask.sum().item())
        if n_time == 0:
            raise ValueError(f"{nc_file}: 没有 year={year} 的时间步。")

        da_year = da.where(year_mask, drop=True)
        da_mean = da_year.mean(dim=time_name, skipna=True)
        da_mean = da_mean.transpose(lat_name, lon_name)

        arr = da_mean.values.astype(np.float32)
        lats = ds[lat_name].values.astype(np.float64)
        lons = ds[lon_name].values.astype(np.float64)

    arr = np.where(np.isfinite(arr), arr, np.nan).astype(np.float32)
    return arr, lats, lons, n_time, var_name


def update_global_grid(
    global_sum: np.ndarray,
    global_count: np.ndarray,
    arr: np.ndarray,
    lat_idx: np.ndarray,
    lon_idx: np.ndarray,
    overlap_method: str,
) -> None:
    """将一个区域的年度平均 CF 填入全球网格。"""
    valid = np.isfinite(arr)

    if overlap_method == "mean":
        global_sum[np.ix_(lat_idx, lon_idx)] += np.where(valid, arr, 0.0).astype(np.float32)
        global_count[np.ix_(lat_idx, lon_idx)] += valid.astype(np.uint16)
        return

    sub_sum = global_sum[np.ix_(lat_idx, lon_idx)]
    sub_count = global_count[np.ix_(lat_idx, lon_idx)]

    if overlap_method == "first":
        mask = valid & (sub_count == 0)
        sub_sum[mask] = arr[mask]
        sub_count[mask] = 1
    elif overlap_method == "last":
        mask = valid
        sub_sum[mask] = arr[mask]
        sub_count[mask] = 1
    elif overlap_method == "max":
        mask_new = valid & (sub_count == 0)
        mask_old = valid & (sub_count > 0)
        sub_sum[mask_new] = arr[mask_new]
        sub_count[mask_new] = 1
        sub_sum[mask_old] = np.maximum(sub_sum[mask_old], arr[mask_old])
    else:
        raise ValueError("--overlap_method 只能是 mean、first、last 或 max。")

    global_sum[np.ix_(lat_idx, lon_idx)] = sub_sum
    global_count[np.ix_(lat_idx, lon_idx)] = sub_count


def merge_regions_to_global(
    cf_root: Path,
    regions: list[str],
    scenario: str,
    year: int,
    technology: str,
    overlap_method: str,
) -> tuple[np.ndarray, np.ndarray, list[RegionResult], str]:
    """读取多个国家/地区，并合并到全球 0.1° 网格。"""
    global_sum = np.zeros((GLOBAL_NLAT, GLOBAL_NLON), dtype=np.float32)
    global_count = np.zeros((GLOBAL_NLAT, GLOBAL_NLON), dtype=np.uint16)
    results: list[RegionResult] = []
    used_var_name = ""

    for i, region in enumerate(regions, 1):
        logger.info(f"[{i}/{len(regions)}] 处理 region={region}")
        region_dir = find_region_dir(cf_root, region)
        nc_file = choose_region_file(region_dir, region, scenario, year, technology)
        arr, lats, lons, n_time, var_name = read_year_mean_cf(nc_file, year, technology)
        used_var_name = used_var_name or var_name

        lat_idx, lon_idx = coord_to_global_indices(lats, lons)
        update_global_grid(global_sum, global_count, arr, lat_idx, lon_idx, overlap_method)

        finite = np.isfinite(arr)
        n_valid = int(np.count_nonzero(finite))
        mean_cf = float(np.nanmean(arr)) if n_valid else float("nan")
        lons0360 = lon_to_0360(lons)
        results.append(
            RegionResult(
                region=region,
                file=nc_file,
                n_time=n_time,
                n_valid=n_valid,
                mean_cf=mean_cf,
                lat_min=float(np.nanmin(lats)),
                lat_max=float(np.nanmax(lats)),
                lon_min_0360=float(np.nanmin(lons0360)),
                lon_max_0360=float(np.nanmax(lons0360)),
            )
        )
        logger.info(f"  file={nc_file.name}; time_steps={n_time}; valid_cells={n_valid}; mean_cf={mean_cf:.4f}")

    with np.errstate(invalid="ignore", divide="ignore"):
        if overlap_method == "mean":
            global_cf = np.where(global_count > 0, global_sum / global_count, np.nan).astype(np.float32)
        else:
            global_cf = np.where(global_count > 0, global_sum, np.nan).astype(np.float32)

    filled_cells = int(np.count_nonzero(global_count > 0))
    overlap_cells = int(np.count_nonzero(global_count > 1))
    max_overlap = int(global_count.max()) if filled_cells else 0
    logger.info(f"全球网格填充完成：filled_cells={filled_cells}, overlap_cells={overlap_cells}, max_overlap={max_overlap}")
    return global_cf, global_count, results, used_var_name


def lon0360_to_180(lon: np.ndarray) -> np.ndarray:
    """[0, 360) -> [-180, 180)。"""
    return (((lon.astype(np.float64) + 180.0) % 360.0) - 180.0).astype(np.float32)


def prepare_plot_data(
    global_cf: np.ndarray,
    lats: np.ndarray,
    lons0360: np.ndarray,
    extent_mode: str,
    pad_degree: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[float]]:
    """将全球 0-360 经度数据转成 -180~180，并按有效数据范围裁剪。"""
    lons180 = lon0360_to_180(lons0360)
    order = np.argsort(lons180)
    lons_plot = lons180[order]
    data_plot = global_cf[:, order]
    lats_plot = lats

    valid = np.isfinite(data_plot)
    if extent_mode == "world" or not np.any(valid):
        return data_plot, lats_plot, lons_plot, [-180.0, 180.0, -90.0, 90.0]

    row_any = np.any(valid, axis=1)
    col_any = np.any(valid, axis=0)
    row_idx = np.where(row_any)[0]
    col_idx = np.where(col_any)[0]

    r0, r1 = int(row_idx.min()), int(row_idx.max())
    c0, c1 = int(col_idx.min()), int(col_idx.max())

    # pad 换算成格点数。
    pad_n = int(np.ceil(pad_degree / GLOBAL_RES))
    r0 = max(0, r0 - pad_n)
    r1 = min(len(lats_plot) - 1, r1 + pad_n)
    c0 = max(0, c0 - pad_n)
    c1 = min(len(lons_plot) - 1, c1 + pad_n)

    data_crop = data_plot[r0 : r1 + 1, c0 : c1 + 1]
    lats_crop = lats_plot[r0 : r1 + 1]
    lons_crop = lons_plot[c0 : c1 + 1]

    lon_min = float(max(-180.0, np.nanmin(lons_crop)))
    lon_max = float(min(180.0, np.nanmax(lons_crop)))
    lat_min = float(max(-90.0, np.nanmin(lats_crop)))
    lat_max = float(min(90.0, np.nanmax(lats_crop)))
    return data_crop, lats_crop, lons_crop, [lon_min, lon_max, lat_min, lat_max]



def sanitize_filename_part(text: str, max_len: int = 120) -> str:
    """把任意字符串转换成适合文件名的一段。"""
    text = re.sub(r"[^A-Za-z0-9._+-]+", "-", text.strip())
    text = re.sub(r"-+", "-", text).strip("-._")
    if not text:
        text = "unknown"
    return text[:max_len]


def build_region_tag(regions: list[str], countries_was_specified: bool) -> str:
    """根据区域列表生成简短文件名标签。"""
    if not countries_was_specified:
        return f"all_regions_n{len(regions)}"
    if len(regions) <= 5:
        return sanitize_filename_part("-".join(regions), max_len=100)
    return f"selected_regions_n{len(regions)}"


def build_output_png_path(
    output_dir: Path,
    technology: str,
    model: str,
    scenario: str,
    year: int,
    regions: list[str],
    countries_was_specified: bool,
    extent: str,
    overlap_method: str,
) -> Path:
    """根据关键信息自动生成输出图片路径。"""
    tech = sanitize_filename_part(technology)
    model_tag = sanitize_filename_part(model)
    scenario_tag = sanitize_filename_part(scenario)
    region_tag = build_region_tag(regions, countries_was_specified)
    filename = f"{tech}_CF_{model_tag}_{scenario_tag}_{year}_{region_tag}_extent-{extent}_overlap-{overlap_method}.png"
    return output_dir / filename


def make_plot_cmap(cmap_name: str, zero_color: str) -> matplotlib.colors.Colormap:
    """构造绘图 colormap；NaN/被遮罩的 0 值用单独颜色显示。"""
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(zero_color)
    # cmap.set_bad((1, 1, 1, 0))  # transparent bad values
    return cmap


def mask_zero_for_plot(data: np.ndarray, zero_threshold: float | None) -> np.ndarray:
    """仅绘图时遮罩 0 或接近 0 的值，避免大面积 0 值占据色标最低端颜色。"""
    out = data.astype(np.float32, copy=True)
    if zero_threshold is None or zero_threshold < 0:
        return out
    out[np.isfinite(out) & (out <= zero_threshold)] = np.nan
    return out


def add_cartopy_features(ax, linewidth: float = 0.5) -> None:
    """添加地图要素；若 Natural Earth 数据下载失败，仍继续绘图。"""
    for feature, kwargs in [
        (cfeature.LAND, {"facecolor": "none"}),
        (cfeature.OCEAN, {"facecolor": "none"}),
        (cfeature.COASTLINE, {"linewidth": linewidth}),
        (cfeature.BORDERS, {"linewidth": linewidth * 0.8}),
    ]:
        try:
            ax.add_feature(feature, **kwargs)
        except Exception as exc:
            logger.warning(f"添加 Cartopy 地图要素失败，已跳过：{feature}; {exc}")


def plot_cf_map(
    global_cf: np.ndarray,
    technology: str,
    model: str,
    scenario: str,
    year: int,
    output_png: Path,
    title: str | None,
    extent_mode: str,
    pad_degree: float,
    cmap: str,
    zero_threshold: float | None,
    zero_color: str,
    vmin: float,
    vmax: float,
    dpi: int,
    figsize: tuple[float, float],
    add_features: bool,
) -> None:
    """绘制容量因子地图。"""
    global_cf_for_plot = mask_zero_for_plot(global_cf, zero_threshold)
    data_plot, lats_plot, lons_plot, extent = prepare_plot_data(
        global_cf=global_cf_for_plot,
        lats=GLOBAL_LATS,
        lons0360=GLOBAL_LONS_0360,
        extent_mode=extent_mode,
        pad_degree=pad_degree,
    )

    output_png.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    if add_features:
        add_cartopy_features(ax)

    mesh = ax.pcolormesh(
        lons_plot,
        lats_plot,
        data_plot,
        transform=ccrs.PlateCarree(),
        shading="auto",
        cmap=make_plot_cmap(cmap, zero_color),
        vmin=vmin,
        vmax=vmax,
    )

    cb = plt.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.055, shrink=0.82)
    cb.set_label("Annual mean capacity factor")

    try:
        gl = ax.gridlines(draw_labels=True, linewidth=0.3, linestyle="--", alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
    except Exception:
        pass

    tech_label = {"solar": "Solar PV", "wind": "Onshore wind"}.get(technology, "Capacity factor")
    if title is None:
        title = f"{tech_label} annual mean CF | {model} | {scenario} | {year}"
    ax.set_title(title)

    fig.savefig(output_png, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"✓ 已保存图片：{output_png}")


def save_merged_nc(
    save_nc: Path,
    global_cf: np.ndarray,
    global_count: np.ndarray,
    technology: str,
    model: str,
    scenario: str,
    year: int,
    overlap_method: str,
    regions: list[str],
) -> None:
    """保存合并后的年度平均全球网格。"""
    save_nc.parent.mkdir(parents=True, exist_ok=True)
    var_name = "solar_cf_yearmean" if technology == "solar" else "wind_cf_yearmean" if technology == "wind" else "cf_yearmean"
    ds = xr.Dataset(
        data_vars={
            var_name: (("latitude", "longitude"), global_cf.astype(np.float32)),
            "overlap_count": (("latitude", "longitude"), global_count.astype(np.uint16)),
        },
        coords={
            "latitude": ("latitude", GLOBAL_LATS.astype(np.float32)),
            "longitude": ("longitude", GLOBAL_LONS_0360.astype(np.float32)),
        },
        attrs={
            "description": "Merged annual mean capacity factor on global ERA5-Land 0.1 degree grid.",
            "technology": technology,
            "model": model,
            "scenario": scenario,
            "year": year,
            "overlap_method": overlap_method,
            "regions": ",".join(regions),
            "longitude_convention": "0-360 degrees_east",
        },
    )
    ds[var_name].attrs.update({"units": "1", "long_name": "Annual mean capacity factor"})
    ds["overlap_count"].attrs.update({"units": "1", "long_name": "Number of region rectangles contributing to each grid cell"})
    ds.to_netcdf(save_nc)
    logger.info(f"✓ 已保存合并 nc：{save_nc}")


def write_region_log(results: list[RegionResult], log_csv: Path | None) -> None:
    """保存每个国家/地区使用的文件和统计信息。"""
    if log_csv is None:
        return
    import csv

    log_csv.parent.mkdir(parents=True, exist_ok=True)
    with log_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "region",
                "file",
                "n_time",
                "n_valid",
                "mean_cf",
                "lat_min",
                "lat_max",
                "lon_min_0360",
                "lon_max_0360",
            ]
        )
        for r in results:
            writer.writerow(
                [
                    r.region,
                    str(r.file),
                    r.n_time,
                    r.n_valid,
                    f"{r.mean_cf:.8f}",
                    f"{r.lat_min:.4f}",
                    f"{r.lat_max:.4f}",
                    f"{r.lon_min_0360:.4f}",
                    f"{r.lon_max_0360:.4f}",
                ]
            )
    logger.info(f"✓ 已保存区域日志：{log_csv}")


def positive_int(value: str) -> int:
    out = int(value)
    if out <= 0:
        raise argparse.ArgumentTypeError("必须是正整数")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="用 Cartopy 绘制风光容量因子年度平均图。")
    parser.add_argument(
        "--cf_root",
        required=True,
        type=Path,
        help="容量因子模型目录，例如 output/CFs_of_solar/MPI-ESM1-2-HR。该目录下通常是多个国家/地区子目录。",
    )
    parser.add_argument("--technology", default="auto", choices=["auto", "solar", "wind"], help="容量因子类型，默认 auto。")
    parser.add_argument("--scenario", default="ssp126", help="SSP 情景，默认 ssp126。")
    parser.add_argument("--year", required=True, type=int, help="要绘制哪一年的年度平均容量因子，例如 2030。")
    parser.add_argument(
        "--countries",
        nargs="*",
        default=None,
        help="要绘制的国家/地区列表；默认 cf_root 下所有国家/地区。支持空格或逗号分隔。",
    )
    parser.add_argument(
        "--overlap_method",
        default="mean",
        choices=["mean", "first", "last", "max"],
        help="多个地区矩形重叠时的处理方式。默认 mean。",
    )
    parser.add_argument("--output_dir", default=Path("figs/CF_maps"), type=Path, help="输出图片目录，默认 figs/CF_maps。图片文件名会自动拼接。")
    parser.add_argument("--output_png", default=None, type=Path, help="可选：手动指定完整输出图片路径。通常不需要使用。")
    parser.add_argument("--save_nc", default=None, type=Path, help="可选：保存合并后的全球年度平均 nc。")
    parser.add_argument("--log_csv", default=None, type=Path, help="可选：保存每个 region 使用的文件和统计信息。")
    parser.add_argument("--title", default=None, help="可选：自定义图标题。")
    parser.add_argument("--extent", default="world", choices=["data", "world"], help="绘图范围：world 为全球；data 为仅显示有数据区域。默认 world。")
    parser.add_argument("--pad_degree", default=2.0, type=float, help="extent=data 时，在有效数据边界外扩多少度。默认 2。")
    parser.add_argument("--cmap", default="YlOrRd", help="matplotlib colormap，默认 YlOrRd，低值更接近浅色，避免大面积 0 显示成紫色。")
    parser.add_argument("--zero_threshold", default=1e-6, type=float, help="仅绘图时把 <= 此阈值的 CF 作为无效值显示，默认 1e-6。设为 -1 可不遮罩 0 值。")
    parser.add_argument("--zero_color", default="white", help="被遮罩的 0/近 0 值显示颜色，默认 white。")
    parser.add_argument("--vmin", default=0.0, type=float, help="色标最小值，默认 0。")
    parser.add_argument("--vmax", default=1.0, type=float, help="色标最大值，默认 1。")
    parser.add_argument("--dpi", default=300, type=positive_int, help="图片 DPI，默认 300。")
    parser.add_argument("--fig_width", default=12.0, type=float, help="图片宽度，默认 12 inch。")
    parser.add_argument("--fig_height", default=6.5, type=float, help="图片高度，默认 6.5 inch。")
    parser.add_argument("--no_features", action="store_true", help="不添加海岸线/国界等 Cartopy 地图要素。")
    args = parser.parse_args()

    cf_root = args.cf_root.expanduser().resolve()
    if not cf_root.exists():
        raise FileNotFoundError(f"cf_root 不存在：{cf_root}")

    countries = split_country_args(args.countries)
    regions = discover_regions(cf_root, countries)
    technology = normalize_technology(args.technology, cf_root)
    model = infer_model_name(cf_root)

    logger.info(f"cf_root       = {cf_root}")
    logger.info(f"model         = {model}")
    logger.info(f"technology    = {technology}")
    logger.info(f"scenario      = {args.scenario}")
    logger.info(f"year          = {args.year}")
    logger.info(f"regions       = {len(regions)}")
    logger.info(f"overlap       = {args.overlap_method}")
    logger.info(f"extent        = {args.extent}")
    logger.info(f"cmap          = {args.cmap}")
    logger.info(f"zero_threshold= {args.zero_threshold}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        global_cf, global_count, results, _ = merge_regions_to_global(
            cf_root=cf_root,
            regions=regions,
            scenario=args.scenario,
            year=args.year,
            technology=technology,
            overlap_method=args.overlap_method,
        )

    output_png = args.output_png
    if output_png is None:
        output_png = build_output_png_path(
            output_dir=args.output_dir.expanduser(),
            technology=technology,
            model=model,
            scenario=args.scenario,
            year=args.year,
            regions=regions,
            countries_was_specified=countries is not None,
            extent=args.extent,
            overlap_method=args.overlap_method,
        )
    else:
        output_png = output_png.expanduser()

    logger.info(f"output_png   = {output_png}")

    plot_cf_map(
        global_cf=global_cf,
        technology=technology,
        model=model,
        scenario=args.scenario,
        year=args.year,
        output_png=output_png,
        title=args.title,
        extent_mode=args.extent,
        pad_degree=args.pad_degree,
        cmap=args.cmap,
        zero_threshold=args.zero_threshold,
        zero_color=args.zero_color,
        vmin=args.vmin,
        vmax=args.vmax,
        dpi=args.dpi,
        figsize=(args.fig_width, args.fig_height),
        add_features=not args.no_features,
    )

    if args.save_nc is not None:
        save_merged_nc(
            save_nc=args.save_nc,
            global_cf=global_cf,
            global_count=global_count,
            technology=technology,
            model=model,
            scenario=args.scenario,
            year=args.year,
            overlap_method=args.overlap_method,
            regions=regions,
        )

    write_region_log(results, args.log_csv)


if __name__ == "__main__":
    main()

"""


示例：
python plot_cf.py \
    --cf_root output/CFs_of_solar/MPI-ESM1-2-HR \
    --scenario ssp126 \
    --year 2060 

python plot_cf.py \
    --cf_root output/CFs_of_wind/MPI-ESM1-2-HR \
    --scenario ssp126 \
    --year 2060 

"""