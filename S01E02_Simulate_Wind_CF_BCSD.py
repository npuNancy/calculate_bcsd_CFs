"""
区域逐时/逐3小时风电容量因子计算（BCSD uas + vas）
=================================================

本脚本用于基于 BCSD 降尺度气候数据计算区域陆上风电容量因子。

计算逻辑：
  1. 从 BCSD 数据中读取 uas 和 vas 两个近地面风速分量；
  2. 将 uas 和 vas 合成为 10 m 风速：
       V10 = sqrt(uas^2 + vas^2)
  3. 使用幂律风速外推方法，将 10 m 风速外推到轮毂高度：
       Vhub = V10 * (Hhub / 10)^alpha
     其中：
       Hhub = 100 m
       alpha = 1/7
  4. 使用 windpowerlib 读取 GE120/2500 风机功率曲线；
  5. 将轮毂高度风速代入功率曲线，计算风机输出功率；
  6. 容量因子定义为：
       wind_cf = power / rated_power
  7. 切出风速设置为 25 m s-1；
  8. 不额外乘风场效率系数；
  9. 只计算陆上风电容量因子；
 10. 按时间分块计算，并逐块写入 NetCDF 文件，避免一次性占用过多内存。

风机参数：
  - TURBINE_TYPE = "GE120/2500"
  - HUB_HEIGHT = 100 m
  - REF_HEIGHT = 10 m
  - POWER_LAW_ALPHA = 1/7
  - CUT_OUT = 25 m s-1

输入数据路径格式：
  {data_dir}/{model}/{region}/{model}/*.nc

例如：
  data/bcsd_outputs/MIROC-ES2H/Austria/MIROC-ES2H/
    uas_3h_bcsd_on_0p1deg_Austria_MIROC-ES2H_ssp126_2015-2060.nc
    vas_3h_bcsd_on_0p1deg_Austria_MIROC-ES2H_ssp126_2015-2060.nc

输入变量要求：
  - uas 文件中包含 uas_bcsd 或唯一数据变量；
  - vas 文件中包含 vas_bcsd 或唯一数据变量；
  - 数据维度通常为 (time, lat, lon)。

输出：
  {output_dir}/{model}/{region}/wind_CF_{region}_{model}_{scenario}_{years}_{months}.nc

输出变量：
  - wind_cf(time, lat, lon)
  - 单位：1
  - 数值范围：[0, 1]

注意：
  - 本脚本不读取 roughness_z0_01deg.nc。
  - 本脚本不使用空间粗糙度，也不使用对数律外推。
  - 本脚本不把 3h 数据插值到小时；输出时间分辨率与输入一致。
  - ERA5-Land 或 BCSD 陆地区域数据通常不包含海上网格，因此本脚本只处理陆上风电。
  - region 可指定单个区域，也可设置为 all 批量处理 model 下所有有效区域。
"""

from __future__ import annotations

import argparse
import gc
import glob
import logging
from pathlib import Path
from typing import Iterable

import numpy as np
import xarray as xr
from netCDF4 import Dataset
from tqdm import tqdm
from windpowerlib import WindTurbine
from global_land_mask import globe

# ── 日志配置 ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── 风机与风速外推参数 ───────────────────────────────────────────────────────
TURBINE_TYPE = "GE120/2500"  # GE 2.5-120；windpowerlib 数据库中与论文 GE 2.5 MW 相近的型号
HUB_HEIGHT = 100.0  # m，论文 onshore wind 使用 100 m hub height
REF_HEIGHT = 10.0  # m，uas/vas 为 10 m 风速分量
POWER_LAW_ALPHA = 1.0 / 7.0  # 论文采用的典型幂律指数
CUT_OUT = 25.0  # m s-1，论文风速转换为电力的上限


def parse_months(months: str) -> list[int]:
    """解析月份参数。空字符串表示 1-12 月。"""
    if not months.strip():
        return list(range(1, 13))
    out = [int(m) for m in months.split(",") if m.strip()]
    bad = [m for m in out if m < 1 or m > 12]
    if bad:
        raise ValueError(f"月份必须在 1-12 之间，收到：{bad}")
    return out


def parse_years(years: str) -> tuple[int, int]:
    """解析年份参数，格式 YYYY-YYYY 或 YYYY。"""
    years = years.strip()
    if "-" in years:
        y0, y1 = map(int, years.split("-", 1))
    else:
        y0 = y1 = int(years)
    if y1 < y0:
        raise ValueError(f"years 不合法：{years}")
    return y0, y1


def months_tag(months: str) -> str:
    """输出文件名中的月份标签。"""
    ms = parse_months(months)
    if len(ms) == 12:
        return "allmonths"
    return "m" + "-".join(f"{m:02d}" for m in ms)


def is_valid_region_name(region: str) -> bool:
    """过滤无效 region 目录。"""
    return bool(region) and not region.startswith("_") and not region.endswith("_") and not region.endswith("_repeated")


def discover_regions(data_dir: str | Path, model: str) -> list[str]:
    """发现某个 model 下的所有有效 region。"""
    root = Path(data_dir) / model
    regions = []
    for p in sorted(root.glob("*")):
        if p.is_dir() and is_valid_region_name(p.name):
            regions.append(p.name)
    return regions


def find_bcsd_file(data_dir: str | Path, model: str, region: str, scenario: str, var: str) -> Path:
    """使用 glob 查找指定变量文件。"""
    base = Path(data_dir) / model / region / model
    patterns = [
        str(base / f"{var}_3h_bcsd_on_0p1deg_{region}_{model}_{scenario}_*.nc"),
        str(base / f"{var}_*_{region}_{model}_{scenario}_*.nc"),
        str(base / f"{var}_*{scenario}*.nc"),
    ]
    files: list[str] = []
    for pat in patterns:
        files = sorted(glob.glob(pat))
        if files:
            break
    if not files:
        raise FileNotFoundError(f"找不到 {var} 文件。已尝试：\n  " + "\n  ".join(patterns))
    if len(files) > 1:
        logger.warning(f"{var} 匹配到多个文件，将使用第一个：{files[0]}")
    return Path(files[0])


def get_var_name(ds: xr.Dataset, preferred: str) -> str:
    """获取变量名。优先使用 {preferred}_bcsd，否则在唯一 data_var 时使用该变量。"""
    candidate = f"{preferred}_bcsd"
    if candidate in ds.data_vars:
        return candidate
    if preferred in ds.data_vars:
        return preferred
    if len(ds.data_vars) == 1:
        return list(ds.data_vars)[0]
    raise KeyError(f"无法在数据集中确定变量 {preferred}，可用变量：{list(ds.data_vars)}")


def find_coord_name(ds: xr.Dataset, candidates: Iterable[str]) -> str:
    """兼容 time/lat/lon 或 valid_time/latitude/longitude 等命名。"""
    for c in candidates:
        if c in ds.coords or c in ds.dims:
            return c
    raise KeyError(f"找不到坐标名，候选：{list(candidates)}；数据维度：{list(ds.dims)}")


def rename_coords_to_match(
    ds: xr.Dataset,
    target_time: str,
    target_lat: str,
    target_lon: str,
) -> xr.Dataset:
    """把另一个数据集的 time/lat/lon 坐标名改成与主数据集一致。"""
    src_time = find_coord_name(ds, ["time", "valid_time"])
    src_lat = find_coord_name(ds, ["lat", "latitude"])
    src_lon = find_coord_name(ds, ["lon", "longitude"])

    rename_map = {}
    if src_time != target_time:
        rename_map[src_time] = target_time
    if src_lat != target_lat:
        rename_map[src_lat] = target_lat
    if src_lon != target_lon:
        rename_map[src_lon] = target_lon

    if rename_map:
        ds = ds.rename(rename_map)

    return ds


def build_land_mask(lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    """根据输出网格 lat/lon 生成陆地掩膜 (n_lat, n_lon)，True=陆地。

    global_land_mask 要求经度 ∈ [-180, 180]、纬度 ∈ [-90, 90]，
    因此对 [0, 360) 经度做归一化，并夹紧纬度避免浮点越界。
    """
    lats = np.asarray(lats, dtype=np.float64)
    lons = np.asarray(lons, dtype=np.float64)

    lons_conv = ((lons + 180.0) % 360.0) - 180.0  # [0,360) -> [-180,180)
    lats_clip = np.clip(lats, -90.0, 90.0)

    lon2d, lat2d = np.meshgrid(lons_conv, lats_clip)  # (n_lat, n_lon)
    return globe.is_land(lat2d, lon2d).astype(bool)


def wind_to_ms(arr: np.ndarray, units: str | None) -> np.ndarray:
    """把风速转换为 m s-1。BCSD uas/vas 通常已经是 m s-1。"""
    x = arr.astype(np.float32, copy=False)
    units_l = (units or "").lower().replace(" ", "")
    if "km/h" in units_l or "kmh-1" in units_l or "kmhr-1" in units_l:
        x = x / 3.6
    # 对 m s-1、m/s、空单位等情况不做额外转换。
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    return x.astype(np.float32)


def build_time_index(time_da: xr.DataArray, years: str, months: str) -> np.ndarray:
    """根据 years/months 生成时间索引。"""
    y0, y1 = parse_years(years)
    month_list = parse_months(months)

    mask = (time_da.dt.year >= y0) & (time_da.dt.year <= y1) & time_da.dt.month.isin(month_list)
    idx = np.where(mask.values)[0]
    if idx.size == 0:
        raise ValueError(f"没有匹配到时间步：years={years}, months={months or 'all'}")
    return idx


def power_law_ratio(
    hub_height: float = HUB_HEIGHT, ref_height: float = REF_HEIGHT, alpha: float = POWER_LAW_ALPHA
) -> np.float32:
    """论文幂律风速外推系数。"""
    return np.float32((hub_height / ref_height) ** alpha)


def _get_power_curve_arrays() -> tuple[np.ndarray, np.ndarray, float]:
    """从 windpowerlib 数据库加载风机功率曲线，返回 (风速, 功率 kW, 额定功率 kW)。"""
    turbine = WindTurbine(turbine_type=TURBINE_TYPE, hub_height=HUB_HEIGHT)
    pc = turbine.power_curve
    ws = pc["wind_speed"].values.astype(np.float32)
    pw_kw = pc["value"].values.astype(np.float32) / 1000.0  # W → kW
    rated_kw = float(turbine.nominal_power / 1000.0)  # W → kW

    # 按风速升序，并只保留 cut-out 以内的功率曲线。
    order = np.argsort(ws)
    ws = ws[order]
    pw_kw = pw_kw[order]
    keep = ws <= CUT_OUT
    ws = ws[keep]
    pw_kw = pw_kw[keep]

    if ws.size == 0:
        raise RuntimeError(f"风机 {TURBINE_TYPE} 的功率曲线在 CUT_OUT={CUT_OUT} 以内为空。")

    # 若数据库功率曲线没有覆盖到 25 m/s，则把最后一个功率值延伸到 cut-out。
    if ws[-1] < CUT_OUT:
        last_pw = pw_kw[-1]
        extra_ws = np.arange(ws[-1] + 0.5, CUT_OUT + 0.001, 0.5, dtype=np.float32)
        if extra_ws.size and extra_ws[-1] < CUT_OUT:
            extra_ws = np.append(extra_ws, np.float32(CUT_OUT))
        elif extra_ws.size == 0:
            extra_ws = np.array([CUT_OUT], dtype=np.float32)
        extra_pw = np.full_like(extra_ws, last_pw, dtype=np.float32)
        ws = np.concatenate([ws, extra_ws])
        pw_kw = np.concatenate([pw_kw, extra_pw])
    elif ws[-1] > CUT_OUT:
        # 理论上不会发生，因为已经 keep <= CUT_OUT；保留防御性处理。
        ws = np.append(ws[ws < CUT_OUT], np.float32(CUT_OUT))
        pw_kw = np.append(pw_kw[: len(ws) - 1], pw_kw[-1])

    return ws.astype(np.float32), pw_kw.astype(np.float32), rated_kw


def apply_power_curve(ws_hub: np.ndarray, ws_curve: np.ndarray, pw_curve: np.ndarray, rated_kw: float) -> np.ndarray:
    """
    向量化的功率曲线插值，返回容量因子 CF。

    不再乘 FARM_EFFICIENCY。容量因子仅由轮毂高度风速和单机功率曲线决定。
    """
    # 在 0 m/s 以下和 CUT_OUT 以上插值为 0。
    ws_ext = np.concatenate([[0.0], ws_curve, [CUT_OUT + 0.01]]).astype(np.float32)
    pw_ext = np.concatenate([[0.0], pw_curve, [0.0]]).astype(np.float32)

    power = np.interp(ws_hub, ws_ext, pw_ext).astype(np.float32)  # kW
    cf = power / np.float32(rated_kw)
    cf = np.clip(cf, 0.0, 1.0)
    return cf.astype(np.float32)


def compute_wind_cf_chunk(
    uas: np.ndarray,
    vas: np.ndarray,
    ws_curve: np.ndarray,
    pw_curve: np.ndarray,
    rated_kw: float,
    extrap_ratio: np.float32,
) -> np.ndarray:
    """计算一个时间块的陆上风电容量因子。"""
    ws10 = np.sqrt(uas.astype(np.float32) ** 2 + vas.astype(np.float32) ** 2)
    ws_hub = ws10 * extrap_ratio
    cf = apply_power_curve(ws_hub, ws_curve, pw_curve, rated_kw)
    return cf.astype(np.float32)


def create_output_file(
    out_file: Path,
    template_file: Path,
    time_idx: np.ndarray,
    lat_values: np.ndarray,
    lon_values: np.ndarray,
    time_name: str,
    lat_name: str,
    lon_name: str,
    chunk_time: int,
    attrs: dict,
    compress_level: int = 4,
) -> Dataset:
    """创建 NetCDF 输出文件，并写入坐标。"""
    out_file.parent.mkdir(parents=True, exist_ok=True)

    # 用 decode_times=False 复制原始 time 数值和属性，避免 cftime 写出问题。
    raw = xr.open_dataset(template_file, decode_times=False)
    raw_time = raw[time_name].values[time_idx]
    time_attrs = dict(raw[time_name].attrs)
    lat_attrs = dict(raw[lat_name].attrs) if lat_name in raw else {}
    lon_attrs = dict(raw[lon_name].attrs) if lon_name in raw else {}
    raw.close()

    nc = Dataset(out_file, "w", format="NETCDF4")
    nc.createDimension("time", len(time_idx))
    nc.createDimension("lat", len(lat_values))
    nc.createDimension("lon", len(lon_values))

    tvar = nc.createVariable("time", raw_time.dtype, ("time",))
    tvar[:] = raw_time
    for k, v in time_attrs.items():
        try:
            setattr(tvar, k, v)
        except Exception:
            pass

    latvar = nc.createVariable("lat", "f4", ("lat",))
    lonvar = nc.createVariable("lon", "f4", ("lon",))
    latvar[:] = lat_values.astype(np.float32)
    lonvar[:] = lon_values.astype(np.float32)
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

    chunksizes = (min(chunk_time, len(time_idx)), len(lat_values), len(lon_values))
    cf = nc.createVariable(
        "wind_cf",
        "f4",
        ("time", "lat", "lon"),
        zlib=True,
        complevel=compress_level,
        chunksizes=chunksizes,
        fill_value=np.float32(np.nan),
    )
    cf.long_name = "Onshore wind capacity factor"
    cf.units = "1"
    cf.description = (
        "Onshore wind capacity factor computed from uas/vas using power-law extrapolation "
        f"Vhub=V10*(Hhub/10)^alpha with Hhub={HUB_HEIGHT:g} m, alpha={POWER_LAW_ALPHA:.6f}, "
        f"turbine={TURBINE_TYPE}, cut_out={CUT_OUT:g} m s-1. No farm-efficiency multiplier is applied."
    )

    for k, v in attrs.items():
        try:
            setattr(nc, k, v)
        except Exception:
            pass

    return nc


def compute_region_wind_cf(
    data_dir: str,
    model: str,
    region: str,
    scenario: str,
    years: str,
    months: str,
    output_dir: str,
    chunk_time: int = 512,
    overwrite: bool = False,
    compress_level: int = 4,
) -> Path:
    """计算单个 model-region-scenario 的陆上风电容量因子。"""
    if not is_valid_region_name(region):
        raise ValueError(f"region 名称不合法或被过滤：{region}")

    uas_file = find_bcsd_file(data_dir, model, region, scenario, "uas")
    vas_file = find_bcsd_file(data_dir, model, region, scenario, "vas")

    logger.info(f"uas: {uas_file}")
    logger.info(f"vas: {vas_file}")

    out_file = (
        Path(output_dir) / model / region / f"wind_CF_{region}_{model}_{scenario}_{years}_{months_tag(months)}.nc"
    )
    if out_file.exists() and not overwrite:
        logger.info(f"✓ 已存在，跳过：{out_file}")
        return out_file

    ds_uas = xr.open_dataset(uas_file)
    ds_vas = xr.open_dataset(vas_file)

    uas_var = get_var_name(ds_uas, "uas")
    vas_var = get_var_name(ds_vas, "vas")

    time_name = find_coord_name(ds_uas, ["time", "valid_time"])
    lat_name = find_coord_name(ds_uas, ["lat", "latitude"])
    lon_name = find_coord_name(ds_uas, ["lon", "longitude"])

    # 让 vas 的坐标名与 uas 保持一致，避免一个叫 time、另一个叫 valid_time。
    ds_vas = rename_coords_to_match(
        ds_vas,
        target_time=time_name,
        target_lat=lat_name,
        target_lon=lon_name,
    )

    # 记录对齐前尺寸，便于日志检查。
    sizes_before = {
        "uas": dict(ds_uas[uas_var].sizes),
        "vas": dict(ds_vas[vas_var].sizes),
    }

    # 关键修复：uas/vas 先按公共 time/lat/lon 取交集。
    # 这样即使某个变量少了最后一个时间步，也不会出现 isel 越界。
    ds_uas, ds_vas = xr.align(
        ds_uas,
        ds_vas,
        join="inner",
        copy=False,
    )

    sizes_after = {
        "uas": dict(ds_uas[uas_var].sizes),
        "vas": dict(ds_vas[vas_var].sizes),
    }

    if sizes_before != sizes_after:
        logger.warning(
            f"uas/vas 时间或空间坐标不完全一致，已按公共坐标取交集：before={sizes_before}, after={sizes_after}"
        )

    # 对齐后再检查尺寸。
    if ds_uas[uas_var].sizes != ds_vas[vas_var].sizes:
        raise ValueError(f"uas 与 vas 对齐后尺寸仍不一致：uas={ds_uas[uas_var].sizes}, " f"vas={ds_vas[vas_var].sizes}")

    # 固定变量维度顺序，避免某些文件维度顺序不同。
    uas_da = ds_uas[uas_var].transpose(time_name, lat_name, lon_name)
    vas_da = ds_vas[vas_var].transpose(time_name, lat_name, lon_name)

    # 对齐之后，再生成 time_idx。
    time_idx = build_time_index(ds_uas[time_name], years, months)

    lats = ds_uas[lat_name].values.astype(np.float32)
    lons = ds_uas[lon_name].values.astype(np.float32)

    n_time = len(time_idx)
    nlat = len(lats)
    nlon = len(lons)
    logger.info(f"匹配时间步: {n_time}；空间维度: lat={nlat}, lon={nlon}")

    land_mask = build_land_mask(lats, lons)  # (n_lat, n_lon) bool
    logger.info(f"陆地格点占比：{100.0 * land_mask.mean():.1f}%")

    ws_curve, pw_curve, rated_kw = _get_power_curve_arrays()
    extrap_ratio = power_law_ratio()
    logger.info(
        f"风速外推: Vhub=V10*({HUB_HEIGHT:.1f}/{REF_HEIGHT:.1f})^{POWER_LAW_ALPHA:.6f}，ratio={float(extrap_ratio):.4f}"
    )
    logger.info(f"风机: {TURBINE_TYPE}；rated={rated_kw:.1f} kW；cut-out={CUT_OUT:.1f} m/s")

    attrs = {
        "model": model,
        "region": region,
        "scenario": scenario,
        "years": years,
        "months": months if months.strip() else "1-12",
        "source_files": ", ".join(str(p) for p in [uas_file, vas_file]),
        "technology": "onshore wind only",
        "turbine_type": TURBINE_TYPE,
        "hub_height_m": HUB_HEIGHT,
        "reference_height_m": REF_HEIGHT,
        "power_law_alpha": POWER_LAW_ALPHA,
        "power_law_ratio": float(extrap_ratio),
        "cut_out_speed_ms": CUT_OUT,
        "farm_efficiency_multiplier": "not applied",
        "notes": "ERA5-Land/BCSD land data are used; offshore wind is not calculated in this script.",
        "ocean_handling": "ocean grid cells set to NaN via global-land-mask",
    }

    nc = create_output_file(
        out_file=out_file,
        template_file=uas_file,
        time_idx=time_idx,
        lat_values=lats,
        lon_values=lons,
        time_name=time_name,
        lat_name=lat_name,
        lon_name=lon_name,
        chunk_time=chunk_time,
        attrs=attrs,
        compress_level=compress_level,
    )
    cf_var = nc.variables["wind_cf"]

    uas_units = ds_uas[uas_var].attrs.get("units", "")
    vas_units = ds_vas[vas_var].attrs.get("units", "")

    total_sum = 0.0
    total_count = 0
    positive_sum = 0.0
    positive_count = 0
    max_cf = 0.0

    logger.info(f"开始分块计算，chunk_time={chunk_time}")
    out_start = 0
    try:
        for start in tqdm(range(0, n_time, chunk_time), desc=f"{region}-{scenario}", unit="块"):
            end = min(start + chunk_time, n_time)
            idx_chunk = time_idx[start:end]

            uas_raw = ds_uas[uas_var].isel({time_name: idx_chunk}).values
            vas_raw = ds_vas[vas_var].isel({time_name: idx_chunk}).values

            uas_ms = wind_to_ms(uas_raw, uas_units)
            vas_ms = wind_to_ms(vas_raw, vas_units)

            cf_chunk = compute_wind_cf_chunk(
                uas=uas_ms,
                vas=vas_ms,
                ws_curve=ws_curve,
                pw_curve=pw_curve,
                rated_kw=rated_kw,
                extrap_ratio=extrap_ratio,
            )

            # 海洋格点恢复为 NaN（修复 nan_to_num 把海洋填 0 的问题）
            cf_chunk = np.where(land_mask[None, :, :], cf_chunk, np.nan).astype(np.float32)

            cf_var[out_start : out_start + (end - start), :, :] = cf_chunk
            out_start += end - start

            finite = np.isfinite(cf_chunk)
            if np.any(finite):
                vals = cf_chunk[finite]
                total_sum += float(vals.sum(dtype=np.float64))
                total_count += int(vals.size)
                if vals.size:
                    max_cf = max(max_cf, float(vals.max()))
                pos = vals[vals > 0]
                if pos.size:
                    positive_sum += float(pos.sum(dtype=np.float64))
                    positive_count += int(pos.size)

            del uas_raw, vas_raw, uas_ms, vas_ms, cf_chunk
            gc.collect()
    except BaseException:
        if nc.isopen():
            nc.close()
        if out_file.exists():
            logger.error(f"计算出错，删除不完整的文件：{out_file}")
            out_file.unlink()
        raise
    finally:
        if nc.isopen():
            nc.close()
        ds_uas.close()
        ds_vas.close()

    logger.info(f"✓ 已保存：{out_file}")
    if total_count:
        logger.info(f"  整体平均 CF   : {total_sum / total_count:.4f}")
        logger.info(f"  全局最大 CF   : {max_cf:.4f}")
        logger.info(f"  有风时段占比  : {100.0 * positive_count / total_count:.2f}%")
    if positive_count:
        logger.info(f"  有风时段平均 CF: {positive_sum / positive_count:.4f}")

    return out_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description="BCSD uas/vas 陆上风电容量因子计算（幂律外推 + windpowerlib 功率曲线）"
    )
    parser.add_argument("--data_dir", default="data/bcsd_outputs", help="数据目录，例如 data/bcsd_outputs")
    parser.add_argument("--model", required=True, help="模型名称，例如 MIROC-ES2H")
    parser.add_argument("--region", required=True, help="区域名称，例如 Austria；也可传 all 处理全部有效 region")
    parser.add_argument("--scenario", required=True, help="情景名称，例如 ssp126")
    parser.add_argument("--years", default="2015-2060", help="年段 YYYY-YYYY，包含两端")
    parser.add_argument("--months", default="", help="逗号分隔月份；默认全年 1-12")
    parser.add_argument("--output_dir", default="output/CFs_of_wind", help="输出目录")
    parser.add_argument("--chunk_time", type=int, default=512, help="每块处理的时间步数；3h 数据默认 512")
    parser.add_argument("--compress_level", type=int, default=4, help="NetCDF 压缩级别 0-9")
    parser.add_argument("--overwrite", action="store_true", help="若输出文件已存在，则覆盖重算")
    args = parser.parse_args()

    if args.region.lower() == "all":
        regions = discover_regions(args.data_dir, args.model)
        if not regions:
            raise RuntimeError(f"没有发现有效 region：{Path(args.data_dir) / args.model}")
    else:
        regions = [args.region]

    logger.info(f"待处理 region 数量：{len(regions)}")
    for i, region in enumerate(regions, 1):
        logger.info(f"\n{'=' * 80}")
        logger.info(f"[{i}/{len(regions)}] model={args.model}, region={region}, scenario={args.scenario}")
        logger.info(f"{'=' * 80}")
        compute_region_wind_cf(
            data_dir=args.data_dir,
            model=args.model,
            region=region,
            scenario=args.scenario,
            years=args.years,
            months=args.months,
            output_dir=args.output_dir,
            chunk_time=args.chunk_time,
            overwrite=args.overwrite,
            compress_level=args.compress_level,
        )


if __name__ == "__main__":
    main()

"""

python S01E02_Simulate_Wind_CF_BCSD.py \
  --data_dir data/bcsd_outputs \
  --model MIROC-ES2H \
  --region all \
  --scenario ssp126 \
  --years 2015-2060


"""
