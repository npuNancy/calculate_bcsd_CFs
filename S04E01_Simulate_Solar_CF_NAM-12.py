"""
CMIP6-CORDEX NAM-12 逐小时光伏容量因子计算（rsds + tas + uas/vas）
=================================================================

本脚本用于基于 CMIP6-CORDEX NAM-12 逐小时数据计算光伏容量因子。

计算逻辑：
  1. 从 CORDEX 数据中读取 rsds、tas、uas、vas 四个变量；
  2. rsds 作为 GHI，单位自动转换为 kW m-2；
  3. tas 自动转换为摄氏度；
  4. 以 rsds 的半点时间轴为主，将瞬时变量 tas/uas/vas 线性插值到 rsds 时间；
  5. uas 和 vas 合成为近地面风速，用于光伏组件温度修正；
  6. 根据 UTC 时间、纬度和经度计算太阳天顶角（使用二维 lat/lon）；
  7. 使用 Erbs 模型将 GHI 分解为直接辐射和散射辐射；
  8. 假设光伏组件采用双轴跟踪；
  9. 使用论文中的温度修正方法；
 10. 使用 SYS_coef = 0.8056；
 11. 计算光伏容量因子，限制在 [0, 1]；
 12. 按年循环处理，每年内按时间分块计算，逐块写入 NetCDF。

时间对齐：
  - rsds 为时间平均（cell_methods: area: time: mean），时间标签在 :30；
  - tas/uas/vas 为瞬时值（cell_methods: area: mean time: point），时间标签在 :00；
  - 两者相差 30 分钟，采用方案4：保留 rsds 时间轴，线性插值 tas/uas/vas。

输入数据路径格式：
  {data_dir}/{gcm_model}/{realization}/{rcm_model}/{scenario}/{var}/
  例如：data/CORDEX-CMIP6/NAM-12/1hr/MPI-ESM1-2-LR/r1i1p1f1/CRCM5/ssp126/rsds/

输入变量要求：
  - rsds/tas/uas/vas 文件中包含同名变量；
  - 数据维度为 (time, rlat, rlon)，旋转极网格；
  - 每年一个文件。

输出：
  {output_dir}/{gcm_model}/{realization}/
    solar_CF_NAM-12_{gcm_model}_{realization}_{rcm_model}_{scenario}_{years}_{months_tag}.nc

输出变量：
  - solar_cf(time, rlat, rlon)
  - 单位：1
  - 数值范围：[0, 1]

示例：
  python S04E01_Simulate_Solar_CF_NAM-12.py \\
    --gcm_model MPI-ESM1-2-LR --realization r1i1p1f1 \\
    --scenario ssp126 --years 2020-2060
"""

from __future__ import annotations

import argparse
import gc
import glob
import logging
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import xarray as xr
from netCDF4 import Dataset
from tqdm import tqdm

# ── 日志配置 ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── Erbs 辐射分解参数 ────────────────────────────────────────────────────────
SOLAR_CONSTANT_KW = 1.367  # kW m-2
MIN_COS_ZENITH = 1e-4

# ── PV 参数 ──────────────────────────────────────────────────────────────────
SYS_COEF = 0.8056
T_STC = 25.0  # °C
GAMMA_TEMP = 0.005  # °C-1
C1 = 4.3  # °C
C2 = 0.943
C3 = 0.028  # °C m2 W-1
C4 = 1.528  # °C s m-1


# ─────────────────────────────────────────────────────────────────────────────
# 0. 通用解析与数据定位
# ─────────────────────────────────────────────────────────────────────────────


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


def find_coord_name(ds: xr.Dataset, candidates: Iterable[str]) -> str:
    """兼容不同坐标命名。"""
    for c in candidates:
        if c in ds.coords or c in ds.dims:
            return c
    raise KeyError(f"找不到坐标名，候选：{list(candidates)}；数据维度：{list(ds.dims)}")


def get_var_name(ds: xr.Dataset, preferred: str) -> str:
    """获取变量名。"""
    if preferred in ds.data_vars:
        return preferred
    if len(ds.data_vars) == 1:
        return list(ds.data_vars)[0]
    raise KeyError(f"无法在数据集中确定变量 {preferred}，可用变量：{list(ds.data_vars)}")


def prepare_dataarray(
    ds: xr.Dataset,
    preferred_var: str,
    time_name: str,
    rlat_name: str,
    rlon_name: str,
) -> xr.DataArray:
    """获取变量，压缩标量维度，转置为 (time, rlat, rlon)。"""
    var = get_var_name(ds, preferred_var)
    da = ds[var]
    for dim in list(da.dims):
        if dim not in {time_name, rlat_name, rlon_name}:
            if da.sizes[dim] == 1:
                da = da.isel({dim: 0}, drop=True)
            else:
                raise ValueError(f"变量 {var} 存在非单例额外维度 {dim}={da.sizes[dim]}")
    return da.transpose(time_name, rlat_name, rlon_name)


# ─────────────────────────────────────────────────────────────────────────────
# 1. CORDEX 文件发现
# ─────────────────────────────────────────────────────────────────────────────


def cordex_var_files(
    data_dir: str | Path,
    gcm_model: str,
    realization: str,
    rcm_model: str,
    scenario: str,
    var: str,
    years_tuple: tuple[int, int],
) -> dict[int, Path]:
    """发现 CORDEX 某变量的所有年度文件，返回 {年份: 路径}。"""
    var_dir = Path(data_dir) / gcm_model / realization / rcm_model / scenario / var
    if not var_dir.is_dir():
        raise FileNotFoundError(f"变量目录不存在：{var_dir}")
    pattern = str(var_dir / f"{var}_NAM-12_{gcm_model}_*.nc")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"找不到 {var} 文件。已尝试：{pattern}")

    y0, y1 = years_tuple
    result: dict[int, Path] = {}
    for f in files:
        m = re.search(r"_(\d{4})\d{8}-\d{12}\.nc$", Path(f).name)
        if m:
            year = int(m.group(1))
            if y0 <= year <= y1:
                result[year] = Path(f)
    if not result:
        raise FileNotFoundError(f"在 {var_dir} 中找不到 {y0}-{y1} 范围内的 {var} 文件")
    return result


# ─────────────────────────────────────────────────────────────────────────────
# 2. 单位转换
# ─────────────────────────────────────────────────────────────────────────────


def to_kw_m2(arr: np.ndarray, units: str | None) -> np.ndarray:
    """把辐照度转换为 kW m-2。"""
    x = arr.astype(np.float32, copy=False)
    units_l = (units or "").lower().replace(" ", "")
    if "kw" in units_l:
        out = x
    elif "w" in units_l and "m" in units_l:
        out = x / 1000.0
    else:
        finite = x[np.isfinite(x)]
        if finite.size and np.nanpercentile(finite, 99) > 10:
            out = x / 1000.0
        else:
            out = x
    out = np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)
    return np.maximum(out, 0.0).astype(np.float32)


def tas_to_celsius(arr: np.ndarray, units: str | None) -> np.ndarray:
    """把气温转换为摄氏度。"""
    x = arr.astype(np.float32, copy=False)
    units_l = (units or "").lower()
    if "k" in units_l and "deg" not in units_l:
        out = x - 273.15
    elif "c" in units_l or "celsius" in units_l:
        out = x
    else:
        finite = x[np.isfinite(x)]
        if finite.size and np.nanmedian(finite) > 100:
            out = x - 273.15
        else:
            out = x
    return out.astype(np.float32)


def wind_to_ms(arr: np.ndarray, units: str | None) -> np.ndarray:
    """把风速分量转换为 m s-1。"""
    x = arr.astype(np.float32, copy=False)
    units_l = (units or "").lower().replace(" ", "")
    if "km/h" in units_l or "kmh-1" in units_l:
        x = x / 3.6
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    return x.astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# 3. 辐射分解与太阳天顶角（二维 lat/lon）
# ─────────────────────────────────────────────────────────────────────────────


def _erbs_diffuse_fraction(kt: np.ndarray) -> np.ndarray:
    """Erbs (1982) 模型：由晴空指数 kt 估计漫射比例。"""
    fd = np.where(
        kt <= 0.22,
        1.0 - 0.09 * kt,
        np.where(
            kt <= 0.80,
            0.9511 - 0.1604 * kt + 4.388 * kt**2 - 16.638 * kt**3 + 12.336 * kt**4,
            0.165,
        ),
    )
    fd = np.clip(fd, 0.0, 1.0)
    fd = np.where(kt <= 0.0, 1.0, fd)
    return fd.astype(np.float32)


def solar_cos_zenith_2d(
    doy: np.ndarray,
    hour_decimal_utc: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
) -> np.ndarray:
    """
    计算太阳天顶角余弦 cos(SZA)，适配二维 lat/lon。

    lats/lons shape = (nrlat, nrlon)；
    doy/hour_decimal_utc shape = (ntime,)；
    输出 shape = (ntime, nrlat, nrlon)。
    """
    doy = doy.astype(np.float64)
    hour_decimal_utc = hour_decimal_utc.astype(np.float64)
    lats = lats.astype(np.float64)
    lons = lons.astype(np.float64)

    b = 2.0 * np.pi * (doy - 1.0) / 365.0

    decl = (
        0.006918
        - 0.399912 * np.cos(b)
        + 0.070257 * np.sin(b)
        - 0.006758 * np.cos(2 * b)
        + 0.000907 * np.sin(2 * b)
        - 0.002697 * np.cos(3 * b)
        + 0.001480 * np.sin(3 * b)
    )

    eot = 229.18 * (
        0.000075 + 0.001868 * np.cos(b) - 0.032077 * np.sin(b) - 0.014615 * np.cos(2 * b) - 0.040849 * np.sin(2 * b)
    )

    # 二维广播：time 维在最前，lat/lon 是二维数组
    true_solar_minutes = hour_decimal_utc[:, None, None] * 60.0 + eot[:, None, None] + 4.0 * lons[None, :, :]
    hour_angle = np.deg2rad(true_solar_minutes / 4.0 - 180.0)

    lat_rad = np.deg2rad(lats)
    cos_sza = np.sin(lat_rad)[None, :, :] * np.sin(decl)[:, None, None] + np.cos(lat_rad)[None, :, :] * np.cos(decl)[
        :, None, None
    ] * np.cos(hour_angle)
    return np.maximum(cos_sza, 0.0).astype(np.float32)


def extraterrestrial_horizontal_irradiance(
    doy: np.ndarray,
    cos_sza: np.ndarray,
) -> np.ndarray:
    """大气层外水平面辐照度，单位 kW m-2。"""
    b = 2.0 * np.pi * (doy.astype(np.float64) - 1.0) / 365.0
    e0 = 1.000110 + 0.034221 * np.cos(b) + 0.001280 * np.sin(b) + 0.000719 * np.cos(2 * b) + 0.000077 * np.sin(2 * b)
    et = SOLAR_CONSTANT_KW * e0
    return (et[:, None, None] * cos_sza).astype(np.float32)


def compute_diffuse_fraction_erbs(
    ghi_kw: np.ndarray,
    doy: np.ndarray,
    hour_decimal_utc: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Erbs 辐射分解。返回 diffuse fraction 和 cos(SZA)。"""
    cos_sza = solar_cos_zenith_2d(doy, hour_decimal_utc, lats, lons)
    etr = extraterrestrial_horizontal_irradiance(doy, cos_sza)

    with np.errstate(invalid="ignore", divide="ignore"):
        kt = np.where(etr > 1e-6, ghi_kw / etr, 0.0)
    kt = np.clip(kt, 0.0, 1.0)
    fd = _erbs_diffuse_fraction(kt)
    return fd, cos_sza


def plane_irradiance_two_axis_erbs(
    ghi_kw: np.ndarray,
    fd: np.ndarray,
    cos_sza: np.ndarray,
) -> np.ndarray:
    """双轴跟踪下的组件平面辐照度，单位 kW m-2。"""
    dhi = ghi_kw * fd
    beam_horiz = np.maximum(ghi_kw - dhi, 0.0)

    with np.errstate(invalid="ignore", divide="ignore"):
        dni = np.where(cos_sza > MIN_COS_ZENITH, beam_horiz / cos_sza, 0.0)

    plane_diffuse = dhi * (1.0 + cos_sza) / 2.0
    plane = dni + plane_diffuse
    plane = np.where(cos_sza > 0.0, plane, 0.0)
    plane = np.nan_to_num(plane, nan=0.0, posinf=0.0, neginf=0.0)
    return np.maximum(plane, 0.0).astype(np.float32)


def paper_temperature_coefficient(
    plane_irradiance_kw: np.ndarray,
    tas_c: np.ndarray,
    wind_speed: np.ndarray,
) -> np.ndarray:
    """论文温度修正。"""
    irradiance_w = plane_irradiance_kw * 1000.0
    t_cell = C1 + C2 * tas_c + C3 * irradiance_w - C4 * wind_speed
    tem_coef = 1.0 - GAMMA_TEMP * (t_cell - T_STC)
    tem_coef = np.nan_to_num(tem_coef, nan=0.0, posinf=0.0, neginf=0.0)
    return np.maximum(tem_coef, 0.0).astype(np.float32)


def compute_solar_cf_chunk(
    rsds_kw: np.ndarray,
    tas_c: np.ndarray,
    uas: np.ndarray,
    vas: np.ndarray,
    doy: np.ndarray,
    hour_decimal_utc: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
) -> np.ndarray:
    """计算一个时间块的光伏容量因子。"""
    fd, cos_sza = compute_diffuse_fraction_erbs(rsds_kw, doy, hour_decimal_utc, lats, lons)
    plane_kw = plane_irradiance_two_axis_erbs(rsds_kw, fd, cos_sza)

    wind_speed = np.sqrt(uas.astype(np.float32) ** 2 + vas.astype(np.float32) ** 2)
    temp_coef = paper_temperature_coefficient(plane_kw, tas_c, wind_speed)

    cf = plane_kw * temp_coef * SYS_COEF
    cf = np.clip(cf, 0.0, 1.0)
    return cf.astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# 4. 时间索引与时间插值
# ─────────────────────────────────────────────────────────────────────────────


def build_time_index(
    time_da: xr.DataArray,
    years: str,
    months: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """根据 years/months 生成时间索引，返回 (idx, doy, hour_decimal)。"""
    y0, y1 = parse_years(years)
    month_list = parse_months(months)

    mask = (time_da.dt.year >= y0) & (time_da.dt.year <= y1) & time_da.dt.month.isin(month_list)
    idx = np.where(mask.values)[0]
    if idx.size == 0:
        raise ValueError(f"没有匹配到时间步：years={years}, months={months or 'all'}")

    selected = time_da.isel({time_da.dims[0]: idx})
    doy = selected.dt.dayofyear.values.astype(np.float32)
    hour = selected.dt.hour.values.astype(np.float32)
    minute = selected.dt.minute.values.astype(np.float32) if hasattr(selected.dt, "minute") else 0.0
    second = selected.dt.second.values.astype(np.float32) if hasattr(selected.dt, "second") else 0.0
    hour_decimal = hour + minute / 60.0 + second / 3600.0
    return idx, doy, hour_decimal


def _datetime64_to_ns(values: np.ndarray) -> np.ndarray:
    """将时间坐标转换为数值（ns 整数或 cftime 天数），用于 searchsorted。"""
    import cftime

    arr = np.asarray(values)
    if np.issubdtype(arr.dtype, np.datetime64):
        return arr.astype("datetime64[ns]").astype(np.int64)
    if arr.dtype == object and arr.size > 0 and isinstance(arr.flat[0], cftime.datetime):
        # cftime 类型：转为距离第一个时间点的微秒偏移量
        base = arr.flat[0]
        return np.array([(v - base).total_seconds() for v in arr.flat], dtype=np.float64)
    return arr.astype(np.float64)


def interp_instantaneous_to_target_times_chunk(
    da: xr.DataArray,
    time_name: str,
    target_times: np.ndarray,
    *,
    fill_boundary: str = "nearest",
) -> np.ndarray:
    """
    将瞬时变量按时间线性插值到目标时间轴，返回 numpy 数组。

    方案4核心：target_times 来自 rsds.time（:30），
    da 来自 tas/uas/vas（:00），线性插值到 rsds 时间。
    """
    src_times = da[time_name].values
    src_num = _datetime64_to_ns(src_times)
    tgt_num = _datetime64_to_ns(target_times)

    if src_num.ndim != 1:
        raise ValueError(f"{da.name or 'variable'} 的时间坐标不是一维。")
    if src_num.size < 1:
        raise ValueError(f"{da.name or 'variable'} 的时间坐标为空。")
    if np.any(np.diff(src_num) <= 0):
        raise ValueError(f"{da.name or 'variable'} 的时间坐标不是严格递增。")

    nsrc = src_num.size
    right = np.searchsorted(src_num, tgt_num, side="left")
    right = np.clip(right, 0, nsrc - 1)
    left = np.clip(right - 1, 0, nsrc - 1)

    exact = src_num[right] == tgt_num
    left[exact] = right[exact]

    before = tgt_num <= src_num[0]
    after = tgt_num >= src_num[-1]
    if fill_boundary == "nearest":
        left[before] = 0
        right[before] = 0
        left[after] = nsrc - 1
        right[after] = nsrc - 1
    elif fill_boundary == "nan":
        pass
    else:
        raise ValueError(f"不支持 fill_boundary={fill_boundary!r}")

    i0 = int(min(left.min(), right.min()))
    i1 = int(max(left.max(), right.max())) + 1
    src_block = da.isel({time_name: slice(i0, i1)}).values.astype(np.float32)

    left_rel = left - i0
    right_rel = right - i0

    left_t = src_num[left].astype(np.float64)
    right_t = src_num[right].astype(np.float64)
    tgt_t = tgt_num.astype(np.float64)
    denom = right_t - left_t
    with np.errstate(invalid="ignore", divide="ignore"):
        w = np.where(denom != 0, (tgt_t - left_t) / denom, 0.0).astype(np.float32)
    w = w[:, None, None]

    out = src_block[left_rel] * (1.0 - w) + src_block[right_rel] * w
    if fill_boundary == "nan":
        out[before | after] = np.nan
    return out.astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# 5. 输出文件创建（旋转极网格）
# ─────────────────────────────────────────────────────────────────────────────


def create_output_file_rotated(
    out_file: Path,
    template_file: Path,
    total_time_steps: int,
    compress_level: int = 4,
    chunk_time: int = 24,
    attrs: dict | None = None,
) -> Dataset:
    """创建旋转极网格的 NetCDF 输出文件。"""
    out_file.parent.mkdir(parents=True, exist_ok=True)

    ds_tpl = xr.open_dataset(template_file)

    rlat_name = find_coord_name(ds_tpl, ["rlat"])
    rlon_name = find_coord_name(ds_tpl, ["rlon"])

    rlat_vals = ds_tpl[rlat_name].values
    rlon_vals = ds_tpl[rlon_name].values
    rlat_attrs = dict(ds_tpl[rlat_name].attrs)
    rlon_attrs = dict(ds_tpl[rlon_name].attrs)

    lat_2d = ds_tpl["lat"].values
    lon_2d = ds_tpl["lon"].values
    lat_attrs = dict(ds_tpl["lat"].attrs)
    lon_attrs = dict(ds_tpl["lon"].attrs)

    crs_attrs = dict(ds_tpl["crs"].attrs) if "crs" in ds_tpl.coords else {}

    raw = xr.open_dataset(template_file, decode_times=False)
    raw_time = raw["time"].values
    time_attrs = dict(raw["time"].attrs)
    raw.close()
    ds_tpl.close()

    nrlat = len(rlat_vals)
    nrlon = len(rlon_vals)

    nc = Dataset(out_file, "w", format="NETCDF4")

    if attrs:
        for k, v in attrs.items():
            try:
                setattr(nc, k, v)
            except Exception:
                pass

    nc.createDimension("time", total_time_steps)
    nc.createDimension(rlat_name, nrlat)
    nc.createDimension(rlon_name, nrlon)

    tvar = nc.createVariable("time", raw_time.dtype, ("time",))
    for k, v in time_attrs.items():
        try:
            setattr(tvar, k, v)
        except Exception:
            pass

    rlat_var = nc.createVariable(rlat_name, "f8", (rlat_name,))
    rlat_var[:] = rlat_vals
    for k, v in rlat_attrs.items():
        try:
            setattr(rlat_var, k, v)
        except Exception:
            pass

    rlon_var = nc.createVariable(rlon_name, "f8", (rlon_name,))
    rlon_var[:] = rlon_vals
    for k, v in rlon_attrs.items():
        try:
            setattr(rlon_var, k, v)
        except Exception:
            pass

    lat_var = nc.createVariable("lat", "f8", (rlat_name, rlon_name))
    lat_var[:] = lat_2d
    for k, v in lat_attrs.items():
        try:
            setattr(lat_var, k, v)
        except Exception:
            pass

    lon_var = nc.createVariable("lon", "f8", (rlat_name, rlon_name))
    lon_var[:] = lon_2d
    for k, v in lon_attrs.items():
        try:
            setattr(lon_var, k, v)
        except Exception:
            pass

    crs_var = nc.createVariable("crs", "i1")
    for k, v in crs_attrs.items():
        try:
            setattr(crs_var, k, v)
        except Exception:
            pass

    chunksizes = (min(chunk_time, total_time_steps), nrlat, nrlon)
    cf = nc.createVariable(
        "solar_cf",
        "f4",
        ("time", rlat_name, rlon_name),
        zlib=True,
        complevel=compress_level,
        chunksizes=chunksizes,
        fill_value=np.float32(np.nan),
    )
    cf.long_name = "Solar photovoltaic capacity factor"
    cf.units = "1"
    cf.grid_mapping = "crs"
    cf.description = (
        "PV capacity factor from CORDEX rsds using Erbs diffuse fraction, "
        "two-axis tracking, paper temperature correction, SYScoef=0.8056."
    )

    return nc


# ─────────────────────────────────────────────────────────────────────────────
# 6. 合并年度文件与单年计算
# ─────────────────────────────────────────────────────────────────────────────


def _merge_yearly_nc_files(
    yearly_files: list[Path],
    out_file: Path,
    cf_var_name: str,
    compress_level: int = 4,
    chunk_time: int = 24,
) -> Path:
    """合并逐年 NC 文件为一个文件。"""
    total_steps = 0
    for f in yearly_files:
        with Dataset(str(f), "r") as nc:
            total_steps += nc.dimensions["time"].size
    logger.info(f"合并 {len(yearly_files)} 个年度文件，总时间步：{total_steps}")

    nc_out = create_output_file_rotated(
        out_file=out_file,
        template_file=yearly_files[0],
        total_time_steps=total_steps,
        compress_level=compress_level,
        chunk_time=chunk_time,
    )
    cf_var_out = nc_out.variables[cf_var_name]
    tvar_out = nc_out.variables["time"]

    try:
        offset = 0
        for f in yearly_files:
            with Dataset(str(f), "r") as nc_in:
                n = nc_in.dimensions["time"].size
                tvar_out[offset : offset + n] = nc_in.variables["time"][:]
                cf_var_out[offset : offset + n, :, :] = nc_in.variables[cf_var_name][:, :, :]
                offset += n

        nc_out.close()
    except BaseException:
        if nc_out.isopen():
            nc_out.close()
        if out_file.exists():
            logger.error(f"合并出错，删除不完整的文件：{out_file}")
            out_file.unlink()
        raise
    return out_file


def _compute_solar_cf_year(
    yr: int,
    rsds_file: Path,
    tas_file: Path,
    uas_file: Path,
    vas_file: Path,
    yr_out_file: Path,
    lat_2d: np.ndarray,
    lon_2d: np.ndarray,
    months: str,
    chunk_time: int = 24,
    compress_level: int = 4,
) -> Path:
    """计算单年光伏 CF 并保存到独立 NC 文件。"""
    ds_rsds = xr.open_dataset(rsds_file)
    ds_tas = xr.open_dataset(tas_file)
    ds_uas = xr.open_dataset(uas_file)
    ds_vas = xr.open_dataset(vas_file)

    time_name = find_coord_name(ds_rsds, ["time"])
    rlat_name = find_coord_name(ds_rsds, ["rlat"])
    rlon_name = find_coord_name(ds_rsds, ["rlon"])

    rsds_da = prepare_dataarray(ds_rsds, "rsds", time_name, rlat_name, rlon_name)
    tas_da = prepare_dataarray(ds_tas, "tas", time_name, rlat_name, rlon_name)
    uas_da = prepare_dataarray(ds_uas, "uas", time_name, rlat_name, rlon_name)
    vas_da = prepare_dataarray(ds_vas, "vas", time_name, rlat_name, rlon_name)

    logger.info(f"  rsds time: {ds_rsds[time_name].values[0]} ~ {ds_rsds[time_name].values[-1]}")
    logger.info(f"  tas  time: {ds_tas[time_name].values[0]} ~ {ds_tas[time_name].values[-1]}")

    rsds_units = ds_rsds[rsds_da.name].attrs.get("units", "")
    tas_units = ds_tas[tas_da.name].attrs.get("units", "")
    uas_units = ds_uas[uas_da.name].attrs.get("units", "")
    vas_units = ds_vas[vas_da.name].attrs.get("units", "")

    n_time = ds_rsds.sizes[time_name]
    time_idx, doy_all, hour_all = build_time_index(ds_rsds[time_name], str(yr), months)
    n_selected = len(time_idx)
    logger.info(f"  匹配时间步：{n_selected}/{n_time}")

    nc = create_output_file_rotated(
        out_file=yr_out_file,
        template_file=rsds_file,
        total_time_steps=n_selected,
        compress_level=compress_level,
        chunk_time=chunk_time,
    )
    cf_var = nc.variables["solar_cf"]
    tvar = nc.variables["time"]

    try:
        raw = xr.open_dataset(rsds_file, decode_times=False)
        tvar[:n_selected] = raw["time"].values[time_idx]
        raw.close()

        out_offset = 0
        for start in tqdm(range(0, n_selected, chunk_time), desc=f"solar-{yr}", unit="块"):
            end = min(start + chunk_time, n_selected)
            idx_chunk = time_idx[start:end]

            target_times = ds_rsds[time_name].isel({time_name: idx_chunk}).values

            rsds_raw = rsds_da.isel({time_name: idx_chunk}).values
            tas_raw = interp_instantaneous_to_target_times_chunk(
                tas_da, time_name, target_times, fill_boundary="nearest"
            )
            uas_raw = interp_instantaneous_to_target_times_chunk(
                uas_da, time_name, target_times, fill_boundary="nearest"
            )
            vas_raw = interp_instantaneous_to_target_times_chunk(
                vas_da, time_name, target_times, fill_boundary="nearest"
            )

            rsds_kw = to_kw_m2(rsds_raw, rsds_units)
            tas_c = tas_to_celsius(tas_raw, tas_units)
            uas_ms = wind_to_ms(uas_raw, uas_units)
            vas_ms = wind_to_ms(vas_raw, vas_units)

            cf_chunk = compute_solar_cf_chunk(
                rsds_kw=rsds_kw,
                tas_c=tas_c,
                uas=uas_ms,
                vas=vas_ms,
                doy=doy_all[start:end],
                hour_decimal_utc=hour_all[start:end],
                lats=lat_2d,
                lons=lon_2d,
            )

            cf_var[out_offset : out_offset + (end - start), :, :] = cf_chunk
            out_offset += end - start

            del rsds_raw, tas_raw, uas_raw, vas_raw, rsds_kw, tas_c, uas_ms, vas_ms, cf_chunk
            gc.collect()
    except BaseException:
        if nc.isopen():
            nc.close()
        if yr_out_file.exists():
            logger.error(f"{yr} 计算出错，删除：{yr_out_file}")
            yr_out_file.unlink()
        ds_rsds.close()
        ds_tas.close()
        ds_uas.close()
        ds_vas.close()
        raise

    nc.close()
    ds_rsds.close()
    ds_tas.close()
    ds_uas.close()
    ds_vas.close()
    gc.collect()
    return yr_out_file


# ─────────────────────────────────────────────────────────────────────────────
# 7. 主计算函数
# ─────────────────────────────────────────────────────────────────────────────


def compute_nam12_solar_cf(
    data_dir: str,
    gcm_model: str,
    realization: str,
    rcm_model: str,
    scenario: str,
    years: str,
    months: str,
    output_dir: str,
    chunk_time: int = 24,
    overwrite: bool = False,
    compress_level: int = 4,
    merge: bool = False,
) -> Path:
    """计算 CORDEX NAM-12 光伏容量因子（逐年独立计算 + 合并）。"""
    y0, y1 = parse_years(years)

    rsds_files = cordex_var_files(data_dir, gcm_model, realization, rcm_model, scenario, "rsds", (y0, y1))
    tas_files = cordex_var_files(data_dir, gcm_model, realization, rcm_model, scenario, "tas", (y0, y1))
    uas_files = cordex_var_files(data_dir, gcm_model, realization, rcm_model, scenario, "uas", (y0, y1))
    vas_files = cordex_var_files(data_dir, gcm_model, realization, rcm_model, scenario, "vas", (y0, y1))

    common_years = sorted(set(rsds_files) & set(tas_files) & set(uas_files) & set(vas_files))
    if not common_years:
        raise FileNotFoundError(f"四个变量没有共同年份（{y0}-{y1}）")

    for name, fdict in [("rsds", rsds_files), ("tas", tas_files), ("uas", uas_files), ("vas", vas_files)]:
        missing = set(common_years) - set(fdict)
        if missing:
            logger.warning(f"{name} 缺少年份：{sorted(missing)}")

    logger.info(f"待处理年份：{common_years[0]}-{common_years[-1]}，共 {len(common_years)} 年")

    out_file = (
        Path(output_dir)
        / gcm_model
        / realization
        / scenario
        / f"solar_CF_NAM-12_{gcm_model}_{realization}_{rcm_model}_{scenario}_{years}_{months_tag(months)}.nc"
    )
    if out_file.exists() and not overwrite:
        logger.info(f"已存在，跳过：{out_file}")
        return out_file

    # 年度文件目录
    yearly_dir = out_file.parent / "yearly"
    yearly_dir.mkdir(parents=True, exist_ok=True)

    # 加载网格坐标
    ds_tpl = xr.open_dataset(rsds_files[common_years[0]])
    lat_2d = ds_tpl["lat"].values.astype(np.float32)
    lon_2d = ds_tpl["lon"].values.astype(np.float32)
    ds_tpl.close()
    logger.info(f"网格大小：lat_2d {lat_2d.shape}, lon_2d {lon_2d.shape}")

    # 逐年计算
    yearly_files: list[Path] = []
    for yr in common_years:
        yr_file = (
            yearly_dir
            / f"solar_CF_NAM-12_{gcm_model}_{realization}_{rcm_model}_{scenario}_{yr}_{months_tag(months)}.nc"
        )
        if yr_file.exists() and not overwrite:
            logger.info(f"  {yr} 已存在，跳过")
            yearly_files.append(yr_file)
            continue

        logger.info(f"处理 {yr}...")
        try:
            _compute_solar_cf_year(
                yr=yr,
                rsds_file=rsds_files[yr],
                tas_file=tas_files[yr],
                uas_file=uas_files[yr],
                vas_file=vas_files[yr],
                yr_out_file=yr_file,
                lat_2d=lat_2d,
                lon_2d=lon_2d,
                months=months,
                chunk_time=chunk_time,
                compress_level=compress_level,
            )
        except Exception:
            if yr_file.exists():
                logger.error(f"{yr} 计算出错，删除：{yr_file}")
                yr_file.unlink()
            raise
        yearly_files.append(yr_file)

    if merge:
        # 合并年度文件
        logger.info(f"开始合并 {len(yearly_files)} 个年度文件...")
        _merge_yearly_nc_files(yearly_files, out_file, "solar_cf", compress_level, chunk_time)
        logger.info(f"已合并保存：{out_file}")
        logger.info(f"可删除年度文件目录：{yearly_dir}")
    else:
        logger.info(f"逐年文件已保存到：{yearly_dir}")
        logger.info(f"未执行合并，单年文件保留。")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CMIP6-CORDEX NAM-12 光伏容量因子计算（Erbs + 双轴跟踪 + 论文温度修正）"
    )
    parser.add_argument("--data_dir", default="data/CORDEX-CMIP6/NAM-12/1hr", help="数据根目录")
    parser.add_argument("--gcm_model", required=True, help="GCM 模型，例如 MPI-ESM1-2-LR")
    parser.add_argument("--realization", required=True, help="实现编号，例如 r1i1p1f1")
    parser.add_argument("--rcm_model", default="CRCM5", help="RCM 模型，默认 CRCM5")
    parser.add_argument("--scenario", required=True, help="情景，例如 ssp126")
    parser.add_argument("--years", default="2020-2060", help="年段 YYYY-YYYY")
    parser.add_argument("--months", default="", help="逗号分隔月份；默认全年")
    parser.add_argument("--output_dir", default="output/CFs_of_solar_NAM-12", help="输出目录")
    parser.add_argument("--chunk_time", type=int, default=24, help="每块时间步数")
    parser.add_argument("--compress_level", type=int, default=4, help="NetCDF 压缩级别 0-9")
    parser.add_argument("--overwrite", action="store_true", help="覆盖已有输出")
    parser.add_argument("--merge", action="store_true", help="执行合并, 默认不合并, 逐年输出")
    args = parser.parse_args()

    compute_nam12_solar_cf(
        data_dir=args.data_dir,
        gcm_model=args.gcm_model,
        realization=args.realization,
        rcm_model=args.rcm_model,
        scenario=args.scenario,
        years=args.years,
        months=args.months,
        output_dir=args.output_dir,
        chunk_time=args.chunk_time,
        overwrite=args.overwrite,
        compress_level=args.compress_level,
        merge=args.merge,
    )


if __name__ == "__main__":
    main()

"""

python S04E01_Simulate_Solar_CF_NAM-12.py \\
  --gcm_model MPI-ESM1-2-LR --realization r1i1p1f1 \\
  --scenario ssp126 --years 2020-2060

"""
