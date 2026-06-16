"""
全球逐时光伏容量因子计算（ERA5-Land ssrd/t2m/u10/v10）
============================================================

本脚本用于从 ERA5-Land 月度逐小时数据计算全球光伏容量因子。
计算方法与 S01E01_Simulate_Solar_CF_BCSD.py 保持一致：
  1. 将 ERA5-Land ssrd 的日内逐小时累计量（J m-2）转换为逐小时 GHI（kW m-2）；
  2. t2m 转为摄氏度；
  3. u10/v10 合成为近地面风速，用于光伏组件温度修正；
  4. 根据 UTC 时间、纬度、经度计算太阳天顶角；
  5. 使用 Erbs 模型将 GHI 分解为直接辐射和散射辐射；
  6. 使用与 BCSD 脚本一致的双轴跟踪近似：直射分量按 DNI 入射，散射分量采用 isotropic diffuse；
  7. 使用论文中的组件温度修正和系统效率系数 SYS_COEF = 0.8056；
  8. 按时间分块计算，并逐块写入 NetCDF，避免一次性保存全月数组。

ERA5-Land 累积变量说明
----------------------
ERA5-Land 的 ssrd 是日内逐小时累计的能量，单位为 J m-2。
本脚本按 [d 01:00, (d+1) 00:00] 的日内累积周期处理：
  - 若时间为 01:00，则该小时增量等于当前累计值；
  - 其他小时增量等于当前累计值减去上一时刻累计值；
  - 月初 00:00 需要使用前一个月文件最后一个时刻作为差分起点。
若找不到前一个月文件，默认把第一个无法差分的小时设为 0，并给出警告。

输入数据路径格式
----------------
  {data_dir}/ERA5_land/global/ssrd/ssrd_YYYY_MM.nc
  {data_dir}/ERA5_land/global/t2m/t2m_YYYY_MM.nc
  {data_dir}/ERA5_land/global/u10/u10_YYYY_MM.nc
  {data_dir}/ERA5_land/global/v10/v10_YYYY_MM.nc

输出
----
  {output_dir}/solar_cf_YYYY_MM.nc

输出变量
--------
  solar_cf(time, lat, lon), 单位 1, 数值范围 [0, 1]

示例
----
python S02E01_Simulate_Solar_CF_ERA5Land.py \
  --data_dir data \
  --years 2020 \
  --months 1 \
  --output_dir output/CFs_of_solar \
  --chunk_time 24
"""

from __future__ import annotations

import argparse
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

# ── 与 BCSD 脚本一致的 PV 参数 ───────────────────────────────────────────────
SYS_COEF = 0.8056  # system performance coefficient = 1 - 19.44% grid/system loss
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


def iter_year_months(years: str, months: str) -> list[tuple[int, int]]:
    """生成待处理的 (year, month) 列表。"""
    y0, y1 = parse_years(years)
    month_list = parse_months(months)
    return [(y, m) for y in range(y0, y1 + 1) for m in month_list]


def previous_year_month(year: int, month: int) -> tuple[int, int]:
    """返回前一个年月。"""
    if month == 1:
        return year - 1, 12
    return year, month - 1


def era5land_file(data_dir: str | Path, var: str, year: int, month: int) -> Path:
    """构建 ERA5-Land 月文件路径。"""
    root = Path(data_dir) / "ERA5_land" / "global"
    return root / var / f"{var}_{year:04d}_{month:02d}.nc"


def get_var_name(ds: xr.Dataset, preferred: str) -> str:
    """获取变量名。优先使用 preferred，否则在唯一 data_var 时使用该变量。"""
    if preferred in ds.data_vars:
        return preferred
    if len(ds.data_vars) == 1:
        return list(ds.data_vars)[0]
    raise KeyError(f"无法在数据集中确定变量 {preferred}，可用变量：{list(ds.data_vars)}")


def find_coord_name(ds: xr.Dataset, candidates: Iterable[str]) -> str:
    """兼容 time/valid_time、lat/latitude、lon/longitude 等命名。"""
    for c in candidates:
        if c in ds.coords or c in ds.dims:
            return c
    raise KeyError(f"找不到坐标名，候选：{list(candidates)}；数据维度：{list(ds.dims)}")


def build_land_mask(lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    """根据输出网格 lat/lon 生成陆地掩膜 (n_lat, n_lon)，True=陆地。

    global_land_mask 要求经度 ∈ [-180, 180]、纬度 ∈ [-90, 90]，
    因此对 ERA5-Land [0, 360) 经度做归一化，并夹紧纬度避免浮点越界。
    """
    lats = np.asarray(lats, dtype=np.float64)
    lons = np.asarray(lons, dtype=np.float64)

    lons_conv = ((lons + 180.0) % 360.0) - 180.0  # [0,360) -> [-180,180)
    lats_clip = np.clip(lats, -90.0, 90.0)

    lon2d, lat2d = np.meshgrid(lons_conv, lats_clip)  # (n_lat, n_lon)
    return globe.is_land(lat2d, lon2d).astype(bool)


def prepare_dataarray(
    ds: xr.Dataset,
    preferred_var: str,
    time_name: str,
    lat_name: str,
    lon_name: str,
) -> xr.DataArray:
    """获取变量，并压缩长度为 1 的额外维度，然后转为 (time, lat, lon)。"""
    var = get_var_name(ds, preferred_var)
    da = ds[var]
    for dim in list(da.dims):
        if dim not in {time_name, lat_name, lon_name}:
            if da.sizes[dim] == 1:
                da = da.isel({dim: 0}, drop=True)
            else:
                raise ValueError(f"变量 {var} 存在非单例额外维度 {dim}={da.sizes[dim]}，无法自动处理。")
    return da.transpose(time_name, lat_name, lon_name)


def rename_coords_to_match(
    ds: xr.Dataset,
    target_time: str,
    target_lat: str,
    target_lon: str,
) -> xr.Dataset:
    """将另一个 ERA5-Land 数据集的坐标名改成与 ssrd 数据集一致。"""
    src_time = find_coord_name(ds, ["time", "valid_time"])
    src_lat = find_coord_name(ds, ["lat", "latitude"])
    src_lon = find_coord_name(ds, ["lon", "longitude"])

    rename_map: dict[str, str] = {}
    if src_time != target_time:
        rename_map[src_time] = target_time
    if src_lat != target_lat:
        rename_map[src_lat] = target_lat
    if src_lon != target_lon:
        rename_map[src_lon] = target_lon
    if rename_map:
        ds = ds.rename(rename_map)
    return ds


def validate_same_grid(
    ref: xr.Dataset, others: dict[str, xr.Dataset], time_name: str, lat_name: str, lon_name: str
) -> None:
    """检查输入变量是否在同一时间和空间网格上。"""
    ref_time = ref[time_name].values
    ref_lat = ref[lat_name].values
    ref_lon = ref[lon_name].values
    for name, ds in others.items():
        if not np.array_equal(ref_time, ds[time_name].values):
            raise ValueError(f"{name} 与 ssrd 的时间坐标不一致。ERA5-Land 月文件应逐小时完全一致。")
        if not np.allclose(ref_lat, ds[lat_name].values, rtol=0.0, atol=1e-6, equal_nan=True):
            raise ValueError(f"{name} 与 ssrd 的纬度网格不一致。")
        if not np.allclose(ref_lon, ds[lon_name].values, rtol=0.0, atol=1e-6, equal_nan=True):
            raise ValueError(f"{name} 与 ssrd 的经度网格不一致。")


# ─────────────────────────────────────────────────────────────────────────────
# 1. 单位转换与 ERA5-Land 累积量差分
# ─────────────────────────────────────────────────────────────────────────────


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
    """把风速分量转换为 m s-1。ERA5-Land u10/v10 通常已经是 m s-1。"""
    x = arr.astype(np.float32, copy=False)
    units_l = (units or "").lower().replace(" ", "")
    if "km/h" in units_l or "kmh-1" in units_l or "kmhr-1" in units_l:
        x = x / 3.6
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    return x.astype(np.float32)


def load_previous_accum_boundary(
    data_dir: str | Path,
    var: str,
    year: int,
    month: int,
    lat_values: np.ndarray,
    lon_values: np.ndarray,
    *,
    time_name: str,
    lat_name: str,
    lon_name: str,
) -> np.ndarray | None:
    """读取前一个月同一变量最后一个时刻，作为月初 00:00 差分起点。"""
    py, pm = previous_year_month(year, month)
    prev_file = era5land_file(data_dir, var, py, pm)
    if not prev_file.exists():
        logger.warning(f"找不到前一个月 {var} 文件：{prev_file}。月初第一个无法差分的小时将按策略处理。")
        return None

    with xr.open_dataset(prev_file) as prev_ds0:
        prev_ds = rename_coords_to_match(prev_ds0, time_name, lat_name, lon_name)
        prev_da = prepare_dataarray(prev_ds, var, time_name, lat_name, lon_name)
        if prev_da.sizes[lat_name] != len(lat_values) or prev_da.sizes[lon_name] != len(lon_values):
            raise ValueError(f"前一个月 {var} 文件空间尺寸与当前文件不一致：{prev_file}")
        if not np.allclose(prev_ds[lat_name].values, lat_values, rtol=0.0, atol=1e-6, equal_nan=True):
            # raise ValueError(f"前一个月 {var} 文件纬度与当前文件不一致：{prev_file}")
            print(f"警告：前一个月 {var} 文件纬度与当前文件不完全一致，已继续处理：{prev_file}")
        if not np.allclose(prev_ds[lon_name].values, lon_values, rtol=0.0, atol=1e-6, equal_nan=True):
            # raise ValueError(f"前一个月 {var} 文件经度与当前文件不一致：{prev_file}")
            print(f"警告：前一个月 {var} 文件经度与当前文件不完全一致，已继续处理：{prev_file}")

        return prev_da.isel({time_name: -1}).values.astype(np.float32)


def accumulated_to_hourly_increment_chunk(
    da: xr.DataArray,
    time_name: str,
    start: int,
    end: int,
    prev_boundary: np.ndarray | None,
    *,
    missing_boundary_policy: str = "zero",
) -> np.ndarray:
    """
    将 ERA5-Land 日内累计变量转换为逐小时增量。

    对目标时间 [start, end) 中的每个时刻 t：
      - 若 t.hour == 1，则增量 = 当前累计值；
      - 否则增量 = 当前累计值 - 上一时刻累计值；
      - 当 start == 0 且第一个时刻不是 01:00 时，上一时刻来自前一个月最后一小时。
    """
    if end <= start:
        raise ValueError(f"非法时间块：start={start}, end={end}")

    cur = da.isel({time_name: slice(start, end)}).values.astype(np.float32)
    time_da = da[time_name].isel({time_name: slice(start, end)})
    hours = time_da.dt.hour.values.astype(np.int16)

    if start > 0:
        prev0 = da.isel({time_name: start - 1}).values.astype(np.float32)
        has_prev0 = True
    else:
        prev0 = prev_boundary
        has_prev0 = prev_boundary is not None

    inc = np.empty_like(cur, dtype=np.float32)

    if has_prev0:
        diff0 = cur[0] - prev0
        inc[0] = cur[0] if hours[0] == 1 else diff0
    else:
        if hours[0] == 1:
            inc[0] = cur[0]
        elif missing_boundary_policy == "zero":
            logger.warning(f"{str(time_da.values[0])} 缺少前一小时累计值，已将该小时增量设为 0。")
            inc[0] = 0.0
        elif missing_boundary_policy == "nan":
            logger.warning(f"{str(time_da.values[0])} 缺少前一小时累计值，已将该小时增量设为 NaN。")
            inc[0] = np.nan
        elif missing_boundary_policy == "current":
            logger.warning(f"{str(time_da.values[0])} 缺少前一小时累计值，已直接使用当前累计值。")
            inc[0] = cur[0]
        else:
            raise ValueError(f"不支持 missing_boundary_policy={missing_boundary_policy!r}")

    if cur.shape[0] > 1:
        diff = cur[1:] - cur[:-1]
        reset = hours[1:] == 1
        inc[1:] = np.where(reset[:, None, None], cur[1:], diff)

    # 防御性处理：日内 reset 或数据异常导致负差分时，用当前累计值；剩余负值截断为 0。
    negative = inc < 0
    if np.any(negative):
        inc = np.where(negative, cur, inc)
    inc = np.maximum(inc, 0.0).astype(np.float32)
    return inc


def ssrd_increment_to_ghi_kw(inc_j_m2: np.ndarray) -> np.ndarray:
    """ssrd 每小时能量增量 J m-2 -> GHI kW m-2。"""
    ghi = inc_j_m2 / np.float32(3_600_000.0)
    ghi = np.nan_to_num(ghi, nan=0.0, posinf=0.0, neginf=0.0)
    return np.maximum(ghi, 0.0).astype(np.float32)


def tp_increment_to_mm_h(inc_m: np.ndarray) -> np.ndarray:
    """tp 每小时水深增量 m -> mm h-1。当前光伏计算不使用该函数。"""
    tp_mm_h = inc_m * np.float32(1000.0)
    tp_mm_h = np.nan_to_num(tp_mm_h, nan=0.0, posinf=0.0, neginf=0.0)
    return np.maximum(tp_mm_h, 0.0).astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# 2. 光伏容量因子计算核心：与 BCSD 脚本保持一致
# ─────────────────────────────────────────────────────────────────────────────


def _erbs_diffuse_fraction(kt: np.ndarray) -> np.ndarray:
    """Erbs 经验模型：由晴空指数 kt 估计漫射比例。"""
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


def solar_cos_zenith(
    doy: np.ndarray,
    hour_decimal_utc: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
) -> np.ndarray:
    """
    计算太阳天顶角余弦 cos(SZA)，shape = (time, lat, lon)。
    这里显式加入经度：true solar time = UTC minutes + equation_of_time + 4 * longitude。
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

    true_solar_minutes = hour_decimal_utc[:, None, None] * 60.0 + eot[:, None, None] + 4.0 * lons[None, None, :]
    hour_angle = np.deg2rad(true_solar_minutes / 4.0 - 180.0)

    lat_rad = np.deg2rad(lats)
    cos_sza = np.sin(lat_rad)[None, :, None] * np.sin(decl)[:, None, None] + np.cos(lat_rad)[None, :, None] * np.cos(
        decl
    )[:, None, None] * np.cos(hour_angle)
    return np.maximum(cos_sza, 0.0).astype(np.float32)


def extraterrestrial_horizontal_irradiance(doy: np.ndarray, cos_sza: np.ndarray) -> np.ndarray:
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
    """GHI -> 漫射比例，并返回 cos(SZA)。"""
    cos_sza = solar_cos_zenith(doy, hour_decimal_utc, lats, lons)
    etr = extraterrestrial_horizontal_irradiance(doy, cos_sza)

    with np.errstate(invalid="ignore", divide="ignore"):
        kt = np.where(etr > 1e-6, ghi_kw / etr, 0.0)
    kt = np.clip(kt, 0.0, 1.0)
    fd = _erbs_diffuse_fraction(kt)
    return fd, cos_sza


def plane_irradiance_two_axis_erbs(ghi_kw: np.ndarray, fd: np.ndarray, cos_sza: np.ndarray) -> np.ndarray:
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
    """
    论文温度修正：
      TEMcoef = 1 - gamma * (Tcell - T_STC)
      Tcell = c1 + c2*Tair + c3*I - c4*V
    其中 I 使用 W m-2，V 使用 m s-1。
    """
    irradiance_w = plane_irradiance_kw * 1000.0
    t_cell = C1 + C2 * tas_c + C3 * irradiance_w - C4 * wind_speed
    tem_coef = 1.0 - GAMMA_TEMP * (t_cell - T_STC)
    tem_coef = np.nan_to_num(tem_coef, nan=0.0, posinf=0.0, neginf=0.0)
    return np.maximum(tem_coef, 0.0).astype(np.float32)


def compute_solar_cf_chunk(
    ghi_kw: np.ndarray,
    tas_c: np.ndarray,
    u10: np.ndarray,
    v10: np.ndarray,
    doy: np.ndarray,
    hour_decimal_utc: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
) -> np.ndarray:
    """计算一个时间块的光伏容量因子。"""
    fd, cos_sza = compute_diffuse_fraction_erbs(ghi_kw, doy, hour_decimal_utc, lats, lons)
    plane_kw = plane_irradiance_two_axis_erbs(ghi_kw, fd, cos_sza)

    wind_speed = np.sqrt(u10.astype(np.float32) ** 2 + v10.astype(np.float32) ** 2)
    temp_coef = paper_temperature_coefficient(plane_kw, tas_c, wind_speed)

    cf = plane_kw * temp_coef * SYS_COEF
    cf = np.clip(cf, 0.0, 1.0)
    return cf.astype(np.float32)


def time_features(time_da: xr.DataArray) -> tuple[np.ndarray, np.ndarray]:
    """由 xarray 时间坐标得到 doy 和 UTC 小时小数。"""
    doy = time_da.dt.dayofyear.values.astype(np.float32)
    hour = time_da.dt.hour.values.astype(np.float32)
    minute = time_da.dt.minute.values.astype(np.float32)
    second = time_da.dt.second.values.astype(np.float32)
    hour_decimal = hour + minute / 60.0 + second / 3600.0
    return doy, hour_decimal


# ─────────────────────────────────────────────────────────────────────────────
# 3. NetCDF 输出
# ─────────────────────────────────────────────────────────────────────────────


def create_output_file(
    out_file: Path,
    template_file: Path,
    n_time: int,
    lat_values: np.ndarray,
    lon_values: np.ndarray,
    time_name: str,
    lat_name: str,
    lon_name: str,
    chunk_time: int,
    attrs: dict,
    compress_level: int = 4,
) -> Dataset:
    """创建 NetCDF 输出文件，并写入坐标。输出维度统一为 time/lat/lon。"""
    out_file.parent.mkdir(parents=True, exist_ok=True)

    raw = xr.open_dataset(template_file, decode_times=False)
    raw_time = raw[time_name].values[:n_time]
    time_attrs = dict(raw[time_name].attrs)
    lat_attrs = dict(raw[lat_name].attrs) if lat_name in raw else {}
    lon_attrs = dict(raw[lon_name].attrs) if lon_name in raw else {}
    raw.close()

    nc = Dataset(out_file, "w", format="NETCDF4")
    nc.createDimension("time", n_time)
    nc.createDimension("lat", len(lat_values))
    nc.createDimension("lon", len(lon_values))

    tvar = nc.createVariable("time", raw_time.dtype, ("time",))
    tvar[:] = raw_time
    for k, v in time_attrs.items():
        try:
            setattr(tvar, k, v)
        except Exception:
            pass
    tvar.standard_name = "time"

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

    chunksizes = (min(chunk_time, n_time), len(lat_values), len(lon_values))
    cf = nc.createVariable(
        "solar_cf",
        "f4",
        ("time", "lat", "lon"),
        zlib=True,
        complevel=compress_level,
        chunksizes=chunksizes,
        fill_value=np.float32(np.nan),
    )
    cf.long_name = "Solar photovoltaic capacity factor"
    cf.units = "1"
    cf.description = (
        "Solar CF computed from ERA5-Land ssrd/t2m/u10/v10. "
        "ssrd cumulative J m-2 is converted to hourly GHI; "
        "CF uses Erbs diffuse fraction, two-axis tracking approximation, "
        "paper temperature correction, and SYS_COEF=0.8056."
    )

    for k, v in attrs.items():
        try:
            setattr(nc, k, v)
        except Exception:
            pass

    return nc


# ─────────────────────────────────────────────────────────────────────────────
# 4. 主计算函数
# ─────────────────────────────────────────────────────────────────────────────


def compute_month_solar_cf(
    data_dir: str,
    year: int,
    month: int,
    output_dir: str,
    chunk_time: int = 24,
    overwrite: bool = False,
    compress_level: int = 4,
    missing_boundary_policy: str = "zero",
) -> Path:
    """计算单个月份的 ERA5-Land 光伏容量因子。"""
    tag = f"{year:04d}_{month:02d}"
    ssrd_file = era5land_file(data_dir, "ssrd", year, month)
    t2m_file = era5land_file(data_dir, "t2m", year, month)
    u10_file = era5land_file(data_dir, "u10", year, month)
    v10_file = era5land_file(data_dir, "v10", year, month)

    for p in [ssrd_file, t2m_file, u10_file, v10_file]:
        if not p.exists():
            raise FileNotFoundError(f"找不到输入文件：{p}")

    out_file = Path(output_dir) / f"solar_cf_{tag}.nc"
    tmp_file = out_file.with_name(f".{out_file.name}.tmp")
    if out_file.exists() and not overwrite:
        logger.info(f"✓ 已存在，跳过：{out_file}")
        return out_file
    if tmp_file.exists():
        tmp_file.unlink()

    logger.info(f"ssrd: {ssrd_file}")
    logger.info(f"t2m : {t2m_file}")
    logger.info(f"u10 : {u10_file}")
    logger.info(f"v10 : {v10_file}")

    ds_ssrd = xr.open_dataset(ssrd_file)
    ds_t2m = xr.open_dataset(t2m_file)
    ds_u10 = xr.open_dataset(u10_file)
    ds_v10 = xr.open_dataset(v10_file)

    time_name = find_coord_name(ds_ssrd, ["valid_time", "time"])
    lat_name = find_coord_name(ds_ssrd, ["latitude", "lat"])
    lon_name = find_coord_name(ds_ssrd, ["longitude", "lon"])

    ds_t2m = rename_coords_to_match(ds_t2m, time_name, lat_name, lon_name)
    ds_u10 = rename_coords_to_match(ds_u10, time_name, lat_name, lon_name)
    ds_v10 = rename_coords_to_match(ds_v10, time_name, lat_name, lon_name)
    validate_same_grid(ds_ssrd, {"t2m": ds_t2m, "u10": ds_u10, "v10": ds_v10}, time_name, lat_name, lon_name)

    ssrd_da = prepare_dataarray(ds_ssrd, "ssrd", time_name, lat_name, lon_name)
    t2m_da = prepare_dataarray(ds_t2m, "t2m", time_name, lat_name, lon_name)
    u10_da = prepare_dataarray(ds_u10, "u10", time_name, lat_name, lon_name)
    v10_da = prepare_dataarray(ds_v10, "v10", time_name, lat_name, lon_name)

    n_time = ssrd_da.sizes[time_name]
    lats = ds_ssrd[lat_name].values.astype(np.float32)
    lons = ds_ssrd[lon_name].values.astype(np.float32)
    logger.info(f"维度：time={n_time}, lat={len(lats)}, lon={len(lons)}")
    logger.info(f"输出坐标统一为：time/lat/lon；变量名：solar_cf")

    land_mask = build_land_mask(lats, lons)  # (n_lat, n_lon) bool
    logger.info(f"陆地格点占比：{100.0 * land_mask.mean():.1f}%")

    prev_ssrd = load_previous_accum_boundary(
        data_dir,
        "ssrd",
        year,
        month,
        lats,
        lons,
        time_name=time_name,
        lat_name=lat_name,
        lon_name=lon_name,
    )

    attrs = {
        "source": "ERA5-Land",
        "year": year,
        "month": month,
        "source_files": ", ".join(str(p) for p in [ssrd_file, t2m_file, u10_file, v10_file]),
        "ssrd_conversion": "ERA5-Land accumulated ssrd J m-2 -> hourly increment -> kW m-2; daily reset at 01:00 UTC",
        "previous_month_boundary_used": bool(prev_ssrd is not None),
        "missing_boundary_policy": missing_boundary_policy,
        "pv_method": "Erbs diffuse fraction + two-axis tracking approximation + paper temperature correction",
        "SYS_COEF": SYS_COEF,
        "GAMMA_TEMP": GAMMA_TEMP,
        "temperature_model": f"Tcell={C1}+{C2}*Tair+{C3}*I_Wm2-{C4}*wind_ms",
        "ocean_handling": "ocean grid cells set to NaN via global-land-mask",
    }

    nc = create_output_file(
        out_file=tmp_file,
        template_file=ssrd_file,
        n_time=n_time,
        lat_values=lats,
        lon_values=lons,
        time_name=time_name,
        lat_name=lat_name,
        lon_name=lon_name,
        chunk_time=chunk_time,
        attrs=attrs,
        compress_level=compress_level,
    )
    cf_var = nc.variables["solar_cf"]

    t2m_units = t2m_da.attrs.get("units", "")
    u10_units = u10_da.attrs.get("units", "")
    v10_units = v10_da.attrs.get("units", "")

    total_sum = 0.0
    total_count = 0
    positive_sum = 0.0
    positive_count = 0
    max_cf = 0.0

    try:
        logger.info(f"开始分块计算，chunk_time={chunk_time}")
        for start in tqdm(range(0, n_time, chunk_time), desc=f"solar-{tag}", unit="块"):
            end = min(start + chunk_time, n_time)

            ssrd_inc = accumulated_to_hourly_increment_chunk(
                ssrd_da,
                time_name,
                start,
                end,
                prev_ssrd,
                missing_boundary_policy=missing_boundary_policy,
            )
            ghi_kw = ssrd_increment_to_ghi_kw(ssrd_inc)
            del ssrd_inc

            t2m_raw = t2m_da.isel({time_name: slice(start, end)}).values
            u10_raw = u10_da.isel({time_name: slice(start, end)}).values
            v10_raw = v10_da.isel({time_name: slice(start, end)}).values

            tas_c = tas_to_celsius(t2m_raw, t2m_units)
            u10_ms = wind_to_ms(u10_raw, u10_units)
            v10_ms = wind_to_ms(v10_raw, v10_units)
            del t2m_raw, u10_raw, v10_raw

            time_chunk_da = ssrd_da[time_name].isel({time_name: slice(start, end)})
            doy, hour_decimal = time_features(time_chunk_da)

            cf_chunk = compute_solar_cf_chunk(
                ghi_kw=ghi_kw,
                tas_c=tas_c,
                u10=u10_ms,
                v10=v10_ms,
                doy=doy,
                hour_decimal_utc=hour_decimal,
                lats=lats,
                lons=lons,
            )

            # 海洋格点恢复为 NaN（修复 nan_to_num 把海洋填 0 的问题）
            cf_chunk = np.where(land_mask[None, :, :], cf_chunk, np.nan).astype(np.float32)

            cf_var[start:end, :, :] = cf_chunk

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

            del ghi_kw, tas_c, u10_ms, v10_ms, cf_chunk
            gc.collect()
    finally:
        nc.close()
        ds_ssrd.close()
        ds_t2m.close()
        ds_u10.close()
        ds_v10.close()

    os.replace(tmp_file, out_file)
    logger.info(f"✓ 已保存：{out_file}")
    if total_count:
        logger.info(f"  整体平均 CF   : {total_sum / total_count:.4f}")
        logger.info(f"  全局最大 CF   : {max_cf:.4f}")
        logger.info(f"  出力时段占比  : {100.0 * positive_count / total_count:.2f}%")
    if positive_count:
        logger.info(f"  出力时段平均 CF: {positive_sum / positive_count:.4f}")
    return out_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ERA5-Land ssrd/t2m/u10/v10 光伏容量因子计算（ssrd 累积量差分 + Erbs + 温度修正）"
    )
    parser.add_argument("--data_dir", default="data", help="数据根目录，例如 data")
    parser.add_argument("--years", default="2015-2024", help="年段 YYYY-YYYY 或 YYYY，包含两端")
    parser.add_argument("--months", default="", help="逗号分隔月份；默认全年 1-12")
    parser.add_argument("--output_dir", default="output/CFs_of_solar_ERA5Land", help="输出目录")
    parser.add_argument("--chunk_time", type=int, default=24, help="每块处理的小时数；全球 0.1° 数据建议 12-48")
    parser.add_argument("--compress_level", type=int, default=4, help="NetCDF 压缩级别 0-9")
    parser.add_argument("--overwrite", action="store_true", help="若输出文件已存在，则覆盖重算")
    parser.add_argument(
        "--missing_boundary_policy",
        default="zero",
        choices=["zero", "nan", "current"],
        help="月初缺少前一个月文件时，第一个无法差分小时的处理方式。默认 zero。",
    )
    args = parser.parse_args()

    year_months = iter_year_months(args.years, args.months)
    logger.info(f"共 {len(year_months)} 个年月待处理：{year_months[0]} -> {year_months[-1]}")
    for i, (year, month) in enumerate(year_months, 1):
        logger.info(f"\n{'=' * 80}")
        logger.info(f"[{i}/{len(year_months)}] 处理 {year:04d}-{month:02d}")
        logger.info(f"{'=' * 80}")
        compute_month_solar_cf(
            data_dir=args.data_dir,
            year=year,
            month=month,
            output_dir=args.output_dir,
            chunk_time=args.chunk_time,
            overwrite=args.overwrite,
            compress_level=args.compress_level,
            missing_boundary_policy=args.missing_boundary_policy,
        )

    logger.info("\n✓ 全部完成！")


if __name__ == "__main__":
    main()

"""
使用方法示例：
python S02E01_Simulate_Solar_CF_ERA5Land.py --data_dir data --years 2015-2025 
"""
