"""
区域逐时/逐3小时光伏容量因子计算（BCSD rsds + tas + uas/vas）
=============================================================

本脚本用于基于 BCSD 降尺度气候数据计算区域光伏容量因子。

计算逻辑：
  1. 从 BCSD 数据中读取 rsds、tas、uas、vas 四个变量；
  2. rsds 作为全球水平短波辐射 GHI，单位自动转换为 kW m-2；
  3. tas 自动转换为摄氏度；
  4. 以 rsds 的半点时间轴为主，将瞬时变量 tas/uas/vas 线性插值到 rsds 时间；
  5. uas 和 vas 合成为近地面风速，用于光伏组件温度修正；
  6. 根据 UTC 时间、纬度和经度计算太阳天顶角；
  7. 使用 Erbs 模型将 GHI 分解为直接辐射和散射辐射；
  8. 假设光伏组件采用双轴跟踪：
       - 直射辐射按 DNI 入射到组件平面；
       - 散射辐射采用 isotropic diffuse 近似；
  9. 使用论文中的温度修正方法计算组件温度和温度系数；
 10. 使用论文中的系统效率系数 SYS_coef = 0.8056；
 11. 计算光伏容量因子，并限制在 [0, 1] 范围内；
 12. 按时间分块计算，并逐块写入 NetCDF 文件，避免一次性占用过多内存。

输入数据路径格式：
  {data_dir}/{model}/{region}/{model}/*.nc

例如：
  data/bcsd_outputs/MIROC-ES2H/Austria/MIROC-ES2H/
    rsds_3h_bcsd_on_0p1deg_Austria_MIROC-ES2H_ssp126_2015-2060.nc
    tas_3h_bcsd_on_0p1deg_Austria_MIROC-ES2H_ssp126_2015-2060.nc
    uas_3h_bcsd_on_0p1deg_Austria_MIROC-ES2H_ssp126_2015-2060.nc
    vas_3h_bcsd_on_0p1deg_Austria_MIROC-ES2H_ssp126_2015-2060.nc

输入变量要求：
  - rsds 文件中包含 rsds_bcsd 或唯一数据变量；
  - tas 文件中包含 tas_bcsd 或唯一数据变量；
  - uas 文件中包含 uas_bcsd 或唯一数据变量；
  - vas 文件中包含 vas_bcsd 或唯一数据变量；
  - 数据维度通常为 (time, lat, lon)。

输出：
  {output_dir}/{model}/{region}/solar_CF_{region}_{model}_{scenario}_{years}_{months}.nc

输出变量：
  - solar_cf(time, lat, lon)
  - 单位：1
  - 数值范围：[0, 1]

注意：
  - 本脚本不把 3h 数据插值到小时；输出时间分辨率与 rsds 输入一致。
  - 时间对齐采用方案4：保留 rsds/pr 的半点时间轴，将瞬时 tas/uas/vas 插值到 rsds 时间。
  - rsds 通常是 W m-2 的瞬时或平均通量，不是 ERA5-Land ssrd 的累计 J m-2。
  - 若输入变量没有单位属性，脚本会根据数值大小进行简单判断。
  - region 可指定单个区域，也可设置为 all 批量处理 model 下所有有效区域。
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

# ── 论文中的 PV 参数 ─────────────────────────────────────────────────────────
SYS_COEF = 0.8056  # system performance coefficient
T_STC = 25.0  # °C
GAMMA_TEMP = 0.005  # °C-1, monocrystalline silicon temperature coefficient
C1 = 4.3  # °C
C2 = 0.943
C3 = 0.028  # °C m2 W-1
C4 = 1.528  # °C s m-1


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


def _parse_year_range_from_filename(name: str) -> tuple[int | None, int | None]:
    """从文件名解析年份范围，如 ..._ssp126_2015-2060.nc -> (2015, 2060)，
    单年 ..._ssp126_2015.nc -> (2015, 2015)。无法解析返回 (None, None)。"""
    m = re.search(r"_(\d{4})-(\d{4})\.nc$", name)
    if m:
        return int(m.group(1)), int(m.group(2))
    m = re.search(r"_(\d{4})\.nc$", name)
    if m:
        y = int(m.group(1))
        return y, y
    return None, None


def _select_best_bcsd_file(files: list[str], years: str | None) -> Path:
    """从多个匹配文件中选出最合适的一个。

    选择策略：
      1. 优先选择文件名年份范围“完全覆盖”请求 years、且范围最紧的文件；
      2. 否则选择文件名年份范围最大（覆盖最广）的文件；
      3. 若文件名都不含可解析的年份范围，回退到排序后的第一个文件。

    这样可避免把只覆盖部分年份（如 2015-2016）的旧/测试文件误选为主文件。
    """
    parsed = [(f, *_parse_year_range_from_filename(Path(f).name)) for f in files]

    req: tuple[int, int] | None = None
    if years:
        try:
            req = parse_years(years)
        except ValueError:
            req = None

    if req is not None:
        y0, y1 = req
        covering = [(f, a, b) for f, a, b in parsed if a is not None and a <= y0 and b >= y1]
        if covering:
            # 范围最紧（区间最短）优先，并列时取起始年较小的。
            covering.sort(key=lambda t: (t[2] - t[1], t[1]))
            return Path(covering[0][0])

    # 兜底：选文件名中年份范围最大的文件。
    ranged = [(f, a, b) for f, a, b in parsed if a is not None]
    if ranged:
        ranged.sort(key=lambda t: (-(t[2] - t[1]), t[1]))
        return Path(ranged[0][0])

    return Path(files[0])


def find_bcsd_file(
    data_dir: str | Path,
    model: str,
    region: str,
    scenario: str,
    var: str,
    years: str | None = None,
) -> Path:
    """使用 glob 查找指定变量文件。

    当匹配到多个文件时，按 years 参数选择“覆盖请求年份范围且范围最紧”的文件；
    若无法按年份筛选，则选择年份范围最大的文件，避免误选只覆盖部分年份的文件
    （例如残留的 2015-2016 测试文件导致结果仅覆盖 2015-2016）。
    """
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
        chosen = _select_best_bcsd_file(files, years)
        logger.warning(
            f"{var} 匹配到 {len(files)} 个文件，按年份范围({years or '未指定'})选择：{chosen}\n"
            + "  全部匹配文件：\n    " + "\n    ".join(files)
        )
        return chosen
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


def to_kw_m2(arr: np.ndarray, units: str | None) -> np.ndarray:
    """把辐照度转换为 kW m-2。"""
    x = arr.astype(np.float32, copy=False)
    units_l = (units or "").lower().replace(" ", "")
    if "kw" in units_l:
        out = x
    elif "w" in units_l and "m" in units_l:
        out = x / 1000.0
    else:
        # 没有单位时，按典型数量级判断：rsds 若最大值 > 10，通常是 W m-2。
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


def _erbs_diffuse_fraction(kt: np.ndarray) -> np.ndarray:
    """Erbs (1982) 经验模型：由晴空指数 kt 估计漫射比例。"""
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

    与原代码不同，这里显式加入经度：
      true solar time = UTC minutes + equation_of_time + 4 * longitude
    经度单位为 degree east。
    """
    doy = doy.astype(np.float64)
    hour_decimal_utc = hour_decimal_utc.astype(np.float64)
    lats = lats.astype(np.float64)
    lons = lons.astype(np.float64)

    b = 2.0 * np.pi * (doy - 1.0) / 365.0

    # Spencer/NOAA 近似太阳赤纬，单位 rad。
    decl = (
        0.006918
        - 0.399912 * np.cos(b)
        + 0.070257 * np.sin(b)
        - 0.006758 * np.cos(2 * b)
        + 0.000907 * np.sin(2 * b)
        - 0.002697 * np.cos(3 * b)
        + 0.001480 * np.sin(3 * b)
    )

    # equation of time, minutes。
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
    """
    保留原代码的 Erbs 辐射分解：GHI -> diffuse fraction。
    返回 diffuse fraction 和 cos(SZA)。
    """
    cos_sza = solar_cos_zenith(doy, hour_decimal_utc, lats, lons)
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
    """
    双轴跟踪下的组件平面辐照度，单位 kW m-2。

    保留当前代码的 direct + diffuse 分解思想：
      - DHI = GHI * fd
      - beam horizontal = GHI - DHI
      - 双轴跟踪时，直射分量按 DNI 入射，即 beam_horizontal / cos(SZA)
      - 散射分量仍采用原代码的 isotropic diffuse 近似。

    未额外加入地表反射项，以保持“当前辐射分解方法”。
    """
    dhi = ghi_kw * fd
    beam_horiz = np.maximum(ghi_kw - dhi, 0.0)

    with np.errstate(invalid="ignore", divide="ignore"):
        dni = np.where(cos_sza > MIN_COS_ZENITH, beam_horiz / cos_sza, 0.0)

    # 双轴跟踪：组件法线指向太阳，tilt beta = SZA，因此 cos(beta)=cos(SZA)。
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

    # CF = I_plane(kW m-2) * TEMcoef * SYScoef
    cf = plane_kw * temp_coef * SYS_COEF
    cf = np.clip(cf, 0.0, 1.0)
    return cf.astype(np.float32)


def build_time_index(
    time_da: xr.DataArray,
    years: str,
    months: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """根据 years/months 生成时间索引，并返回 doy/hour_decimal。"""
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



def _coord_values_close(a: np.ndarray, b: np.ndarray, name: str) -> bool:
    """检查空间坐标是否一致。"""
    if a.shape != b.shape:
        logger.error(f"{name} 坐标长度不一致：{a.shape} vs {b.shape}")
        return False
    if np.issubdtype(a.dtype, np.number) and np.issubdtype(b.dtype, np.number):
        ok = np.allclose(a, b, rtol=0.0, atol=1e-6, equal_nan=True)
    else:
        ok = np.array_equal(a, b)
    if not ok:
        logger.error(f"{name} 坐标值不一致，当前脚本不做空间插值。")
    return ok


def validate_same_spatial_grid(
    ds_ref: xr.Dataset,
    others: dict[str, xr.Dataset],
    lat_name: str,
    lon_name: str,
) -> None:
    """确保 rsds/tas/uas/vas 在同一空间网格上。"""
    ref_lat = ds_ref[lat_name].values
    ref_lon = ds_ref[lon_name].values
    for name, ds in others.items():
        if lat_name not in ds.coords and lat_name not in ds.dims:
            raise KeyError(f"{name} 文件中找不到纬度坐标 {lat_name}")
        if lon_name not in ds.coords and lon_name not in ds.dims:
            raise KeyError(f"{name} 文件中找不到经度坐标 {lon_name}")
        if not _coord_values_close(ref_lat, ds[lat_name].values, f"{name}.{lat_name}"):
            raise ValueError(f"{name} 与 rsds 的纬度网格不一致。")
        if not _coord_values_close(ref_lon, ds[lon_name].values, f"{name}.{lon_name}"):
            raise ValueError(f"{name} 与 rsds 的经度网格不一致。")


def _datetime64_to_ns(values: np.ndarray) -> np.ndarray:
    """将时间坐标转换为数值，便于 searchsorted。支持 datetime64 和 cftime。"""
    import cftime

    arr = np.asarray(values)
    if np.issubdtype(arr.dtype, np.datetime64):
        return arr.astype("datetime64[ns]").astype(np.int64)
    if arr.dtype == object and arr.size > 0 and isinstance(arr.flat[0], cftime.datetime):
        # cftime 类型（如 CANESM5 的 noleap 日历）：转为距第一个时间点的秒偏移量。
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

    这是方案4的核心：
      - target_times 来自 rsds.time，例如 01:30, 04:30, ...；
      - da 来自 tas/uas/vas，例如 03:00, 06:00, ...；
      - 对时间轴内部目标点做线性插值；
      - 对超出源时间范围的边界点，默认用最近的源时间值填充。

    注意：只沿时间维插值，不做空间插值。
    """
    src_times = da[time_name].values
    src_num = _datetime64_to_ns(src_times)
    tgt_num = _datetime64_to_ns(target_times)

    if src_num.ndim != 1:
        raise ValueError(f"{da.name or 'variable'} 的时间坐标不是一维。")
    if src_num.size < 1:
        raise ValueError(f"{da.name or 'variable'} 的时间坐标为空。")
    if np.any(np.diff(src_num) <= 0):
        raise ValueError(f"{da.name or 'variable'} 的时间坐标不是严格递增，无法插值。")

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
        # 内部仍线性插值，边界后续置 NaN。
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


def create_output_file(
    out_file: Path,
    rsds_file: Path,
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
    raw = xr.open_dataset(rsds_file, decode_times=False)
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
        "PV capacity factor computed from rsds using Erbs diffuse fraction, "
        "two-axis tracking geometry, paper temperature correction, and SYScoef=0.8056."
    )

    for k, v in attrs.items():
        try:
            setattr(nc, k, v)
        except Exception:
            pass

    return nc


def compute_region_solar_cf(
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
    """计算单个 model-region-scenario 的光伏容量因子。"""
    if not is_valid_region_name(region):
        raise ValueError(f"region 名称不合法或被过滤：{region}")

    rsds_file = find_bcsd_file(data_dir, model, region, scenario, "rsds", years)
    tas_file = find_bcsd_file(data_dir, model, region, scenario, "tas", years)
    uas_file = find_bcsd_file(data_dir, model, region, scenario, "uas", years)
    vas_file = find_bcsd_file(data_dir, model, region, scenario, "vas", years)

    logger.info(f"rsds: {rsds_file}")
    logger.info(f"tas : {tas_file}")
    logger.info(f"uas : {uas_file}")
    logger.info(f"vas : {vas_file}")

    out_file = (
        Path(output_dir) / model / region / f"solar_CF_{region}_{model}_{scenario}_{years}_{months_tag(months)}.nc"
    )
    if out_file.exists() and not overwrite:
        logger.info(f"✓ 已存在，跳过：{out_file}")
        return out_file

    ds_rsds = xr.open_dataset(rsds_file)
    ds_tas = xr.open_dataset(tas_file)
    ds_uas = xr.open_dataset(uas_file)
    ds_vas = xr.open_dataset(vas_file)

    rsds_var = get_var_name(ds_rsds, "rsds")
    tas_var = get_var_name(ds_tas, "tas")
    uas_var = get_var_name(ds_uas, "uas")
    vas_var = get_var_name(ds_vas, "vas")

    time_name = find_coord_name(ds_rsds, ["time", "valid_time"])
    lat_name = find_coord_name(ds_rsds, ["lat", "latitude"])
    lon_name = find_coord_name(ds_rsds, ["lon", "longitude"])

    # 方案4：以 rsds 的半点时间轴为主。
    # 不再对四个变量做 join="inner" 的严格时间交集，否则 01:30 与 03:00 不重合会导致时间步丢失。
    # tas/uas/vas 是瞬时变量，后续在每个时间块内线性插值到 rsds.time。
    validate_same_spatial_grid(
        ds_rsds,
        {"tas": ds_tas, "uas": ds_uas, "vas": ds_vas},
        lat_name=lat_name,
        lon_name=lon_name,
    )

    n_time_before = {
        "rsds": ds_rsds.sizes[time_name],
        "tas": ds_tas.sizes[time_name],
        "uas": ds_uas.sizes[time_name],
        "vas": ds_vas.sizes[time_name],
    }
    logger.info(f"输入时间步数量：{n_time_before}")
    logger.info("时间对齐策略：保留 rsds.time，将 tas/uas/vas 线性插值到 rsds.time（边界用最近值）。")

    time_idx, doy_all, hour_all = build_time_index(ds_rsds[time_name], years, months)
    lats = ds_rsds[lat_name].values.astype(np.float32)
    lons = ds_rsds[lon_name].values.astype(np.float32)

    n_time = len(time_idx)
    nlat = len(lats)
    nlon = len(lons)
    logger.info(f"匹配时间步: {n_time}；空间维度: lat={nlat}, lon={nlon}")

    land_mask = build_land_mask(lats, lons)  # (n_lat, n_lon) bool
    logger.info(f"陆地格点占比：{100.0 * land_mask.mean():.1f}%")

    attrs = {
        "model": model,
        "region": region,
        "scenario": scenario,
        "years": years,
        "months": months if months.strip() else "1-12",
        "ocean_handling": "ocean grid cells set to NaN via global-land-mask",
        "source_files": ", ".join(str(p) for p in [rsds_file, tas_file, uas_file, vas_file]),
        "time_alignment": "scheme4: keep rsds time axis; linearly interpolate instantaneous tas/uas/vas to rsds time; boundary uses nearest value",
        "radiation_decomposition": "Erbs diffuse fraction from rsds/GHI; same data requirement as original code",
        "tracking": "two-axis tracking; direct component uses DNI",
        "temperature_correction": "TEMcoef=1-gamma*(Tcell-25), Tcell=4.3+0.943*Tair+0.028*I-1.528*V",
        "system_coefficient": SYS_COEF,
    }

    nc = create_output_file(
        out_file=out_file,
        rsds_file=rsds_file,
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
    cf_var = nc.variables["solar_cf"]

    rsds_units = ds_rsds[rsds_var].attrs.get("units", "")
    tas_units = ds_tas[tas_var].attrs.get("units", "")

    logger.info(f"开始分块计算，chunk_time={chunk_time}")
    out_start = 0
    try:
        for start in tqdm(range(0, n_time, chunk_time), desc=f"{region}-{scenario}", unit="块"):
            end = min(start + chunk_time, n_time)
            idx_chunk = time_idx[start:end]

            # rsds 是主时间轴变量，直接按 rsds.time 读取。
            # tas/uas/vas 是瞬时变量，按方案4插值到当前 rsds 时间块。
            target_times = ds_rsds[time_name].isel({time_name: idx_chunk}).values

            rsds_raw = ds_rsds[rsds_var].isel({time_name: idx_chunk}).values
            tas_raw = interp_instantaneous_to_target_times_chunk(
                ds_tas[tas_var], time_name, target_times, fill_boundary="nearest"
            )
            uas_raw = interp_instantaneous_to_target_times_chunk(
                ds_uas[uas_var], time_name, target_times, fill_boundary="nearest"
            )
            vas_raw = interp_instantaneous_to_target_times_chunk(
                ds_vas[vas_var], time_name, target_times, fill_boundary="nearest"
            )

            rsds_kw = to_kw_m2(rsds_raw, rsds_units)
            tas_c = tas_to_celsius(tas_raw, tas_units)

            cf_chunk = compute_solar_cf_chunk(
                rsds_kw=rsds_kw,
                tas_c=tas_c,
                uas=uas_raw,
                vas=vas_raw,
                doy=doy_all[start:end],
                hour_decimal_utc=hour_all[start:end],
                lats=lats,
                lons=lons,
            )

            # 海洋格点恢复为 NaN（修复 nan_to_num 把海洋填 0 的问题）
            cf_chunk = np.where(land_mask[None, :, :], cf_chunk, np.nan).astype(np.float32)

            cf_var[out_start : out_start + (end - start), :, :] = cf_chunk
            out_start += end - start

            del rsds_raw, tas_raw, uas_raw, vas_raw, rsds_kw, tas_c, cf_chunk
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
        ds_rsds.close()
        ds_tas.close()
        ds_uas.close()
        ds_vas.close()

    logger.info(f"✓ 已保存：{out_file}")
    return out_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description="BCSD rsds/tas/uas/vas 光伏容量因子计算（Erbs + 双轴跟踪 + 论文温度修正）"
    )
    parser.add_argument("--data_dir", default="data/bcsd_outputs", help="数据目录，例如 data/bcsd_outputs")
    parser.add_argument("--model", required=True, help="模型名称，例如 MIROC-ES2H")
    parser.add_argument("--region", required=True, help="区域名称，例如 Austria；也可传 all 处理全部有效 region")
    parser.add_argument("--scenario", required=True, help="情景名称，例如 ssp126")
    parser.add_argument("--years", default="2015-2060", help="年段 YYYY-YYYY，包含两端")
    parser.add_argument("--months", default="", help="逗号分隔月份；默认全年 1-12")
    parser.add_argument("--output_dir", default="output/CFs_of_solar", help="输出目录")
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
        compute_region_solar_cf(
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
python S01E01_Simulate_Solar_CF_BCSD.py \
  --data_dir data/bcsd_outputs \
  --model MIROC-ES2H \
  --region Austria \
  --scenario ssp126 \
  --years 2015-2060

"""
