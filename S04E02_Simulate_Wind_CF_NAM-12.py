"""
CMIP6-CORDEX NAM-12 逐小时陆上风电容量因子计算（uas + vas）
============================================================

本脚本用于基于 CMIP6-CORDEX NAM-12 逐小时数据计算陆上风电容量因子。

计算逻辑：
  1. 从 CORDEX 数据中读取 uas 和 vas 两个近地面风速分量；
  2. 将 uas 和 vas 合成为 10 m 风速：V10 = sqrt(uas^2 + vas^2)
  3. 使用幂律风速外推方法，将 10 m 风速外推到轮毂高度：
       Vhub = V10 * (Hhub / 10)^alpha
     其中 Hhub = 100 m, alpha = 1/7
  4. 使用 windpowerlib 读取 GE120/2500 风机功率曲线；
  5. 将轮毂高度风速代入功率曲线，计算风机输出功率；
  6. 容量因子 = power / rated_power，限制在 [0, 1]；
  7. 按年循环处理，每年内按时间分块计算，逐块写入 NetCDF。

风机参数：
  - TURBINE_TYPE = "GE120/2500"
  - HUB_HEIGHT = 100 m
  - REF_HEIGHT = 10 m
  - POWER_LAW_ALPHA = 1/7
  - CUT_OUT = 25 m s-1

输入数据路径格式：
  {data_dir}/{gcm_model}/{realization}/{rcm_model}/{scenario}/{var}/
  例如：data/CORDEX-CMIP6/NAM-12/1hr/MPI-ESM1-2-LR/r1i1p1f1/CRCM5/ssp126/uas/

输入变量要求：
  - uas 文件中包含 uas 变量；
  - vas 文件中包含 vas 变量；
  - 数据维度为 (time, rlat, rlon)，旋转极网格；
  - 每年一个文件。

输出：
  {output_dir}/{gcm_model}/{realization}/
    wind_CF_NAM-12_{gcm_model}_{realization}_{rcm_model}_{scenario}_{years}_{months_tag}.nc

输出变量：
  - wind_cf(time, rlat, rlon)
  - 单位：1
  - 数值范围：[0, 1]

注意：
  - uas/vas 时间均在整点（:00），不需要时间对齐。
  - 旋转极网格：rlat/rlon 为主维度，lat/lon 为二维辅助坐标。
  - 仅计算陆上风电容量因子。

示例：
  python S04E02_Simulate_Wind_CF_NAM-12.py \\
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
from windpowerlib import WindTurbine

# ── 日志配置 ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── 风机与风速外推参数 ───────────────────────────────────────────────────────
TURBINE_TYPE = "GE120/2500"
HUB_HEIGHT = 100.0  # m
REF_HEIGHT = 10.0  # m
POWER_LAW_ALPHA = 1.0 / 7.0
CUT_OUT = 25.0  # m s-1


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
    """获取变量名。优先使用 preferred，否则在唯一 data_var 时使用该变量。"""
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
    """获取变量，压缩标量维度（如 height），转置为 (time, rlat, rlon)。"""
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


def cordex_var_dir(
    data_dir: str | Path,
    gcm_model: str,
    realization: str,
    rcm_model: str,
    scenario: str,
    var: str,
) -> Path:
    """返回 CORDEX 某变量的目录路径。"""
    return Path(data_dir) / gcm_model / realization / rcm_model / scenario / var


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
    var_dir = cordex_var_dir(data_dir, gcm_model, realization, rcm_model, scenario, var)
    if not var_dir.is_dir():
        raise FileNotFoundError(f"变量目录不存在：{var_dir}")
    pattern = str(var_dir / f"{var}_NAM-12_{gcm_model}_*.nc")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"找不到 {var} 文件。已尝试：{pattern}")

    y0, y1 = years_tuple
    result: dict[int, Path] = {}
    for f in files:
        # 从文件名解析年份，例如 ..._1hr_202001010000-202012312300.nc
        m = re.search(r"_(\d{4})\d{8}-\d{12}\.nc$", Path(f).name)
        if m:
            year = int(m.group(1))
            if y0 <= year <= y1:
                result[year] = Path(f)
    if not result:
        raise FileNotFoundError(f"在 {var_dir} 中找不到 {y0}-{y1} 范围内的 {var} 文件")
    return result


# ─────────────────────────────────────────────────────────────────────────────
# 2. 单位转换与风速外推
# ─────────────────────────────────────────────────────────────────────────────


def wind_to_ms(arr: np.ndarray, units: str | None) -> np.ndarray:
    """把风速转换为 m s-1。"""
    x = arr.astype(np.float32, copy=False)
    units_l = (units or "").lower().replace(" ", "")
    if "km/h" in units_l or "kmh-1" in units_l or "kmhr-1" in units_l:
        x = x / 3.6
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    return x.astype(np.float32)


def power_law_ratio(
    hub_height: float = HUB_HEIGHT,
    ref_height: float = REF_HEIGHT,
    alpha: float = POWER_LAW_ALPHA,
) -> np.float32:
    """幂律风速外推系数。"""
    return np.float32((hub_height / ref_height) ** alpha)


def _get_power_curve_arrays() -> tuple[np.ndarray, np.ndarray, float]:
    """从 windpowerlib 数据库加载风机功率曲线。"""
    turbine = WindTurbine(turbine_type=TURBINE_TYPE, hub_height=HUB_HEIGHT)
    pc = turbine.power_curve
    ws = pc["wind_speed"].values.astype(np.float32)
    pw_kw = pc["value"].values.astype(np.float32) / 1000.0
    rated_kw = float(turbine.nominal_power / 1000.0)

    order = np.argsort(ws)
    ws = ws[order]
    pw_kw = pw_kw[order]
    keep = ws <= CUT_OUT
    ws = ws[keep]
    pw_kw = pw_kw[keep]

    if ws.size == 0:
        raise RuntimeError(f"风机 {TURBINE_TYPE} 的功率曲线在 CUT_OUT={CUT_OUT} 以内为空。")

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

    return ws.astype(np.float32), pw_kw.astype(np.float32), rated_kw


def apply_power_curve(ws_hub: np.ndarray, ws_curve: np.ndarray, pw_curve: np.ndarray, rated_kw: float) -> np.ndarray:
    """向量化的功率曲线插值，返回容量因子。"""
    ws_ext = np.concatenate([[0.0], ws_curve, [CUT_OUT + 0.01]]).astype(np.float32)
    pw_ext = np.concatenate([[0.0], pw_curve, [0.0]]).astype(np.float32)
    power = np.interp(ws_hub, ws_ext, pw_ext).astype(np.float32)
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


# ─────────────────────────────────────────────────────────────────────────────
# 3. 输出文件创建（旋转极网格）
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

    # 用 decode_times=False 复制原始 time 数值和属性
    raw = xr.open_dataset(template_file, decode_times=False)
    raw_time = raw["time"].values
    time_attrs = dict(raw["time"].attrs)
    raw.close()
    ds_tpl.close()

    nrlat = len(rlat_vals)
    nrlon = len(rlon_vals)

    nc = Dataset(out_file, "w", format="NETCDF4")

    # 全局属性
    if attrs:
        for k, v in attrs.items():
            try:
                setattr(nc, k, v)
            except Exception:
                pass

    # 维度
    nc.createDimension("time", total_time_steps)
    nc.createDimension(rlat_name, nrlat)
    nc.createDimension(rlon_name, nrlon)

    # 时间变量（先创建占位，后续逐年写入）
    tvar = nc.createVariable("time", raw_time.dtype, ("time",))
    for k, v in time_attrs.items():
        try:
            setattr(tvar, k, v)
        except Exception:
            pass

    # rlat / rlon
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

    # lat / lon 二维辅助坐标
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

    # crs 网格映射
    crs_var = nc.createVariable("crs", "i1")
    for k, v in crs_attrs.items():
        try:
            setattr(crs_var, k, v)
        except Exception:
            pass

    # wind_cf
    chunksizes = (min(chunk_time, total_time_steps), nrlat, nrlon)
    cf = nc.createVariable(
        "wind_cf",
        "f4",
        ("time", rlat_name, rlon_name),
        zlib=True,
        complevel=compress_level,
        chunksizes=chunksizes,
        fill_value=np.float32(np.nan),
    )
    cf.long_name = "Onshore wind capacity factor"
    cf.units = "1"
    cf.grid_mapping = "crs"
    cf.description = (
        "Onshore wind capacity factor computed from CORDEX uas/vas using power-law extrapolation "
        f"Vhub=V10*(Hhub/10)^alpha with Hhub={HUB_HEIGHT:g} m, alpha={POWER_LAW_ALPHA:.6f}, "
        f"turbine={TURBINE_TYPE}, cut_out={CUT_OUT:g} m s-1."
    )

    return nc


# ─────────────────────────────────────────────────────────────────────────────
# 4. 合并年度文件与单年计算
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


def _compute_wind_cf_year(
    yr: int,
    uas_file: Path,
    vas_file: Path,
    yr_out_file: Path,
    ws_curve: np.ndarray,
    pw_curve: np.ndarray,
    rated_kw: float,
    extrap_ratio: np.float32,
    chunk_time: int = 24,
    compress_level: int = 4,
) -> Path:
    """计算单年风电 CF 并保存到独立 NC 文件。"""
    ds_uas = xr.open_dataset(uas_file)
    ds_vas = xr.open_dataset(vas_file)

    time_name = find_coord_name(ds_uas, ["time"])
    rlat_name = find_coord_name(ds_uas, ["rlat"])
    rlon_name = find_coord_name(ds_uas, ["rlon"])

    # 让 vas 的坐标名与 uas 一致
    vas_time = find_coord_name(ds_vas, ["time"])
    vas_rlat = find_coord_name(ds_vas, ["rlat"])
    vas_rlon = find_coord_name(ds_vas, ["rlon"])
    rename_map: dict[str, str] = {}
    if vas_time != time_name:
        rename_map[vas_time] = time_name
    if vas_rlat != rlat_name:
        rename_map[vas_rlat] = rlat_name
    if vas_rlon != rlon_name:
        rename_map[vas_rlon] = rlon_name
    if rename_map:
        ds_vas = ds_vas.rename(rename_map)

    uas_da = prepare_dataarray(ds_uas, "uas", time_name, rlat_name, rlon_name)
    vas_da = prepare_dataarray(ds_vas, "vas", time_name, rlat_name, rlon_name)

    # 对齐
    ds_uas_aligned, ds_vas_aligned = xr.align(ds_uas, ds_vas, join="inner", copy=False)
    n_time = ds_uas_aligned.sizes[time_name]
    if ds_uas_aligned.sizes[time_name] != ds_uas.sizes[time_name]:
        uas_da = prepare_dataarray(ds_uas_aligned, "uas", time_name, rlat_name, rlon_name)
        vas_da = prepare_dataarray(ds_vas_aligned, "vas", time_name, rlat_name, rlon_name)
        logger.warning(f"{yr} 对齐后时间步变为 {n_time}")

    uas_units = ds_uas[uas_da.name].attrs.get("units", "")
    vas_units = ds_vas[vas_da.name].attrs.get("units", "")

    logger.info(f"  处理 {yr}（{n_time} 步）...")

    nc = create_output_file_rotated(
        out_file=yr_out_file,
        template_file=uas_file,
        total_time_steps=n_time,
        compress_level=compress_level,
        chunk_time=chunk_time,
    )
    cf_var = nc.variables["wind_cf"]
    tvar = nc.variables["time"]

    try:
        # 写入 time 值
        raw = xr.open_dataset(uas_file, decode_times=False)
        tvar[:n_time] = raw["time"].values[:n_time]
        raw.close()

        for start in tqdm(range(0, n_time, chunk_time), desc=f"wind-{yr}", unit="块"):
            end = min(start + chunk_time, n_time)

            uas_raw = uas_da.isel({time_name: slice(start, end)}).values
            vas_raw = vas_da.isel({time_name: slice(start, end)}).values

            uas_ms = wind_to_ms(uas_raw, uas_units)
            vas_ms = wind_to_ms(vas_raw, vas_units)

            cf_chunk = compute_wind_cf_chunk(uas_ms, vas_ms, ws_curve, pw_curve, rated_kw, extrap_ratio)

            cf_var[start:end, :, :] = cf_chunk

            del uas_raw, vas_raw, uas_ms, vas_ms, cf_chunk
            gc.collect()
    except BaseException:
        if nc.isopen():
            nc.close()
        if yr_out_file.exists():
            logger.error(f"{yr} 计算出错，删除：{yr_out_file}")
            yr_out_file.unlink()
        ds_uas.close()
        ds_vas.close()
        raise

    nc.close()
    ds_uas.close()
    ds_vas.close()
    gc.collect()
    return yr_out_file


# ─────────────────────────────────────────────────────────────────────────────
# 5. 主计算函数
# ─────────────────────────────────────────────────────────────────────────────


def compute_nam12_wind_cf(
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
    """计算 CORDEX NAM-12 陆上风电容量因子（逐年独立计算 + 合并）。"""
    y0, y1 = parse_years(years)

    uas_files = cordex_var_files(data_dir, gcm_model, realization, rcm_model, scenario, "uas", (y0, y1))
    vas_files = cordex_var_files(data_dir, gcm_model, realization, rcm_model, scenario, "vas", (y0, y1))

    common_years = sorted(set(uas_files) & set(vas_files))
    if not common_years:
        raise FileNotFoundError(f"uas/vas 没有共同年份（{y0}-{y1}）")
    missing_uas = set(vas_files) - set(uas_files)
    missing_vas = set(uas_files) - set(vas_files)
    if missing_uas:
        logger.warning(f"uas 缺少年份：{sorted(missing_uas)}")
    if missing_vas:
        logger.warning(f"vas 缺少年份：{sorted(missing_vas)}")

    logger.info(f"待处理年份：{common_years[0]}-{common_years[-1]}，共 {len(common_years)} 年")

    out_file = (
        Path(output_dir)
        / gcm_model
        / realization
        / scenario
        / f"wind_CF_NAM-12_{gcm_model}_{realization}_{rcm_model}_{scenario}_{years}_{months_tag(months)}.nc"
    )
    if out_file.exists() and not overwrite:
        logger.info(f"已存在，跳过：{out_file}")
        return out_file

    # 加载功率曲线和外推系数
    ws_curve, pw_curve, rated_kw = _get_power_curve_arrays()
    extrap_ratio = power_law_ratio()
    logger.info(f"风速外推 ratio={float(extrap_ratio):.4f}，风机 {TURBINE_TYPE}，rated={rated_kw:.1f} kW")

    # 年度文件目录
    yearly_dir = out_file.parent / "yearly"
    yearly_dir.mkdir(parents=True, exist_ok=True)

    # 逐年计算
    yearly_files: list[Path] = []
    for yr in common_years:
        yr_file = (
            yearly_dir / f"wind_CF_NAM-12_{gcm_model}_{realization}_{rcm_model}_{scenario}_{yr}_{months_tag(months)}.nc"
        )
        if yr_file.exists() and not overwrite:
            logger.info(f"  {yr} 已存在，跳过")
            yearly_files.append(yr_file)
            continue

        logger.info(f"处理 {yr}...")
        try:
            _compute_wind_cf_year(
                yr=yr,
                uas_file=uas_files[yr],
                vas_file=vas_files[yr],
                yr_out_file=yr_file,
                ws_curve=ws_curve,
                pw_curve=pw_curve,
                rated_kw=rated_kw,
                extrap_ratio=extrap_ratio,
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
        _merge_yearly_nc_files(yearly_files, out_file, "wind_cf", compress_level, chunk_time)
        logger.info(f"已合并保存：{out_file}")
        logger.info(f"可删除年度文件目录：{yearly_dir}")

    else:
        logger.info(f"逐年文件已保存到：{yearly_dir}")
        logger.info(f"未执行合并，单年文件保留。")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CMIP6-CORDEX NAM-12 陆上风电容量因子计算（幂律外推 + windpowerlib 功率曲线）"
    )
    parser.add_argument("--data_dir", default="data/CORDEX-CMIP6/NAM-12/1hr", help="数据根目录")
    parser.add_argument("--gcm_model", required=True, help="GCM 模型，例如 MPI-ESM1-2-LR")
    parser.add_argument("--realization", required=True, help="实现编号，例如 r1i1p1f1")
    parser.add_argument("--rcm_model", default="CRCM5", help="RCM 模型，默认 CRCM5")
    parser.add_argument("--scenario", required=True, help="情景，例如 ssp126")
    parser.add_argument("--years", default="2020-2060", help="年段 YYYY-YYYY")
    parser.add_argument("--months", default="", help="逗号分隔月份；默认全年")
    parser.add_argument("--output_dir", default="output/CFs_of_wind_NAM-12", help="输出目录")
    parser.add_argument("--chunk_time", type=int, default=24, help="每块时间步数")
    parser.add_argument("--compress_level", type=int, default=4, help="NetCDF 压缩级别 0-9")
    parser.add_argument("--overwrite", action="store_true", help="覆盖已有输出")
    parser.add_argument("--merge", action="store_true", help="执行合并, 默认不合并, 逐年输出")
    args = parser.parse_args()

    compute_nam12_wind_cf(
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

python S04E02_Simulate_Wind_CF_NAM-12.py \\
  --gcm_model MPI-ESM1-2-LR --realization r1i1p1f1 \\
  --scenario ssp126 --years 2020-2060

"""
