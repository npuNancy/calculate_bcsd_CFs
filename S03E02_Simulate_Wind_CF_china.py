"""
中国区域逐3小时风电容量因子计算（CMIP6 BCSD sfcWind）
=====================================================

本脚本用于基于 CMIP6 BCSD 降尺度数据计算中国区域陆上风电容量因子。

计算逻辑：
  1. 从 BCSD 数据中读取 sfcWind（10m 合成风速）；
  2. 使用幂律风速外推方法，将 10m 风速外推到轮毂高度：
       Vhub = V10 * (Hhub / 10)^alpha
     其中 Hhub = 100 m, alpha = 1/7；
  3. 使用 windpowerlib 读取 GE120/2500 风机功率曲线；
  4. 将轮毂高度风速代入功率曲线，计算风机输出功率；
  5. 容量因子 = power / rated_power，限制在 [0, 1]；
  6. 切出风速 25 m s-1；
  7. 按时间分块计算，逐块写入 NetCDF 文件。

风机参数：
  - TURBINE_TYPE = "GE120/2500"
  - HUB_HEIGHT = 100 m
  - REF_HEIGHT = 10 m
  - POWER_LAW_ALPHA = 1/7
  - CUT_OUT = 25 m s-1

输入数据路径格式：
  {data_dir}/{model}/{var}_3h_bcsd_on_0.1deg_china_{scenario}_{years}.nc

例如：
  data/cmip6_downscaling_3hr/MIROC-ES2H/
    sfcWind_3h_bcsd_on_0.1deg_china_ssp126_2015-2100.nc

输入变量要求：
  - sfcWind 文件中包含 sfcWind_bcsd 或唯一数据变量；
  - 数据维度为 (time, lat, lon)。

输出：
  {output_dir}/{model}/wind_CF_china_{model}_{scenario}_{years}_{months}.nc

输出变量：
  - wind_cf(time, lat, lon)，单位 1，范围 [0, 1]

注意：
  - 本脚本直接使用 sfcWind（10m 合成风速），而非 uas/vas 风速分量。
  - 输出时间分辨率与 sfcWind 输入一致（3 小时），不做时间插值。
  - 只计算陆上风电容量因子。
  - 各模型 sfcWind 起始时间戳不同（00:00/01:30/03:00），输出时间轴与 sfcWind 一致。
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


def find_china_bcsd_file(data_dir: str | Path, model: str, scenario: str, var: str) -> Path:
    """查找中国区域的 BCSD 降尺度文件。"""
    base = Path(data_dir) / model
    patterns = [
        str(base / f"{var}_3h_bcsd_on_0.1deg_china_{scenario}_*.nc"),
        str(base / f"{var}_*china_{scenario}_*.nc"),
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


def wind_to_ms(arr: np.ndarray, units: str | None) -> np.ndarray:
    """把风速转换为 m s-1。"""
    x = arr.astype(np.float32, copy=False)
    units_l = (units or "").lower().replace(" ", "")
    if "km/h" in units_l or "kmh-1" in units_l or "kmhr-1" in units_l:
        x = x / 3.6
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
    """幂律风速外推系数。"""
    return np.float32((hub_height / ref_height) ** alpha)


def _get_power_curve_arrays() -> tuple[np.ndarray, np.ndarray, float]:
    """从 windpowerlib 加载风机功率曲线，返回 (风速, 功率 kW, 额定功率 kW)。"""
    turbine = WindTurbine(turbine_type=TURBINE_TYPE, hub_height=HUB_HEIGHT)
    pc = turbine.power_curve
    ws = pc["wind_speed"].values.astype(np.float32)
    pw_kw = pc["value"].values.astype(np.float32) / 1000.0  # W -> kW
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
    elif ws[-1] > CUT_OUT:
        ws = np.append(ws[ws < CUT_OUT], np.float32(CUT_OUT))
        pw_kw = np.append(pw_kw[: len(ws) - 1], pw_kw[-1])

    return ws.astype(np.float32), pw_kw.astype(np.float32), rated_kw


def apply_power_curve(ws_hub: np.ndarray, ws_curve: np.ndarray, pw_curve: np.ndarray, rated_kw: float) -> np.ndarray:
    """向量化的功率曲线插值，返回容量因子 CF。"""
    ws_ext = np.concatenate([[0.0], ws_curve, [CUT_OUT + 0.01]]).astype(np.float32)
    pw_ext = np.concatenate([[0.0], pw_curve, [0.0]]).astype(np.float32)

    power = np.interp(ws_hub, ws_ext, pw_ext).astype(np.float32)
    cf = power / np.float32(rated_kw)
    cf = np.clip(cf, 0.0, 1.0)
    return cf.astype(np.float32)


def compute_wind_cf_chunk(
    sfcWind: np.ndarray,
    ws_curve: np.ndarray,
    pw_curve: np.ndarray,
    rated_kw: float,
    extrap_ratio: np.float32,
) -> np.ndarray:
    """计算一个时间块的陆上风电容量因子。"""
    ws10 = sfcWind.astype(np.float32)
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
    cf.long_name = "Onshore wind capacity factor (China)"
    cf.units = "1"
    cf.description = (
        "Onshore wind capacity factor for China computed from sfcWind (10m wind speed) "
        f"using power-law extrapolation Vhub=V10*(Hhub/10)^alpha with Hhub={HUB_HEIGHT:g} m, "
        f"alpha={POWER_LAW_ALPHA:.6f}, turbine={TURBINE_TYPE}, cut_out={CUT_OUT:g} m s-1."
    )

    for k, v in attrs.items():
        try:
            setattr(nc, k, v)
        except Exception:
            pass

    return nc


def compute_china_wind_cf(
    data_dir: str,
    model: str,
    scenario: str,
    years: str,
    months: str,
    output_dir: str,
    chunk_time: int = 512,
    overwrite: bool = False,
    compress_level: int = 4,
) -> Path:
    """计算中国区域陆上风电容量因子。"""
    sfcWind_file = find_china_bcsd_file(data_dir, model, scenario, "sfcWind")

    logger.info(f"sfcWind: {sfcWind_file}")

    out_file = (
        Path(output_dir) / model / f"wind_CF_china_{model}_{scenario}_{years}_{months_tag(months)}.nc"
    )
    if out_file.exists() and not overwrite:
        logger.info(f"✓ 已存在，跳过：{out_file}")
        return out_file

    ds_wind = xr.open_dataset(sfcWind_file)

    wind_var = get_var_name(ds_wind, "sfcWind")

    time_name = find_coord_name(ds_wind, ["time", "valid_time"])
    lat_name = find_coord_name(ds_wind, ["lat", "latitude"])
    lon_name = find_coord_name(ds_wind, ["lon", "longitude"])

    time_idx = build_time_index(ds_wind[time_name], years, months)

    lats = ds_wind[lat_name].values.astype(np.float32)
    lons = ds_wind[lon_name].values.astype(np.float32)

    n_time = len(time_idx)
    nlat = len(lats)
    nlon = len(lons)
    logger.info(f"匹配时间步: {n_time}；空间维度: lat={nlat}, lon={nlon}")

    ws_curve, pw_curve, rated_kw = _get_power_curve_arrays()
    extrap_ratio = power_law_ratio()
    logger.info(
        f"风速外推: Vhub=V10*({HUB_HEIGHT:.1f}/{REF_HEIGHT:.1f})^{POWER_LAW_ALPHA:.6f}，ratio={float(extrap_ratio):.4f}"
    )
    logger.info(f"风机: {TURBINE_TYPE}；rated={rated_kw:.1f} kW；cut-out={CUT_OUT:.1f} m/s")

    attrs = {
        "model": model,
        "region": "china",
        "scenario": scenario,
        "years": years,
        "months": months if months.strip() else "1-12",
        "source_files": str(sfcWind_file),
        "technology": "onshore wind only",
        "turbine_type": TURBINE_TYPE,
        "hub_height_m": HUB_HEIGHT,
        "reference_height_m": REF_HEIGHT,
        "power_law_alpha": POWER_LAW_ALPHA,
        "power_law_ratio": float(extrap_ratio),
        "cut_out_speed_ms": CUT_OUT,
        "farm_efficiency_multiplier": "not applied",
    }

    nc = create_output_file(
        out_file=out_file,
        template_file=sfcWind_file,
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

    wind_units = ds_wind[wind_var].attrs.get("units", "")

    total_sum = 0.0
    total_count = 0
    positive_sum = 0.0
    positive_count = 0
    max_cf = 0.0

    logger.info(f"开始分块计算，chunk_time={chunk_time}")
    out_start = 0
    try:
        for start in tqdm(range(0, n_time, chunk_time), desc=f"{model}-{scenario}", unit="块"):
            end = min(start + chunk_time, n_time)
            idx_chunk = time_idx[start:end]

            sfcWind_raw = ds_wind[wind_var].isel({time_name: idx_chunk}).values
            sfcWind_ms = wind_to_ms(sfcWind_raw, wind_units)

            cf_chunk = compute_wind_cf_chunk(
                sfcWind=sfcWind_ms,
                ws_curve=ws_curve,
                pw_curve=pw_curve,
                rated_kw=rated_kw,
                extrap_ratio=extrap_ratio,
            )

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

            del sfcWind_raw, sfcWind_ms, cf_chunk
            gc.collect()
    finally:
        nc.close()
        ds_wind.close()

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
        description="CMIP6 BCSD 中国区域陆上风电容量因子计算（sfcWind + 幂律外推 + windpowerlib 功率曲线）"
    )
    parser.add_argument("--data_dir", default="data/cmip6_downscaling_3hr", help="数据目录")
    parser.add_argument("--model", required=True, help="模型名称，如 BCC-CSM2-MR, MIROC-ES2H 等")
    parser.add_argument("--scenario", required=True, help="情景名称：historical, ssp126, ssp245, ssp370, ssp585")
    parser.add_argument("--years", default="2015-2100", help="年段 YYYY-YYYY，包含两端")
    parser.add_argument("--months", default="", help="逗号分隔月份；默认全年 1-12")
    parser.add_argument("--output_dir", default="output/CFs_of_wind_china", help="输出目录")
    parser.add_argument("--chunk_time", type=int, default=512, help="每块处理的时间步数")
    parser.add_argument("--compress_level", type=int, default=4, help="NetCDF 压缩级别 0-9")
    parser.add_argument("--overwrite", action="store_true", help="若输出文件已存在，则覆盖重算")
    args = parser.parse_args()

    compute_china_wind_cf(
        data_dir=args.data_dir,
        model=args.model,
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

python S03E02_Simulate_Wind_CF_china.py \
  --data_dir data/cmip6_downscaling_3hr \
  --model BCC-CSM2-MR \
  --scenario ssp126 \
  --years 2015-2100

python S03E02_Simulate_Wind_CF_china.py \
  --data_dir data/cmip6_downscaling_3hr \
  --model MIROC-ES2H \
  --scenario historical \
  --years 1979-2014

"""
