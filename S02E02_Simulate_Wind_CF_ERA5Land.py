"""
全球逐时陆上风电容量因子计算（ERA5-Land u10/v10 + windpowerlib）
=================================================================

本脚本用于从 ERA5-Land 月度逐小时数据计算全球陆上风电容量因子。
计算方法与 S01E02_Simulate_Wind_CF_BCSD.py 保持一致：
  1. 从 ERA5-Land 读取 u10 和 v10；
  2. 合成 10 m 风速 V10 = sqrt(u10^2 + v10^2)；
  3. 使用幂律方法外推至 100 m 轮毂高度：Vhub = V10 * (100 / 10)^(1/7)；
  4. 使用 windpowerlib 中 GE120/2500 的功率曲线；
  5. 切出风速设为 25 m s-1；
  6. 容量因子定义为 power / rated_power，不额外乘风场效率系数；
  7. 按时间分块计算，并逐块写入 NetCDF，避免一次性保存全月数组。

输入数据路径格式
----------------
  {data_dir}/ERA5_land/global/u10/u10_YYYY_MM.nc
  {data_dir}/ERA5_land/global/v10/v10_YYYY_MM.nc

输出
----
  {output_dir}/wind_cf_YYYY_MM.nc

输出变量
--------
  wind_cf(time, lat, lon), 单位 1, 数值范围 [0, 1]

示例
----
python S02E02_Simulate_Wind_CF_ERA5Land.py \
  --data_dir data \
  --years 2020 \
  --months 1 \
  --output_dir output/CFs_of_wind \
  --chunk_time 72
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
from windpowerlib import WindTurbine
from global_land_mask import globe

# ── 日志配置 ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── 与 BCSD 脚本一致的风机与风速外推参数 ───────────────────────────────────
TURBINE_TYPE = "GE120/2500"  # GE 2.5-120；windpowerlib 数据库中与论文 GE 2.5 MW 相近的型号
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


def iter_year_months(years: str, months: str) -> list[tuple[int, int]]:
    """生成待处理的 (year, month) 列表。"""
    y0, y1 = parse_years(years)
    month_list = parse_months(months)
    return [(y, m) for y in range(y0, y1 + 1) for m in month_list]


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


def rename_coords_to_match(
    ds: xr.Dataset,
    target_time: str,
    target_lat: str,
    target_lon: str,
) -> xr.Dataset:
    """将另一个 ERA5-Land 数据集的坐标名改成与主数据集一致。"""
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


def validate_same_grid(
    ref: xr.Dataset, other: xr.Dataset, name: str, time_name: str, lat_name: str, lon_name: str
) -> None:
    """检查两个数据集是否在同一时间和空间网格上。"""
    if not np.array_equal(ref[time_name].values, other[time_name].values):
        raise ValueError(f"{name} 与 u10 的时间坐标不一致。ERA5-Land 月文件应逐小时完全一致。")
    if not np.allclose(ref[lat_name].values, other[lat_name].values, rtol=0.0, atol=1e-6, equal_nan=True):
        raise ValueError(f"{name} 与 u10 的纬度网格不一致。")
    if not np.allclose(ref[lon_name].values, other[lon_name].values, rtol=0.0, atol=1e-6, equal_nan=True):
        raise ValueError(f"{name} 与 u10 的经度网格不一致。")


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


def wind_to_ms(arr: np.ndarray, units: str | None) -> np.ndarray:
    """把风速分量转换为 m s-1。ERA5-Land u10/v10 通常已经是 m s-1。"""
    x = arr.astype(np.float32, copy=False)
    units_l = (units or "").lower().replace(" ", "")
    if "km/h" in units_l or "kmh-1" in units_l or "kmhr-1" in units_l:
        x = x / 3.6
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    return x.astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# 1. 风机功率曲线与容量因子计算
# ─────────────────────────────────────────────────────────────────────────────


def power_law_ratio(
    hub_height: float = HUB_HEIGHT,
    ref_height: float = REF_HEIGHT,
    alpha: float = POWER_LAW_ALPHA,
) -> np.float32:
    """幂律风速外推系数。"""
    return np.float32((hub_height / ref_height) ** alpha)


def _get_power_curve_arrays() -> tuple[np.ndarray, np.ndarray, float]:
    """从 windpowerlib 数据库加载风机功率曲线，返回 (风速, 功率 kW, 额定功率 kW)。"""
    turbine = WindTurbine(turbine_type=TURBINE_TYPE, hub_height=HUB_HEIGHT)
    pc = turbine.power_curve
    ws = pc["wind_speed"].values.astype(np.float32)
    pw_kw = pc["value"].values.astype(np.float32) / 1000.0  # W -> kW
    rated_kw = float(turbine.nominal_power / 1000.0)  # W -> kW

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
    """
    向量化功率曲线插值，返回容量因子。
    不额外乘 FARM_EFFICIENCY；容量因子仅由轮毂高度风速和单机功率曲线决定。
    """
    ws_ext = np.concatenate([[0.0], ws_curve, [CUT_OUT + 0.01]]).astype(np.float32)
    pw_ext = np.concatenate([[0.0], pw_curve, [0.0]]).astype(np.float32)
    power = np.interp(ws_hub, ws_ext, pw_ext).astype(np.float32)
    cf = power / np.float32(rated_kw)
    cf = np.clip(cf, 0.0, 1.0)
    return cf.astype(np.float32)


def compute_wind_cf_chunk(
    u10: np.ndarray,
    v10: np.ndarray,
    ws_curve: np.ndarray,
    pw_curve: np.ndarray,
    rated_kw: float,
    extrap_ratio: np.float32,
) -> np.ndarray:
    """计算一个时间块的陆上风电容量因子。"""
    ws10 = np.sqrt(u10.astype(np.float32) ** 2 + v10.astype(np.float32) ** 2)
    ws_hub = ws10 * extrap_ratio
    cf = apply_power_curve(ws_hub, ws_curve, pw_curve, rated_kw)
    return cf.astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# 2. NetCDF 输出
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
        "Onshore wind CF computed from ERA5-Land u10/v10 using power-law extrapolation "
        f"Vhub=V10*({HUB_HEIGHT:g}/{REF_HEIGHT:g})^{POWER_LAW_ALPHA:.6f}, "
        f"turbine={TURBINE_TYPE}, cut_out={CUT_OUT:g} m s-1. "
        "No farm-efficiency multiplier is applied."
    )

    for k, v in attrs.items():
        try:
            setattr(nc, k, v)
        except Exception:
            pass

    return nc


# ─────────────────────────────────────────────────────────────────────────────
# 3. 主计算函数
# ─────────────────────────────────────────────────────────────────────────────


def compute_month_wind_cf(
    data_dir: str,
    year: int,
    month: int,
    output_dir: str,
    chunk_time: int = 72,
    overwrite: bool = False,
    compress_level: int = 4,
    *,
    ws_curve: np.ndarray,
    pw_curve: np.ndarray,
    rated_kw: float,
    extrap_ratio: np.float32,
) -> Path:
    """计算单个月份的 ERA5-Land 陆上风电容量因子。"""
    tag = f"{year:04d}_{month:02d}"
    u10_file = era5land_file(data_dir, "u10", year, month)
    v10_file = era5land_file(data_dir, "v10", year, month)

    for p in [u10_file, v10_file]:
        if not p.exists():
            raise FileNotFoundError(f"找不到输入文件：{p}")

    out_file = Path(output_dir) / f"wind_cf_{tag}.nc"
    tmp_file = out_file.with_name(f".{out_file.name}.tmp")
    if out_file.exists() and not overwrite:
        logger.info(f"✓ 已存在，跳过：{out_file}")
        return out_file
    if tmp_file.exists():
        tmp_file.unlink()

    logger.info(f"u10: {u10_file}")
    logger.info(f"v10: {v10_file}")

    ds_u10 = xr.open_dataset(u10_file)
    ds_v10 = xr.open_dataset(v10_file)

    time_name = find_coord_name(ds_u10, ["valid_time", "time"])
    lat_name = find_coord_name(ds_u10, ["latitude", "lat"])
    lon_name = find_coord_name(ds_u10, ["longitude", "lon"])

    ds_v10 = rename_coords_to_match(ds_v10, time_name, lat_name, lon_name)
    validate_same_grid(ds_u10, ds_v10, "v10", time_name, lat_name, lon_name)

    u10_da = prepare_dataarray(ds_u10, "u10", time_name, lat_name, lon_name)
    v10_da = prepare_dataarray(ds_v10, "v10", time_name, lat_name, lon_name)

    n_time = u10_da.sizes[time_name]
    lats = ds_u10[lat_name].values.astype(np.float32)
    lons = ds_u10[lon_name].values.astype(np.float32)
    logger.info(f"维度：time={n_time}, lat={len(lats)}, lon={len(lons)}")
    logger.info(
        f"风速外推：Vhub=V10*({HUB_HEIGHT:.1f}/{REF_HEIGHT:.1f})^{POWER_LAW_ALPHA:.6f}，ratio={float(extrap_ratio):.4f}"
    )
    logger.info(f"风机：{TURBINE_TYPE}；rated={rated_kw:.1f} kW；cut-out={CUT_OUT:.1f} m/s")
    logger.info("输出坐标统一为：time/lat/lon；变量名：wind_cf")

    land_mask = build_land_mask(lats, lons)  # (n_lat, n_lon) bool
    logger.info(f"陆地格点占比：{100.0 * land_mask.mean():.1f}%")

    attrs = {
        "source": "ERA5-Land",
        "year": year,
        "month": month,
        "source_files": ", ".join(str(p) for p in [u10_file, v10_file]),
        "technology": "onshore wind only",
        "turbine_type": TURBINE_TYPE,
        "hub_height_m": HUB_HEIGHT,
        "reference_height_m": REF_HEIGHT,
        "power_law_alpha": POWER_LAW_ALPHA,
        "power_law_ratio": float(extrap_ratio),
        "cut_out_speed_ms": CUT_OUT,
        "farm_efficiency_multiplier": "not applied",
        "notes": "This script uses ERA5-Land u10/v10 and computes onshore wind CF only.",
        "ocean_handling": "ocean grid cells set to NaN via global-land-mask",
    }

    nc = create_output_file(
        out_file=tmp_file,
        template_file=u10_file,
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
    cf_var = nc.variables["wind_cf"]

    u10_units = u10_da.attrs.get("units", "")
    v10_units = v10_da.attrs.get("units", "")

    total_sum = 0.0
    total_count = 0
    positive_sum = 0.0
    positive_count = 0
    max_cf = 0.0

    try:
        logger.info(f"开始分块计算，chunk_time={chunk_time}")
        for start in tqdm(range(0, n_time, chunk_time), desc=f"wind-{tag}", unit="块"):
            end = min(start + chunk_time, n_time)

            u10_raw = u10_da.isel({time_name: slice(start, end)}).values
            v10_raw = v10_da.isel({time_name: slice(start, end)}).values

            u10_ms = wind_to_ms(u10_raw, u10_units)
            v10_ms = wind_to_ms(v10_raw, v10_units)
            del u10_raw, v10_raw

            cf_chunk = compute_wind_cf_chunk(
                u10=u10_ms,
                v10=v10_ms,
                ws_curve=ws_curve,
                pw_curve=pw_curve,
                rated_kw=rated_kw,
                extrap_ratio=extrap_ratio,
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

            del u10_ms, v10_ms, cf_chunk
            gc.collect()
    finally:
        nc.close()
        ds_u10.close()
        ds_v10.close()

    os.replace(tmp_file, out_file)
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
        description="ERA5-Land u10/v10 陆上风电容量因子计算（幂律外推 + windpowerlib 功率曲线）"
    )
    parser.add_argument("--data_dir", default="data", help="数据根目录，例如 data")
    parser.add_argument("--years", default="2015-2024", help="年段 YYYY-YYYY 或 YYYY，包含两端")
    parser.add_argument("--months", default="", help="逗号分隔月份；默认全年 1-12")
    parser.add_argument("--output_dir", default="output/CFs_of_wind_ERA5Land", help="输出目录")
    parser.add_argument("--chunk_time", type=int, default=72, help="每块处理的小时数；全球 0.1° 数据建议 24-168")
    parser.add_argument("--compress_level", type=int, default=4, help="NetCDF 压缩级别 0-9")
    parser.add_argument("--overwrite", action="store_true", help="若输出文件已存在，则覆盖重算")
    args = parser.parse_args()

    ws_curve, pw_curve, rated_kw = _get_power_curve_arrays()
    extrap_ratio = power_law_ratio()

    year_months = iter_year_months(args.years, args.months)
    logger.info(f"共 {len(year_months)} 个年月待处理：{year_months[0]} -> {year_months[-1]}")
    for i, (year, month) in enumerate(year_months, 1):
        logger.info(f"\n{'=' * 80}")
        logger.info(f"[{i}/{len(year_months)}] 处理 {year:04d}-{month:02d}")
        logger.info(f"{'=' * 80}")
        compute_month_wind_cf(
            data_dir=args.data_dir,
            year=year,
            month=month,
            output_dir=args.output_dir,
            chunk_time=args.chunk_time,
            overwrite=args.overwrite,
            compress_level=args.compress_level,
            ws_curve=ws_curve,
            pw_curve=pw_curve,
            rated_kw=rated_kw,
            extrap_ratio=extrap_ratio,
        )

    logger.info("\n✓ 全部完成！")


if __name__ == "__main__":
    main()

"""
使用方法示例：
python S02E02_Simulate_Wind_CF_ERA5Land.py --data_dir data --years 2015-2025 
"""
