"""远程实验（在 THU_Spark02 的 climate conda 环境下运行）。

目的：把 "ERA5-Land CF 远高于 BCSD CF" 拆成两个独立因素：
  (A) 时间分辨率：ERA5 逐小时 vs BCSD 逐 3 小时；
  (B) 数据源：ERA5-Land 再分析 vs BCSD 降尺度 GCM 的气象本身差异。

做法（以 Germany 2015 为例，用 ERA5-Land 本国月文件）：
  1. 用仓库原函数（S02E01 solar / S02E02 wind）在逐小时下算 CF；
  2. 把同一份 ERA5 输入降采样到 3 小时（抽样，相位 0/1/2）与 3 小时块平均，
     用同一算法重算 CF —— 模拟 BCSD 的时间分辨率；
  3. 报告 ERA5 Germany 的平均/最大 10m 风速、平均 GHI，用于和 BCSD 对比。

对照（本地 BCSD Germany 2015, ssp126）：
  solar CF=0.132, wind CF=0.0814, mean WS10=2.61 m/s, max WS10=16.9, mean rsds=148 W/m2
"""
import sys
import numpy as np
import xarray as xr

sys.path.insert(0, ".")
import importlib
SOL = importlib.import_module("S02E01_Simulate_Solar_CF_ERA5Land")
WIN = importlib.import_module("S02E02_Simulate_Wind_CF_ERA5Land")

DATA = "data/ERA5_land/Germany"
YEAR = 2015


def load_year(var):
    das = []
    for m in range(1, 13):
        ds = xr.open_dataset(f"{DATA}/{var}/{var}_{YEAR}_{m:02d}.nc")
        tn = next(c for c in ["valid_time", "time"] if c in ds.coords or c in ds.dims)
        latn = next(c for c in ["latitude", "lat"] if c in ds.coords or c in ds.dims)
        lonn = next(c for c in ["longitude", "lon"] if c in ds.coords or c in ds.dims)
        v = list(ds.data_vars)[0]
        da = ds[v].rename({tn: "time", latn: "lat", lonn: "lon"}).transpose("time", "lat", "lon")
        das.append(da)
    return xr.concat(das, dim="time")


print("加载 Germany 2015 逐小时 ssrd/t2m/u10/v10 ...")
ssrd = load_year("ssrd")
t2m = load_year("t2m")
u10 = load_year("u10")
v10 = load_year("v10")
lats = ssrd["lat"].values.astype("float64")
lons = ssrd["lon"].values.astype("float64")
N = ssrd.sizes["time"]
print(f"  时间步 N={N}  网格 {len(lats)}x{len(lons)}")

land = SOL.build_land_mask(lats, lons)  # (lat,lon) True=陆地
nland = int(land.sum())

# ── 去累积得到逐小时 GHI ──
inc = SOL.accumulated_to_hourly_increment_chunk(ssrd, "time", 0, N, None, missing_boundary_policy="zero")
ghi_kw = SOL.ssrd_increment_to_ghi_kw(inc)                 # (N,lat,lon) kW/m2
tas_c = SOL.tas_to_celsius(t2m.values, t2m.attrs.get("units"))
u = u10.values.astype("float32"); v = v10.values.astype("float32")

time_da = ssrd["time"]
doy_h, hod_h = SOL.time_features(time_da)

# 功率曲线（风电）
ws_curve, pw_curve, rated = WIN._get_power_curve_arrays()
ratio = WIN.power_law_ratio()


def spatial_mean_cf(cf_txy):
    """时间平均 -> 陆地掩膜 -> 空间平均。"""
    cell = np.nanmean(cf_txy, axis=0)          # (lat,lon)
    cell = np.where(land, cell, np.nan)
    return float(np.nanmean(cell))


def solar_cf(ghi, tas, uu, vv, doy, hod):
    return SOL.compute_solar_cf_chunk(ghi, tas, uu, vv, doy, hod, lats, lons)


def wind_cf(uu, vv):
    return WIN.compute_wind_cf_chunk(uu, vv, ws_curve, pw_curve, rated, ratio)


# ── 逐小时基准 ──
cf_sol_h = solar_cf(ghi_kw, tas_c, u, v, doy_h, hod_h)
cf_win_h = wind_cf(u, v)
ws10 = np.sqrt(u**2 + v**2)

print("\n========== ERA5-Land Germany 2015 气象 (陆地) ==========")
ws_land = ws10[:, land]
ghi_land = (ghi_kw * 1000.0)[:, land]   # W/m2（含夜间）
print(f"  mean WS10 = {np.nanmean(ws_land):.3f} m/s   max WS10 = {np.nanmax(ws_land):.2f}")
print(f"  mean GHI  = {np.nanmean(ghi_land):.2f} W/m2  max GHI = {np.nanmax(ghi_land):.1f}")
print(f"  (对照 BCSD Germany: mean WS10=2.61, max=16.9, mean rsds=148)")


def degrade_report(name, doy_s, hod_s, ghi_s, tas_s, u_s, v_s):
    cs = solar_cf(ghi_s, tas_s, u_s, v_s, doy_s, hod_s)
    cw = wind_cf(u_s, v_s)
    print(f"  {name:22s} solar={spatial_mean_cf(cs):.4f}   wind={spatial_mean_cf(cw):.4f}")


print("\n========== CF：逐小时 vs 降到 3 小时 ==========")
print(f"  {'hourly (基准)':22s} solar={spatial_mean_cf(cf_sol_h):.4f}   wind={spatial_mean_cf(cf_win_h):.4f}")

# 3 小时抽样，相位 0/1/2
for off in (0, 1, 2):
    idx = np.arange(off, N, 3)
    degrade_report(f"3h-subsample off={off}", doy_h[idx], hod_h[idx],
                   ghi_kw[idx], tas_c[idx], u[idx], v[idx])

# 3 小时块平均（标签取块中点）
M = (N // 3) * 3
def block_mean(a):
    return a[:M].reshape(M // 3, 3, *a.shape[1:]).mean(axis=1)
ghi_b = block_mean(ghi_kw); tas_b = block_mean(tas_c)
u_b = block_mean(u); v_b = block_mean(v)
doy_b = doy_h[:M].reshape(M // 3, 3).mean(axis=1)
hod_b = hod_h[:M].reshape(M // 3, 3).mean(axis=1)
degrade_report("3h-blockmean", doy_b, hod_b, ghi_b, tas_b, u_b, v_b)

print("\n  (对照 BCSD Germany 2015 ssp126: solar=0.132, wind=0.0814)")
