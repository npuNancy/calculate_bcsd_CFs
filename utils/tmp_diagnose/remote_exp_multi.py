"""远程多国实验：对比 ERA5-Land 与 BCSD 的输入气象与逐小时 CF（2015 年）。

已确认时间分辨率不是主因（见 remote_exp_germany.py），故此处只算逐小时基准，
重点对比：ERA5 vs BCSD 的 平均/最大 10m 风速、平均 GHI、以及 solar/wind CF。
"""
import sys, importlib
import numpy as np
import xarray as xr

sys.path.insert(0, ".")
SOL = importlib.import_module("S02E01_Simulate_Solar_CF_ERA5Land")
WIN = importlib.import_module("S02E02_Simulate_Wind_CF_ERA5Land")

YEAR = 2015
# BCSD 名 -> (ERA5 目录名, BCSD: solar_cf, wind_cf, meanWS10, maxWS10, mean_rsds)
COUNTRIES = {
    "Germany":        ("Germany",        0.1318, 0.0814, 2.61, 16.9, 148.0),
    "United-Kingdom": ("United Kingdom", 0.0102, 0.0186, 3.76, 12.66, 123.4),
    "Italy":          ("Italy",          0.0472, 0.0060, 1.63, 9.19, 183.2),
    "Spain":          ("Spain",          0.1105, 0.0867, 2.94, 18.87, 208.9),
    "Austria":        ("Austria",        0.1957, 0.0268, 1.42, 9.58, 166.1),
    "Japan":          ("Japan",          0.0375, 0.0106, 2.05, 12.34, 200.6),
}

ws_curve, pw_curve, rated = WIN._get_power_curve_arrays()
ratio = WIN.power_law_ratio()


def load_year(base, var):
    das = []
    for m in range(1, 13):
        ds = xr.open_dataset(f"{base}/{var}/{var}_{YEAR}_{m:02d}.nc")
        tn = next(c for c in ["valid_time", "time"] if c in ds.coords or c in ds.dims)
        latn = next(c for c in ["latitude", "lat"] if c in ds.coords or c in ds.dims)
        lonn = next(c for c in ["longitude", "lon"] if c in ds.coords or c in ds.dims)
        v = list(ds.data_vars)[0]
        da = ds[v].rename({tn: "time", latn: "lat", lonn: "lon"}).transpose("time", "lat", "lon")
        das.append(da)
    return xr.concat(das, dim="time")


print(f"{'country':15s} | {'WS10 mean(E/B)':>16s} {'WS10 max(E/B)':>15s} {'GHI mean(E/B)':>15s} | "
      f"{'solar(E/B)':>15s} {'wind(E/B)':>15s}")
print("-" * 110)

for bname, (ename, b_sol, b_win, b_ws, b_wsmax, b_rsds) in COUNTRIES.items():
    base = f"data/ERA5_land/{ename}"
    try:
        ssrd = load_year(base, "ssrd"); t2m = load_year(base, "t2m")
        u10 = load_year(base, "u10"); v10 = load_year(base, "v10")
    except Exception as e:
        print(f"{bname:15s} | ERROR {e}")
        continue
    lats = ssrd["lat"].values.astype("float64"); lons = ssrd["lon"].values.astype("float64")
    N = ssrd.sizes["time"]
    land = SOL.build_land_mask(lats, lons)

    inc = SOL.accumulated_to_hourly_increment_chunk(ssrd, "time", 0, N, None, missing_boundary_policy="zero")
    ghi_kw = SOL.ssrd_increment_to_ghi_kw(inc)
    tas_c = SOL.tas_to_celsius(t2m.values, t2m.attrs.get("units"))
    u = u10.values.astype("float32"); v = v10.values.astype("float32")
    doy, hod = SOL.time_features(ssrd["time"])

    cf_sol = SOL.compute_solar_cf_chunk(ghi_kw, tas_c, u, v, doy, hod, lats, lons)
    cf_win = WIN.compute_wind_cf_chunk(u, v, ws_curve, pw_curve, rated, ratio)
    ws10 = np.sqrt(u**2 + v**2)

    def smean(cf):
        cell = np.nanmean(cf, axis=0)
        return float(np.nanmean(np.where(land, cell, np.nan)))

    e_ws = float(np.nanmean(ws10[:, land])); e_wsmax = float(np.nanmax(ws10[:, land]))
    e_ghi = float(np.nanmean((ghi_kw * 1000.0)[:, land]))
    e_sol = smean(cf_sol); e_win = smean(cf_win)

    print(f"{bname:15s} | {e_ws:6.2f}/{b_ws:<6.2f}   {e_wsmax:6.1f}/{b_wsmax:<6.1f}  "
          f"{e_ghi:6.1f}/{b_rsds:<6.1f} | {e_sol:6.3f}/{b_sol:<6.3f}  {e_win:6.3f}/{b_win:<6.3f}")
