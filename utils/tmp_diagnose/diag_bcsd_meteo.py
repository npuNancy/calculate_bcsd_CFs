"""诊断：统计 BCSD 各国 2015 年陆地格点上的平均风速 / 平均 rsds。

目的：判断 ERA5-Land CF 远高于 BCSD CF，是输入气象本身的差异，
还是 CF 算法的差异。这里只看 BCSD 输入气象（本地有数据）。

对每个国家：
  - 读 uas/vas/rsds 的 ssp126 2015 年（前 2920 个 3h 步）
  - 用 global_land_mask 构建陆地掩膜（经度转 [-180,180]）
  - 报告陆地格点上：平均 10m 风速 mean(sqrt(u^2+v^2))、平均 rsds
  - 同时报告陆地格点占比、rsds 的量级（判断单位 W/m2 vs J/m2）
"""
import sys
import numpy as np
import xarray as xr
from global_land_mask import globe

ROOT = "data/bcsd_outputs/NESM3"
SCN = "ssp126"
COUNTRIES = ["Germany", "United-Kingdom", "Denmark", "Italy", "Japan",
             "Austria", "France", "Spain", "Australia", "Egypt"]
N_2015 = 2920  # 3h 步/年（非闰年）


def land_mask(lat, lon):
    lon2 = ((lon + 180.0) % 360.0) - 180.0
    lat2 = np.clip(lat, -90.0, 90.0)
    lon2d, lat2d = np.meshgrid(lon2, lat2)
    return globe.is_land(lat2d, lon2d)  # (nlat, nlon)


import glob
def open_var(country, var):
    cands = sorted(glob.glob(f"{ROOT}/{country}/NESM3/{var}_3h_bcsd_on_0p1deg_{country}_NESM3_{SCN}_*.nc"))
    if not cands:
        raise FileNotFoundError(f"{country}/{var}")
    p = cands[0]
    ds = xr.open_dataset(p)
    da = ds[list(ds.data_vars)[0]]
    # 统一维度顺序 time,lat,lon
    da = da.transpose("time", "lat", "lon")
    return ds, da


print(f"{'country':16s} {'land%':>6s} {'meanWS10':>9s} {'meanRSDS':>9s} {'rsds_max':>9s} {'ws_max':>7s}")
for c in COUNTRIES:
    try:
        dsu, u = open_var(c, "uas")
        dsv, v = open_var(c, "vas")
        dsr, r = open_var(c, "rsds")
        lat = u["lat"].values.astype("float64")
        lon = u["lon"].values.astype("float64")
        lm = land_mask(lat, lon)  # (nlat,nlon)
        nland = int(lm.sum())
        landfrac = nland / lm.size

        us = u.isel(time=slice(0, N_2015)).values
        vs = v.isel(time=slice(0, N_2015)).values
        rs = r.isel(time=slice(0, N_2015)).values
        ws = np.sqrt(us**2 + vs**2)  # (t,lat,lon)

        lm3 = np.broadcast_to(lm, ws.shape)
        ws_land = ws[lm3]
        rs_land = rs[lm3]
        mean_ws = np.nanmean(ws_land)
        mean_rsds = np.nanmean(rs_land)
        rsds_max = np.nanmax(rs_land)
        ws_max = np.nanmax(ws_land)
        print(f"{c:16s} {landfrac*100:6.1f} {mean_ws:9.3f} {mean_rsds:9.2f} {rsds_max:9.1f} {ws_max:7.2f}")
        for d in (dsu, dsv, dsr):
            d.close()
    except Exception as e:
        print(f"{c:16s} ERROR {e}")
