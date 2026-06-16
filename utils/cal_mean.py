from pathlib import Path

import numpy as np
import xarray as xr
import pandas as pd
import dask
from dask.distributed import Client, LocalCluster
from global_land_mask import globe


def build_land_mask(ds):
    """根据数据集的 lat/lon 坐标，用 global-land-mask 构建陆地掩膜。

    返回一个 (lat, lon) 的布尔 DataArray，True 表示陆地。
    注意：global_land_mask 要求经度在 [-180, 180]，
    而 ERA5-Land 经度通常是 [0, 360)，因此需要先转换。
    """
    lats = ds["lat"].values.astype("float64")
    lons = ds["lon"].values.astype("float64")

    # 经度归一化到 [-180, 180]，仅用于查询掩膜，不改变数据本身的坐标。
    lons_conv = ((lons + 180.0) % 360.0) - 180.0

    # 纬度也夹紧到 [-90, 90]，避免浮点误差（如 90.00000001）触发越界。
    lats_clip = np.clip(lats, -90.0, 90.0)

    lon2d, lat2d = np.meshgrid(lons_conv, lats_clip)  # 均为 (n_lat, n_lon)
    land = globe.is_land(lat2d, lon2d)  # (n_lat, n_lon) 布尔数组

    return xr.DataArray(
        land,
        dims=("lat", "lon"),
        coords={"lat": ds["lat"], "lon": ds["lon"]},
        name="land_mask",
    )


def _global_and_land_means(da, land_mask):
    """返回 (全球均值任务, 仅陆地均值任务)。

    - 全球均值：对所有有效格点求平均（海洋格点当前为 0，会被计入，
      与原始口径一致）。
    - 仅陆地均值：先用 land_mask 把海洋格点置为 NaN，再 skipna 求平均。
    """
    dims = ("time", "lat", "lon")
    global_task = da.mean(dim=dims, skipna=True)
    land_task = da.where(land_mask).mean(dim=dims, skipna=True)
    return global_task, land_task


def calc_one_month(file_wind, file_solar, year, month):
    ds_wind = xr.open_dataset(
        file_wind,
        chunks={"time": 72},  # 匹配原始 NetCDF chunk
        cache=False,
    )

    ds_solar = xr.open_dataset(
        file_solar,
        chunks={"time": 72},
        cache=False,
    )

    # 风、光使用各自文件的 lat/lon 构建掩膜（两者网格通常一致）。
    wind_land_mask = build_land_mask(ds_wind)
    solar_land_mask = build_land_mask(ds_solar)

    wind_global_task, wind_land_task = _global_and_land_means(ds_wind["wind_cf"], wind_land_mask)
    solar_global_task, solar_land_task = _global_and_land_means(ds_solar["solar_cf"], solar_land_mask)

    wind_global, wind_land, solar_global, solar_land = dask.compute(
        wind_global_task,
        wind_land_task,
        solar_global_task,
        solar_land_task,
    )

    ds_wind.close()
    ds_solar.close()

    return {
        "year": year,
        "month": month,
        "wind_cf_global": float(wind_global),
        "wind_cf_land": float(wind_land),
        "solar_cf_global": float(solar_global),
        "solar_cf_land": float(solar_land),
    }


def main_1():
    year = 2015

    wind_dir = Path(
        "/data1/yanxiaokai/project_climate/calculate_bcsd_cfs/" "output/mean_CFs_of_wind_ERA5Land/2015_2015/monthly"
    )

    solar_dir = Path(
        "/data1/yanxiaokai/project_climate/calculate_bcsd_cfs/" "output/mean_CFs_of_solar_ERA5Land/2015_2015/monthly"
    )

    # 注意：单个 chunk 解压后约 1.87 GB；
    # skipna=True 会产生额外临时数组，所以不要开太多 worker。
    cluster = LocalCluster(
        n_workers=1,
        threads_per_worker=1,
        memory_limit="24GB",
        dashboard_address=":0",  # 避免 8787 端口冲突
    )

    client = Client(cluster)
    print(client)

    records = []

    for month in range(1, 13):
        file_wind = wind_dir / f"mean_wind_cf_{month:02d}.nc"
        file_solar = solar_dir / f"mean_solar_cf_{month:02d}.nc"

        print(f"Processing {year}-{month:02d}")
        print("  wind :", file_wind)
        print("  solar:", file_solar)

        rec = calc_one_month(file_wind, file_solar, year, month)
        records.append(rec)

        print("  result:", rec)

    df = pd.DataFrame(records)
    out_csv = f"global_monthly_mean_cf_{year}.csv"
    df.to_csv(out_csv, index=False)

    print(df)
    print(f"Saved to {out_csv}")

    client.close()
    cluster.close()


def main_2():
    pass

    file_wind = "output/mean_CFs_of_wind_ERA5Land/2015_2015/mean_wind_cf_2015_2015.nc"
    file_solar = "output/mean_CFs_of_solar_ERA5Land/2015_2015/mean_solar_cf_2015_2015.nc"

    print("  wind :", file_wind)
    print("  solar:", file_solar)

    rec = calc_one_month(file_wind, file_solar, 2015, "full-year")
    print("  result:", rec)
    print(f"  风电 CF: 全球(含海洋0)={rec['wind_cf_global']:.4f}  仅陆地={rec['wind_cf_land']:.4f}")
    print(f"  光伏 CF: 全球(含海洋0)={rec['solar_cf_global']:.4f}  仅陆地={rec['solar_cf_land']:.4f}")


if __name__ == "__main__":
    main_1()
    # main_2()
