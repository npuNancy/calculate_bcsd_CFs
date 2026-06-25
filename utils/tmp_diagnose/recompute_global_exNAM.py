"""重算 ERA5-Land 全球加权均值（剔除 NAM-12），并与 BCSD NESM3 全球均值对比。

ERA5-Land per_country 比 BCSD 多一个 NAM-12 区域。全球均值按陆地格点数加权：
    global = Σ_r mean_cf_r * n_land_cells_r / Σ_r n_land_cells_r
剔除 NAM-12 后重算，使两边区域集合一致（均为 27 国 + china... 实际同一批国家）。
"""
import csv
from collections import defaultdict

ERA5 = {
    "solar": "output/annual_mean_cf_ERA5Land/solar/per_country_annual_cf_ERA5Land.csv",
    "wind":  "output/annual_mean_cf_ERA5Land/wind/per_country_annual_cf_ERA5Land.csv",
}
BCSD_GLOBAL = "output/annual_mean_cf/NESM3/global_annual_cf_NESM3.csv"
EXCLUDE = {"NAM-12"}


def era5_global(path, exclude):
    # 列: energy,source,region,year,mean_cf,n_timesteps,n_land_cells
    num = defaultdict(float)  # year -> Σ cf*n
    den = defaultdict(float)  # year -> Σ n
    regions = set()
    excl_num = defaultdict(float); excl_den = defaultdict(float)
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            reg = row["region"].strip()
            year = int(row["year"])
            cf = float(row["mean_cf"]); n = float(row["n_land_cells"])
            if reg in exclude:
                excl_num[year] += cf*n; excl_den[year] += n
                continue
            regions.add(reg)
            num[year] += cf*n; den[year] += n
    years = sorted(num)
    out = {}
    for y in years:
        g_ex = num[y]/den[y]
        g_all = (num[y]+excl_num.get(y,0))/(den[y]+excl_den.get(y,0))
        out[y] = (g_all, g_ex, int(den[y]))
    return out, sorted(regions)


def bcsd_global(path, energy, scenario):
    out = {}
    with open(path) as f:
        r = csv.reader(f)
        header = next(r)
        for row in r:
            row = [c.strip() for c in row]
            if row[0] == energy and row[2] == scenario:
                out[int(row[3])] = (float(row[4]), int(row[6]))
    return out


for energy in ["solar", "wind"]:
    e, regs = era5_global(ERA5[energy], EXCLUDE)
    b = bcsd_global(BCSD_GLOBAL, energy, "ssp126")
    print("="*78)
    print(f"{energy.upper()}   ERA5-Land 区域数(剔NAM)={len(regs)}")
    print(f"{'year':>5} {'ERA5_all':>9} {'ERA5_exNAM':>11} {'BCSD_ssp126':>12} {'ERA5exNAM/BCSD':>15}")
    for y in sorted(e):
        g_all, g_ex, n = e[y]
        bc = b.get(y, (float('nan'),0))[0]
        ratio = g_ex/bc if bc==bc and bc>0 else float('nan')
        print(f"{y:>5} {g_all:9.4f} {g_ex:11.4f} {bc:12.4f} {ratio:15.2f}")
