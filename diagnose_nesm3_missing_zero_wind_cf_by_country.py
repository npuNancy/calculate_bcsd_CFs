#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diagnose pseudo-zero wind CF caused by missing BCSD meteorology.

This script is intentionally narrow in model/scenario/technology, but can run
multiple regions:
  - model: NESM3
  - technology: wind
  - default scenario: ssp126

It compares two means from the same existing wind CF output:
  1. current_mean_cf: current CF mean, where source BCSD NaNs may have already
     become zeros during CF simulation.
  2. source_valid_mean_cf: the same CF array after dropping grid-time points
     where source uas or vas is NaN.

No large NetCDF is written. The output is a small CSV.
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


ROOT = Path(__file__).parent
DEFAULT_BCSD_DIR = ROOT / "data" / "bcsd_outputs"
DEFAULT_CF_DIR = ROOT / "output" / "CFs_of_wind"
DEFAULT_OUT = ROOT / "output" / "diagnostics" / "nesm3_missing_zero_wind_cf_by_country.csv"

MODEL = "NESM3"
DEFAULT_SCENARIOS = ("ssp126",)
DEFAULT_REGIONS = (
    "Turkey",
    "Italy",
    "United-Kingdom",
    "Germany",
    "Denmark",
    "Portugal",
)


def parse_years(text: str | None) -> list[int] | None:
    if not text:
        return None
    text = text.strip()
    if "-" in text:
        y0, y1 = map(int, text.split("-", 1))
        if y1 < y0:
            raise ValueError(f"Invalid --years: {text}")
        return list(range(y0, y1 + 1))
    return [int(x) for x in text.split(",") if x.strip()]


def find_one(pattern: str) -> Path:
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No file matched: {pattern}")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple files matched {pattern}: {matches}")
    return matches[0]


def data_var(ds: xr.Dataset, preferred: str) -> str:
    candidate = f"{preferred}_bcsd"
    if candidate in ds.data_vars:
        return candidate
    if preferred in ds.data_vars:
        return preferred
    if len(ds.data_vars) == 1:
        return next(iter(ds.data_vars))
    raise KeyError(f"Cannot identify variable {preferred}; data_vars={list(ds.data_vars)}")


def parse_csv_list(text: str) -> list[str]:
    return [x.strip() for x in text.split(",") if x.strip()]


def open_inputs(bcsd_dir: Path, cf_dir: Path, region: str, scenario: str, chunk_time: int):
    bcsd_base = bcsd_dir / MODEL / region / MODEL
    uas_path = find_one(str(bcsd_base / f"uas_3h_bcsd_on_0p1deg_{region}_{MODEL}_{scenario}_*.nc"))
    vas_path = find_one(str(bcsd_base / f"vas_3h_bcsd_on_0p1deg_{region}_{MODEL}_{scenario}_*.nc"))
    cf_path = find_one(str(cf_dir / MODEL / region / f"wind_CF_{region}_{MODEL}_{scenario}_*_allmonths.nc"))

    chunks = {"time": chunk_time}
    ds_uas = xr.open_dataset(uas_path, chunks=chunks)
    ds_vas = xr.open_dataset(vas_path, chunks=chunks)
    ds_cf = xr.open_dataset(cf_path, chunks=chunks)
    return ds_uas, ds_vas, ds_cf, uas_path, vas_path, cf_path


def scalar(value) -> float:
    out = value.compute() if hasattr(value, "compute") else value
    return float(out.values if hasattr(out, "values") else out)


def int_scalar(value) -> int:
    return int(round(scalar(value)))


def summarize_period(
    region: str,
    scenario: str,
    period: str,
    cf: xr.DataArray,
    source_valid: xr.DataArray,
    uas: xr.DataArray,
    vas: xr.DataArray,
) -> dict:
    cf_domain = np.isfinite(cf)
    source_missing_in_cf_domain = cf_domain & (~source_valid)
    zero_cf = cf_domain & (cf == 0)
    zero_from_source_missing = zero_cf & (~source_valid)

    n_cf = int_scalar(cf_domain.sum())
    n_source_missing = int_scalar(source_missing_in_cf_domain.sum())
    n_zero = int_scalar(zero_cf.sum())
    n_zero_source_missing = int_scalar(zero_from_source_missing.sum())

    current_mean = scalar(cf.mean(skipna=True))
    source_valid_mean = scalar(cf.where(source_valid).mean(skipna=True))
    missing_uas_frac = scalar((~np.isfinite(uas)).mean())
    missing_vas_frac = scalar((~np.isfinite(vas)).mean())

    return {
        "region": region,
        "model": MODEL,
        "scenario": scenario,
        "period": period,
        "current_mean_cf": current_mean,
        "source_valid_mean_cf": source_valid_mean,
        "delta_cf": source_valid_mean - current_mean,
        "relative_increase_pct": (source_valid_mean / current_mean - 1.0) * 100.0
        if current_mean > 0
        else np.nan,
        "cf_valid_count": n_cf,
        "source_missing_in_cf_domain_count": n_source_missing,
        "source_missing_in_cf_domain_frac": n_source_missing / n_cf if n_cf else np.nan,
        "zero_cf_count": n_zero,
        "zero_cf_frac": n_zero / n_cf if n_cf else np.nan,
        "zero_from_source_missing_count": n_zero_source_missing,
        "zero_from_source_missing_frac_of_zero": n_zero_source_missing / n_zero if n_zero else np.nan,
        "uas_nan_frac_all_gridtime": missing_uas_frac,
        "vas_nan_frac_all_gridtime": missing_vas_frac,
    }


def scenario_rows(
    bcsd_dir: Path,
    cf_dir: Path,
    region: str,
    scenario: str,
    years: list[int] | None,
    chunk_time: int,
) -> list[dict]:
    ds_uas, ds_vas, ds_cf, uas_path, vas_path, cf_path = open_inputs(
        bcsd_dir, cf_dir, region, scenario, chunk_time
    )
    try:
        uas_name = data_var(ds_uas, "uas")
        vas_name = data_var(ds_vas, "vas")
        uas = ds_uas[uas_name]
        vas = ds_vas[vas_name]
        cf = ds_cf["wind_cf"]

        # Some source BCSD files cover a longer period than the produced CF file.
        # Compare only the exact timestamps present in the CF output.
        uas = uas.sel(time=cf.time)
        vas = vas.sel(time=cf.time)

        # Avoid coordinate alignment issues from float32/float64 lat/lon dtypes.
        uas = uas.assign_coords(time=cf.time, lat=cf.lat, lon=cf.lon)
        vas = vas.assign_coords(time=cf.time, lat=cf.lat, lon=cf.lon)
        source_valid = (np.isfinite(uas) & np.isfinite(vas)).assign_coords(
            time=cf.time, lat=cf.lat, lon=cf.lon
        )

        available_years = [int(y) for y in np.unique(cf.time.dt.year.values)]
        selected_years = years if years is not None else available_years

        rows: list[dict] = []
        for year in selected_years:
            if year not in available_years:
                continue
            year_sel = dict(time=str(year))
            rows.append(
                summarize_period(
                    region,
                    scenario,
                    str(year),
                    cf.sel(**year_sel),
                    source_valid.sel(**year_sel),
                    uas.sel(**year_sel),
                    vas.sel(**year_sel),
                )
            )

        if years is None or len(selected_years) > 1:
            rows.append(summarize_period(region, scenario, "all", cf, source_valid, uas, vas))

        for row in rows:
            row["cf_file"] = str(cf_path)
            row["uas_file"] = str(uas_path)
            row["vas_file"] = str(vas_path)
        return rows
    finally:
        ds_uas.close()
        ds_vas.close()
        ds_cf.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare NESM3 wind CF before/after dropping BCSD source NaNs by country."
    )
    parser.add_argument("--bcsd-dir", type=Path, default=DEFAULT_BCSD_DIR)
    parser.add_argument("--cf-dir", type=Path, default=DEFAULT_CF_DIR)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--regions", default=",".join(DEFAULT_REGIONS), help="Comma-separated country/region names.")
    parser.add_argument("--scenarios", default=",".join(DEFAULT_SCENARIOS), help="Comma-separated SSPs.")
    parser.add_argument("--years", default=None, help="YYYY, YYYY-YYYY, or comma-separated years. Default: all.")
    parser.add_argument("--chunk-time", type=int, default=512, help="xarray/dask time chunk size.")
    args = parser.parse_args()

    regions = parse_csv_list(args.regions)
    scenarios = parse_csv_list(args.scenarios)
    years = parse_years(args.years)

    rows: list[dict] = []
    for region in regions:
        for scenario in scenarios:
            print(f"[RUN] {region} {MODEL} {scenario}", flush=True)
            rows.extend(scenario_rows(args.bcsd_dir, args.cf_dir, region, scenario, years, args.chunk_time))

    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    cols = [
        "region",
        "scenario",
        "period",
        "current_mean_cf",
        "source_valid_mean_cf",
        "delta_cf",
        "relative_increase_pct",
        "source_missing_in_cf_domain_frac",
        "zero_cf_frac",
        "zero_from_source_missing_frac_of_zero",
    ]
    print(df[cols].to_string(index=False, float_format=lambda x: f"{x:.6f}"))
    print(f"[OK] wrote {args.out}")


if __name__ == "__main__":
    main()
