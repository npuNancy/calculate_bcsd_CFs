# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Scripts that compute **renewable-energy capacity factors (CF)** — solar PV and onshore
wind — from downscaled climate data (CMIP6 BCSD, CMIP6-CORDEX NAM-12) and ERA5-Land
reanalysis. Output is gridded NetCDF time series of `solar_cf` / `wind_cf` plus
derived multi-year and annual means. Not a package — a flat collection of standalone,
argparse-driven scripts. Code comments, docstrings, and commit messages are in Chinese.

## Environment

The `.venv` is **uv-managed Python 3.9 and has no pip**. Install packages with:

```bash
VIRTUAL_ENV=.venv uv pip install <pkg>
```

(The README says Python 3.10+ / `pip install -r requirements.txt`, but the actual
checked-out venv is 3.9 with uv — trust the venv.) Run scripts via `.venv/bin/python`.

## Script naming and pipeline

Scripts follow `S{stage}E{episode}_...` — stage = data source / processing phase,
episode = solar(E01) / wind(E02) / extra(E03). The big picture:

| Stage | Data source | Grid / notes |
|---|---|---|
| **S01** E01/E02 | CMIP6 **BCSD** per-country, 3-hourly | `(time, lat, lon)`; uses `uas`/`vas` |
| **S02** E01/E02 | **ERA5-Land** global, hourly, monthly files | regular grid; E03 = multi-year mean hourly CF |
| **S03** E01/E02 | BCSD **china** region | uses scalar `sfcWind` (not uas/vas) |
| **S04** E01/E02 | CMIP6-**CORDEX NAM-12**, hourly | rotated pole `(time, rlat, rlon)`, 2-D lat/lon |
| **S05** E01/E02/E03 | ERA5-Land climatology | mean 100 m wind / sunshine hours / temperature |
| **S06** E01/E02 | reduces existing CF files | E01 = BCSD per-country + global; E02 = ERA5-Land annual means |

Each CF script takes `--data_dir --model --region --scenario --years [--months]
[--output_dir] [--chunk_time] [--compress_level] [--overwrite]`. `--region all`
batch-processes every valid region under a model (skips dirs starting/ending in `_`
or ending in `_repeated`). Output defaults to `output/CFs_of_{solar,wind}*/...`.

### Running batches

`create_run_script.py`, `S03_create_run_script.py`, `S04_create_run_script.py`
**generate** the `run_all_{solar,wind}_{model}.sh` shell scripts (which loop over
regions × scenarios). Regenerate after data changes; then `bash run_all_*.sh`.
Edit the `model` / `scenario_list` constants at the top of the generator, not the
generated `.sh` files.

Example single run:

```bash
.venv/bin/python S01E01_Simulate_Solar_CF_BCSD.py \
  --data_dir data/bcsd_outputs --model MIROC-ES2H \
  --region Austria --scenario ssp126 --years 2015-2060
```

## Shared physics / conventions

- **Solar**: `rsds` → GHI (kW m⁻²), Erbs decomposition into DNI/diffuse, **dual-axis
  tracking** plane-of-array, module-temperature correction, system coefficient
  `SYS_COEF = 0.8056`. `tas`/wind are instantaneous and linearly interpolated onto the
  `rsds` half-hour time axis (rsds is a time-mean flux labelled at `:30`). CF clipped to `[0, 1]`.
- **Wind**: `windpowerlib` turbine `GE120/2500`, hub height 100 m, power-law
  extrapolation `alpha = 1/7`, cut-out 25 m s⁻¹, **no** farm-efficiency multiplier.
- Scripts process in **time chunks** and stream chunk-by-chunk into NetCDF to bound
  memory. Output is written to a `.tmp` file then `os.replace`d (atomic; partial/failed
  outputs are cleaned up). Existing outputs are skipped unless `--overwrite`.
- 3-hourly input is **not** interpolated to hourly; output resolution matches input.

## Critical gotcha: ocean cells

CF NetCDFs carry **meaningless values over ocean** (early versions filled NaN with 0
via `np.nan_to_num`; ERA5-Land is ~66% ocean, so unmasked spatial means are badly
diluted). **Always apply a land mask before any spatial mean.** Use
`global_land_mask` (note: it needs longitudes in `[-180, 180]`, but ERA5-Land uses
`[0, 360)` — convert first; see `utils/cal_mean.py: build_land_mask`). S06 masks
on-the-fly by default (`--no-land-mask` to disable); `utils/mask_ocean_era5land_cf.py`
is a post-hoc tool to rewrite ocean cells to NaN in existing files.

## Other helpers

- `era5land_utils.py` — shared ERA5-Land I/O: coord-name resolution
  (`valid_time`/`time`, `latitude`/`lat`, …), grid probing, atomic `save_result`.
- `utils/cal_mean.py` — dask-based land-masked spatial means.
- `utils/save_region_grids.py` — dumps per-region lat/lon grids to `data/grid_of_regions/`.
- `plot_cf.py` — Cartopy CF maps → `figs/CF_maps/`.
- `document/*.md` — Chinese design/bug-fix plans for past changes.

## Data layout

- BCSD input: `data/bcsd_outputs/{model}/{region}/{model}/{var}_3h_bcsd_*.nc`
  (one variable per file, named e.g. `rsds_bcsd`; falls back to the sole data var).
- ERA5-Land input: `{data_dir}/{var}/{var}_{YYYY}_{MM}.nc`.
- Outputs and `data/` are gitignored; only code, grids, and a sample of figures are tracked.


## 要求
- 每次回复我时，都称呼我为 `小凯`，并且在回复的最后加上 `希望对你有帮助，小凯！`(不需要加反引号)。
