# global_bcsd patchify CF SCNet goal

- Canonical generator: `scnet/create_cf_patch_jobs.py`; it only writes scripts and a manifest.
- One unit: `model × scenario × patch × tech`. The production patch list is supplied explicitly from the
  global_bcsd manifest (currently 47 active patches); never infer a 5×12 grid.
- Runtime command: `patchify_station_cf.py`; output is station-only NetCDF plus a JSON sidecar.
- Required sidecar identity: `patch_id` and `variable` must match every input file.
- Runtime environment: `source /work/home/acbpgywfpz/miniconda3/bin/activate climate`.
- Job/log directories are external to the repository. The generator does not call `sbatch`.

