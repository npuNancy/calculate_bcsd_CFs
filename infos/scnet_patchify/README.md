# CF patchify guide

The canonical launcher is `scnet/create_cf_patch_jobs.py`. Run it with explicit models, scenarios, patches, station
catalog, BCSD root, output and external job/log directories. It never calls
`sbatch`; submit and monitor under the shared SCNet lock. Each generated script
executes exactly one station-only CF unit.
Pass the active patch IDs explicitly with `--patches`; the generator intentionally has no
default 5×12 patch list.
