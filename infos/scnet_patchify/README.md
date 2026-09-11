# CF patchify guide

Run `create_cf_patch_jobs.py` with explicit models, scenarios, patches, station
catalog, BCSD root, output and external job/log directories. It never calls
`sbatch`; submit and monitor under the shared SCNet lock. Each generated script
executes exactly one station-only CF unit.
