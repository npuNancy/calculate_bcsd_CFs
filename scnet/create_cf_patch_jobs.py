"""Canonical SCNet entry point; implementation is kept with the patchify guide."""
from pathlib import Path
import runpy

if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).resolve().parents[1] / "infos" / "scnet_patchify" / "create_cf_patch_jobs.py"), run_name="__main__")

