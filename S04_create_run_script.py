"""
生成 CORDEX NAM-12 光伏/风电容量因子批量计算脚本
=================================================

扫描 data/CORDEX-CMIP6/NAM-12/1hr/ 下所有可用的
gcm_model / realization / rcm_model / scenario 组合，
为每种组合生成 solar 和 wind 计算命令，输出为 shell 脚本。

用法：
  python S04_create_run_script.py
  # 生成 run_all_solar_NAM-12.sh 和 run_all_wind_NAM-12.sh
"""

from __future__ import annotations

import os
from pathlib import Path

DATA_DIR = "data/CORDEX-CMIP6/NAM-12/1hr"
YEARS = "2020-2060"
CHUNK_TIME = 24


def discover_combinations(data_dir: str) -> list[tuple[str, str, str, str]]:
    """扫描目录，返回所有 (gcm_model, realization, rcm_model, scenario) 组合。"""
    combos: list[tuple[str, str, str, str]] = []
    root = Path(data_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"数据目录不存在：{data_dir}")

    for gcm in sorted(os.listdir(root)):
        gcm_path = root / gcm
        if not gcm_path.is_dir() or gcm.startswith("_"):
            continue
        for real in sorted(os.listdir(gcm_path)):
            real_path = gcm_path / real
            if not real_path.is_dir():
                continue
            for rcm in sorted(os.listdir(real_path)):
                rcm_path = real_path / rcm
                if not rcm_path.is_dir():
                    continue
                for scenario in sorted(os.listdir(rcm_path)):
                    scenario_path = rcm_path / scenario
                    if not scenario_path.is_dir() or scenario.startswith("_") or scenario.endswith("_"):
                        continue
                    combos.append((gcm, real, rcm, scenario))

    return combos


def main() -> None:
    combos = discover_combinations(DATA_DIR)
    if not combos:
        print("未发现任何数据组合。")
        return

    print(f"发现 {len(combos)} 个组合：")
    for gcm, real, rcm, scenario in combos:
        print(f"  {gcm} / {real} / {rcm} / {scenario}")
    print()

    solar_cmds: list[str] = []
    wind_cmds: list[str] = []

    for gcm, real, rcm, scenario in combos:
        common_args = (
            f"--gcm_model {gcm} --realization {real} "
            f"--rcm_model {rcm} --scenario {scenario} "
            f"--years {YEARS} --chunk_time {CHUNK_TIME}"
        )
        solar_cmds.append(
            f"python S04E01_Simulate_Solar_CF_NAM-12.py {common_args}"
        )
        wind_cmds.append(
            f"python S04E02_Simulate_Wind_CF_NAM-12.py {common_args}"
        )

    # 写入 shell 脚本
    for fname, cmds in [
        ("run_all_solar_NAM-12.sh", solar_cmds),
        ("run_all_wind_NAM-12.sh", wind_cmds),
    ]:
        with open(fname, "w") as f:
            f.write("#!/bin/bash\n\n")
            f.write("set -e\nset -x\n\n")
            f.write(f"# 自动生成，共 {len(cmds)} 个任务\n")
            f.write(f"# data_dir={DATA_DIR}, years={YEARS}, chunk_time={CHUNK_TIME}\n\n")
            for cmd in cmds:
                f.write(cmd + "\n")
        print(f"已生成 {fname}（{len(cmds)} 条命令）")


if __name__ == "__main__":
    main()
