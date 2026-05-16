#  生成 中国区域所有模型×情景 的 光伏/风电 容量因子计算脚本

import os

china_data_dir = "data/cmip6_downscaling_3hr"

china_models = [
    # "BCC-CSM2-MR",
    # "CanESM5",
    "MIROC-ES2H",
    # "MIROC6",
    # "MPI-ESM1-2-HR",
    # "MRI-ESM2-0",
    # "NESM3",
]

# MIROC6 和 NESM3 缺少 ssp370
scenarios = {
    # "historical": {"years": "1979-2014"},
    # "ssp370": {"years": "2015-2060"},
    "ssp126": {"years": "2060-2060"},
    "ssp245": {"years": "2060-2060"},
    "ssp585": {"years": "2060-2060"},
}

# 不含 ssp370 的模型
no_ssp370 = {"MIROC6", "NESM3"}

cmd_wind_list = []
cmd_solar_list = []

for model in china_models:
    for scenario, info in scenarios.items():
        if model in no_ssp370 and scenario == "ssp370":
            continue
        years = info["years"]
        cmd_wind = (
            f"python S03E02_Simulate_Wind_CF_china.py "
            f"--data_dir {china_data_dir} "
            f"--model {model} "
            f"--scenario {scenario} "
            f"--years {years}"
        )
        cmd_solar = (
            f"python S03E01_Simulate_Solar_CF_china.py "
            f"--data_dir {china_data_dir} "
            f"--model {model} "
            f"--scenario {scenario} "
            f"--years {years}"
        )
        cmd_wind_list.append(cmd_wind)
        cmd_solar_list.append(cmd_solar)

with open("run_all_wind_china.sh", "w") as f:
    f.write("#!/bin/bash\n\n")
    f.write("set -e\nset -x\n\n")
    for cmd in cmd_wind_list:
        f.write(cmd + "\n")

with open("run_all_solar_china.sh", "w") as f:
    f.write("#!/bin/bash\n\n")
    f.write("set -e\nset -x\n\n")
    for cmd in cmd_solar_list:
        f.write(cmd + "\n")

print(f"已生成 run_all_wind_china.sh（{len(cmd_wind_list)} 条命令）")
print(f"已生成 run_all_solar_china.sh（{len(cmd_solar_list)} 条命令）")
