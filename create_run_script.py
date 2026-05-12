#  生成 所有 区域×情景 的 光伏/风电 容量因子 计算

import os


def collect_regions(root: str, model: str) -> list[str]:
    # 收集指定模型下的所有区域
    # region 不以 '_' 开头/结尾，不以 '_repeated' 结尾
    model_dir = os.path.join(root, model)
    if not os.path.isdir(model_dir):
        return []
    regions = []
    for entry in sorted(os.listdir(model_dir)):
        if entry.startswith("_") or entry.endswith("_"):
            continue
        if entry.endswith("_repeated"):
            continue
        region_path = os.path.join(model_dir, entry, model)
        if os.path.isdir(region_path):
            regions.append(entry)
    return regions


data_dir = "data/bcsd_outputs"
model = "MIROC-ES2H"
scenario_list = ["ssp126", "ssp245", "ssp585"]

regions = collect_regions(data_dir, model)

cmd_wind_list = []
cmd_solar_list = []

# 遍历所有区域和情景，生成命令
for region in regions:
    for scenario in scenario_list:
        cmd_wind = f"python S01E02_Simulate_Wind_CF_BCSD.py --data_dir {data_dir} --model {model} --region {region} --scenario {scenario} --years 2015-2060"
        cmd_solar = f"python S01E01_Simulate_Solar_CF_BCSD.py --data_dir {data_dir} --model {model} --region {region} --scenario {scenario} --years 2015-2060"
        cmd_wind_list.append(cmd_wind)
        cmd_solar_list.append(cmd_solar)

with open("run_all_wind.sh", "w") as f_wind:
    f_wind.write("#!/bin/bash\n\n")
    f_wind.write("set -e\nset -x\n\n")
    for cmd in cmd_wind_list:
        f_wind.write(cmd + "\n")

with open("run_all_solar.sh", "w") as f_solar:
    f_solar.write("#!/bin/bash\n\n")
    f_solar.write("set -e\nset -x\n\n")
    for cmd in cmd_solar_list:
        f_solar.write(cmd + "\n")
