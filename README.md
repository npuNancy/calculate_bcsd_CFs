# BCSD 风光容量因子计算脚本

本目录包含两个用于计算 BCSD 气候数据风光容量因子的 Python 脚本：

- `S01E01_Simulate_Solar_CF_BCSD.py`：计算光伏容量因子。
- `S01E02_Simulate_Wind_CF_BCSD.py`：计算陆上风电容量因子。

## 1. Python 环境

建议使用 Python 3.10 或更高版本。

安装依赖：

```bash
pip install -r requirements.txt -i https://mirrors.aliyun.com/pypi/simple/
```

## 2. 输入数据目录

两个脚本都假设输入数据采用如下结构：

```text
data/bcsd_outputs/{model}/{region}/{model}/*.nc
```

例如：

```text
data/bcsd_outputs/MIROC-ES2H/Austria/MIROC-ES2H/
├── pr_3h_bcsd_on_0p1deg_Austria_MIROC-ES2H_ssp126_2015-2060.nc
├── rsds_3h_bcsd_on_0p1deg_Austria_MIROC-ES2H_ssp126_2015-2060.nc
├── tas_3h_bcsd_on_0p1deg_Austria_MIROC-ES2H_ssp126_2015-2060.nc
├── uas_3h_bcsd_on_0p1deg_Austria_MIROC-ES2H_ssp126_2015-2060.nc
└── vas_3h_bcsd_on_0p1deg_Austria_MIROC-ES2H_ssp126_2015-2060.nc
```

每个 NetCDF 文件中通常只有一个变量，例如：

```text
pr_bcsd, rsds_bcsd, tas_bcsd, uas_bcsd, vas_bcsd
```

## 3. 光伏容量因子计算

光伏脚本使用：

- `rsds`：向下短波辐射。
- `tas`：近地面气温。
- `uas`、`vas`：近地面风速分量，用于组件温度修正。

运行示例：

```bash
python S01E01_Simulate_Solar_CF_BCSD.py \
  --data_dir data/bcsd_outputs \
  --model MIROC-ES2H \
  --region Austria \
  --scenario ssp126 \
  --years 2020-2020 
```

输出示例：

```text
output/CFs_of_solar/MIROC-ES2H/Austria/solar_CF_Austria_MIROC-ES2H_ssp126_2015-2060_allmonths.nc
```

输出变量：

```text
solar_cf
```

## 4. 风电容量因子计算

风电脚本使用：

- `uas`：10 m 东西向风速分量。
- `vas`：10 m 南北向风速分量。

当前风电参数：

- 风机类型：`GE120/2500`。
- 轮毂高度：`100 m`。
- 风速外推：幂律外推，`alpha = 1/7`。
- 切出风速：`25 m/s`。
- 不再额外乘风场效率 `0.9`。
- 只计算陆上风电容量因子。

运行示例：

```bash
python S01E02_Simulate_Wind_CF_BCSD.py \
  --data_dir data/bcsd_outputs \
  --model MIROC-ES2H \
  --region Austria \
  --scenario ssp126 \
  --years 2015-2060
```

输出示例：

```text
output/CFs_of_wind/MIROC-ES2H/Austria/wind_CF_Austria_MIROC-ES2H_ssp126_2015-2060_allmonths.nc
```

输出变量：

```text
wind_cf
```

## 5. 常用参数

两个脚本的参数基本一致：

| 参数 | 说明 | 示例 |
|---|---|---|
| `--data_dir` | BCSD 数据根目录 | `data/bcsd_outputs` |
| `--model` | 气候模型名称 | `MIROC-ES2H` |
| `--region` | 区域名称，也可以设为 `all` | `Austria` |
| `--scenario` | 情景名称 | `ssp126` |
| `--years` | 年份范围 | `2015-2060` |
| `--months` | 月份列表，空值表示全年 | `1,2,7,8` |
| `--output_dir` | 输出目录 | `output/CFs_of_wind` |
| `--chunk_time` | 每次处理的时间步数 | `512` |
| `--compress_level` | NetCDF 压缩等级 | `4` |
| `--overwrite` | 覆盖已有输出文件 | 无需额外值 |

处理某个模型下所有有效区域：

```bash
python S01E02_Simulate_Wind_CF_BCSD.py \
  --data_dir data/bcsd_outputs \
  --model MIROC-ES2H \
  --region all \
  --scenario ssp126 \
  --years 2015-2060 \
  --output_dir output/CFs_of_wind
```

## 6. 注意事项

1. `region=all` 时，脚本会自动跳过以下目录：
   - 以 `_` 开头的目录；
   - 以 `_` 结尾的目录；
   - 以 `_repeated` 结尾的目录。

2. 脚本不会把 3 小时数据插值到小时。输出时间分辨率与输入一致。

3. 若输出文件已存在，默认跳过。需要重算时加 `--overwrite`。

4. 如果 `windpowerlib` 找不到 `GE120/2500`，请先尝试升级：

```bash
pip install -U windpowerlib -i https://mirrors.aliyun.com/pypi/simple/
```
