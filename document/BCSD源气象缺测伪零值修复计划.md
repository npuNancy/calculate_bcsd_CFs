# BCSD 源气象缺测伪零值修复计划

## 背景

`S01E01_Simulate_Solar_CF_BCSD.py` 与 `S01E02_Simulate_Wind_CF_BCSD.py` 在单位转换或计算保护阶段会把输入气象中的 `NaN` 转为 `0`。脚本随后只用 `global_land_mask` 把海洋格点恢复为 `NaN`，但 BCSD 区域边界处的陆地缺测格点仍保留为 `CF = 0`。

这些 `0` 不是物理意义上的零出力，而是源气象缺测传播形成的伪零值，会拉低格点年均 CF、场站 CF 和下游 RQ1 统计。

## 修改目标

在保持现有 CF 计算公式、风机功率曲线、光伏温度修正和海洋掩膜逻辑不变的前提下：

1. 计算前记录源气象有效性掩膜。
2. 允许计算过程继续用 `nan_to_num` 防御数值异常。
3. 写出 NetCDF 前，将源气象缺测位置恢复为 `NaN`。
4. 同时继续将海洋格点设为 `NaN`。

## 修改范围

- `ref_code/calculate_bcsd_cfs/S01E01_Simulate_Solar_CF_BCSD.py`
- `ref_code/calculate_bcsd_cfs/S01E02_Simulate_Wind_CF_BCSD.py`

不修改 S02/S03/S04 系列，不修改 RQ1 绘图脚本。

## 具体方案

### 光伏 S01E01

在每个时间块读取或插值得到：

- `rsds_raw`
- `tas_raw`
- `uas_raw`
- `vas_raw`

后，先构造：

```python
source_valid = (
    np.isfinite(rsds_raw)
    & np.isfinite(tas_raw)
    & np.isfinite(uas_raw)
    & np.isfinite(vas_raw)
)
```

计算 `cf_chunk` 后写入前执行：

```python
valid_mask = source_valid & land_mask[None, :, :]
cf_chunk = np.where(valid_mask, cf_chunk, np.nan).astype(np.float32)
```

### 风电 S01E02

在每个时间块读取：

- `uas_raw`
- `vas_raw`

后，先构造：

```python
source_valid = np.isfinite(uas_raw) & np.isfinite(vas_raw)
```

计算 `cf_chunk` 后写入前执行：

```python
valid_mask = source_valid & land_mask[None, :, :]
cf_chunk = np.where(valid_mask, cf_chunk, np.nan).astype(np.float32)
```

## 元数据更新

输出文件属性中更新缺测处理说明，明确：

- ocean grid cells set to NaN via `global-land-mask`
- source-missing BCSD grid-time cells set to NaN before writing CF

## 验证方案

以 `NESM3 / ssp126 / Denmark / 2015` 为测试对象，分别重算光伏和风电到临时输出目录：

```bash
cd ref_code/calculate_bcsd_cfs
./.venv/bin/python S01E01_Simulate_Solar_CF_BCSD.py \
  --data_dir data/bcsd_outputs \
  --model NESM3 \
  --region Denmark \
  --scenario ssp126 \
  --years 2015 \
  --output_dir output/test_missing_nan_fix/CFs_of_solar \
  --overwrite

./.venv/bin/python S01E02_Simulate_Wind_CF_BCSD.py \
  --data_dir data/bcsd_outputs \
  --model NESM3 \
  --region Denmark \
  --scenario ssp126 \
  --years 2015 \
  --output_dir output/test_missing_nan_fix/CFs_of_wind \
  --overwrite
```

测试检查：

1. 输出文件存在且可被 xarray 打开。
2. 新输出中源气象缺测位置应为 `NaN`。
3. 新输出中 `NaN` 比例应接近源气象缺测比例与海洋掩膜共同作用后的比例。
4. 与旧输出对比，`skipna=True` 的均值应高于或等于旧输出均值。

## 全量运行

验证通过后，在 `ref_code/calculate_bcsd_cfs/.venv` 环境下后台启动：

```bash
cd ref_code/calculate_bcsd_cfs
nohup bash -lc 'source .venv/bin/activate && bash run_all_solar_NESM3.sh' > logs/run_all_solar_NESM3.nohup.log 2>&1 &
nohup bash -lc 'source .venv/bin/activate && bash run_all_wind_NESM3.sh' > logs/run_all_wind_NESM3.nohup.log 2>&1 &
```

记录两个进程 PID，并用日志文件跟踪运行状态。
