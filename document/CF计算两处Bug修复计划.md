# 容量因子计算脚本：两处 Bug 修复计划

> 起因：运行
> `python S01E01_Simulate_Solar_CF_BCSD.py --data_dir data/bcsd_outputs --model CANESM5 --region all --scenario ssp126 --years 2015-2060`
> 时报错并遗留失败产物。

本计划修复两个独立的 bug，并对所有 `S*E01 / S*E02` 计算脚本做横向排查。

---

## Bug 1：cftime 时间坐标转 float 报 TypeError

### 现象
```
File ".../S01E01_Simulate_Solar_CF_BCSD.py", line 461, in _datetime64_to_ns
    return arr.astype(np.float64)
TypeError: float() argument must be a string or a real number, not 'cftime._cftime.DatetimeNoLeap'
```

### 根因
`compute_region_solar_cf` 把瞬时变量 `tas/uas/vas` 线性插值到 `rsds.time`
（`interp_instantaneous_to_target_times_chunk`），其中调用 `_datetime64_to_ns`
将时间坐标转为数值用于 `searchsorted`。

CANESM5 等 CMIP6 模式使用非标准日历（如 `noleap`），xarray 解码后 time 坐标是
`cftime.DatetimeNoLeap` 对象数组（`dtype=object`），既不是 `datetime64`，
也不能直接 `astype(np.float64)`，于是落到 `return arr.astype(np.float64)` 抛错。

S01E01 当前的 `_datetime64_to_ns` 缺少 cftime 分支：
```python
def _datetime64_to_ns(values):
    arr = np.asarray(values)
    if not np.issubdtype(arr.dtype, np.datetime64):
        return arr.astype(np.float64)          # ← cftime 在这里崩
    return arr.astype("datetime64[ns]").astype(np.int64)
```

### 历史修复
此 bug 早先在 `S04E01_Simulate_Solar_CF_NAM-12.py` 修过
（commit `d962005 fix: cftime 时间差转换使用 total_seconds()`）。
当前 S04E01 的实现为参考标准：
```python
def _datetime64_to_ns(values):
    import cftime
    arr = np.asarray(values)
    if np.issubdtype(arr.dtype, np.datetime64):
        return arr.astype("datetime64[ns]").astype(np.int64)
    if arr.dtype == object and arr.size > 0 and isinstance(arr.flat[0], cftime.datetime):
        base = arr.flat[0]
        return np.array([(v - base).total_seconds() for v in arr.flat], dtype=np.float64)
    return arr.astype(np.float64)
```

### 修复方案
将 S01E01 的 `_datetime64_to_ns` 替换为 S04E01 的实现（增加 cftime 分支，
用 `(v - base).total_seconds()` 得到相对第一个时间点的秒偏移，单调递增，
满足 `searchsorted` 需求）。

---

## Bug 2：计算失败后遗留不完整产物

### 现象
- `output/.../Australia/solar_CF_Australia_..._allmonths.nc` 是上次失败遗留的产物，
  本次因 `out_file.exists()` 被误判为"已存在，跳过"。
- 本次 Austria 在分块写入中途崩溃，又遗留了一个不完整的
  `solar_CF_Austria_..._allmonths.nc`。

### 根因
`compute_region_solar_cf` 先用 `create_output_file` 创建并打开目标 nc 文件，
再在 `try/finally` 分块循环中逐块写入。`finally` 只 `close` 文件和数据集，
**异常时不删除半成品**，于是不完整的 `.nc` 留在最终路径，下次运行被当作已完成跳过。

### 历史修复
`S04E01/S04E02` 的逐年计算函数与合并函数已有此逻辑
（commit `1334593 合并出现异常则删除合并文件`，以及逐年 `except` 删除），
参考 S04E01：
```python
except BaseException:
    if nc.isopen():
        nc.close()
    if yr_out_file.exists():
        logger.error(f"{yr} 计算出错，删除：{yr_out_file}")
        yr_out_file.unlink()
    ...
    raise
```

### 修复方案
把直写脚本的 `try/finally` 改为 `try/except BaseException/finally`：
- 新增 `except BaseException`：先 `if nc.isopen(): nc.close()`，再
  `if out_file.exists(): logger.error(...); out_file.unlink()`，最后 `raise`；
- `finally` 保留，但 `nc.close()` 加 `if nc.isopen()` 守卫避免重复 close 报错，
  仍负责关闭输入数据集。

> 注：ERA5-Land 系列（S02E01/S02E02）采用"写临时文件 `.xxx.nc.tmp` →
> `os.replace` 原子改名"模式，并在开头清理残留 tmp，最终路径只会出现完整文件，
> 因此无此 bug，无需改动。

---

## 横向排查结论（S*E01 / S*E02）

| 脚本 | Bug 1（cftime） | Bug 2（失败产物） | 处理 |
|------|----------------|------------------|------|
| **S01E01** Solar BCSD | ❌ 受影响（本次报错） | ❌ 受影响 | 两处都改 |
| **S01E02** Wind BCSD | ✅ 不涉及（uas/vas 同一时间轴直读，无插值，不调用 `_datetime64_to_ns`） | ❌ 受影响 | 仅改 Bug 2 |
| **S02E01** Solar ERA5Land | ✅ 不涉及（真实日历 datetime64，无 cftime） | ✅ 已用 tmp+替换原子写 | 不改 |
| **S02E02** Wind ERA5Land | ✅ 不涉及 | ✅ 已用 tmp+替换原子写 | 不改 |
| **S03E01** Solar China | ✅ 已有 cftime 分支（`date2num`） | ❌ 受影响 | 仅改 Bug 2 |
| **S03E02** Wind China | ✅ 不涉及（风速直读无插值） | ❌ 受影响 | 仅改 Bug 2 |
| **S04E01** Solar NAM-12 | ✅ 已修 | ✅ 已修 | 不改 |
| **S04E02** Wind NAM-12 | ✅ 已修 | ✅ 已修 | 不改 |

需要改动的文件：
- Bug 1：`S01E01`
- Bug 2：`S01E01`、`S01E02`、`S03E01`、`S03E02`

---

## 执行步骤

1. `S01E01`：替换 `_datetime64_to_ns`（Bug 1）+ 直写循环加 `except` 删半成品（Bug 2）。
2. `S01E02`、`S03E01`、`S03E02`：直写循环加 `except` 删半成品（Bug 2）。
3. 删除已遗留的失败产物（Australia、Austria 等不完整文件）。
4. 用 `python -c "import ast; ast.parse(open(f).read())"` 做语法校验；重跑命令验证。
