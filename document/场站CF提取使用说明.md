# 全网格 CF → 场站 CF 使用说明

## 输入与范围

实现入口是 `prepare_station_cf.py` 和 `extract_station_cf.py`，科学输入仅为已完成的网格 CF 年份 NC。场站 CF 保留源时间、日历、float32 数值和缺测，不乘装机容量，不重新计算天气或风光物理过程。

生产准备固定检查 4 GCM × 3 Climate × 3 Station × 2 tech × 47 patch，共 3,384 个 unit，2015–2060 年八段。`--sample` 仅允许缩小验证范围，发布索引明确标记 sample，不能当正式全量结果。

依赖使用仓库现有 `.venv`，无需新增包。本说明不包含超算部署或调度。

## 1. 准备

复制并填写 [配置示例](station_cf_config.example.json)。`grid_cf_root` 指向含 `inputs/`、`outputs/`、`runtime/completion.json` 的已完成网格运行根。

`station_reference_files` **必填**：指向对应 `Optimization_10km_ssp*` 批次独立导出的 `stations_10km.csv`。准备器按 `(year,type,lon,lat)` 排序，对比坐标集合、记录数和逐位置 `capacity_gw`，确认待使用的远程 CSV 与上游快照一致。参考文件不能直接指向同一个输入 CSV。当前实现接受上游已导出的 CSV，不直接读取优化 NPZ；NPZ 需先通过上游已有 `plot_10km_optimal_stations.py` 导出。

装机语义冻结为 `snapshot_total`。2030/2040/2050 各为累计快照，绝不跨年累加。`first_snapshot_year` 只表示首次出现的快照年。

```bash
.venv/bin/python prepare_station_cf.py \
  --config /path/to/station_cf_config.json \
  --output-root /path/to/cf_stations/run1
```

配置内相对 `grid_cf_root`、`stations_root`、参考 CSV 路径以配置文件所在目录解析；`station_files` 相对 `stations_root` 解析。输出根不能与源 CF 根重叠。

准备器串行生成 catalog、全局最近邻 mapping、源清单和任务表，最后发布 `prepared.json`。过程只读源 CF 的时间/坐标/domain 元数据及 sidecar，不扫描全部 CF 数值，不 hash 大型 CF NC。

最近邻在全局规则格网中选择，再查来源 patch；半格等距时按规范化坐标升序决策。最近格点域外时保留 fill，无 source patch 的位置保留覆盖原因。准备器要求源年份分片的空间 chunk 布局一致，保证站点写出顺序一致；不满足时明确报错。

重复准备只复用身份、文件 hash/stat 和实现均一致的产物。改配置、代码或源输入使用新输出根。已发布 prepared 冻结了代码 SHA 和实际源码 hash，未提交源码也能追溯。

## 2. 提取单个双情景 unit

```bash
.venv/bin/python extract_station_cf.py \
  --prepared /path/to/cf_stations/run1/prepared.json \
  --model CANESM5 \
  --climate-scenario ssp126 --station-scenario ssp585 \
  --tech wind --patch R02C09 \
  --processes 8 --time-chunk 240 --station-chunk 1024 --compress-level 2
```

`--patch` 是映射后 CF 格点的 `source_patch`。两个情景必须显式给出；Station 的 ssp560 别名归一化为 ssp585，Climate 不接受 ssp560。

进程池使用 spawn，一层最多八个进程，各自只读一个年份源文件并独占输出文件。时间/空间按 chunk 流式读取，在内存中配对 gather；默认不合并八段。`processes=1` 可用于基线或调试。

结果目录：

```text
outputs/<model>/climate_<C>/station_<S>/<source_patch>/<tech>/
  manifest.json
  audit.json
  blocks/<unit_identity>/cf_<start>-<end>.nc[.json]
```

输出含 CF、station_id、经纬度、来源行列、domain 状态、距离、first_snapshot_year，以及逐站 `valid_count/missing_count`。全局属性及 sidecar 保存源/映射/目录/实现身份。sidecar 记录读取/写出时间、读取数组字节数、RSS、有效值范围和样本审核方式，可用于后续真实样本性能评估。

同一 unit 同时只能有一个写入者。成功分片复用；孤立 NC 或缺 sidecar 的分片保留原文件，新尝试写入 `attempt_<uuid>/`，manifest 指向新成功版本。更改 time/station chunk 或压缩等级会改变编码身份，需新输出根；更改进程数不会改变科学结果身份。

没有来源分配站点的 unit 发布 `EMPTY_NO_STATIONS`，不写空 NC；有站但全域外的 unit 写完整 fill 文件。无来源站点不丢弃，由全局 reader 补缺测。

## 3. 审核并发布全局索引

全部任务处理后，使用 8 个 spawn 进程并行验收，主进程统一检查覆盖并发布索引：

```bash
.venv/bin/python extract_station_cf.py \
  --prepared /path/to/cf_stations/run1/prepared.json --publish --processes 8
```

检查所有 unit/八段、NC/sidecar 身份、时间、站点、stat、固定种子源值抽样、catalog 与覆盖守恒。每个流块在提取时已检查有效范围；独立审核为抽样，不是对全部源数组二次扫描。

发布 `index/authoritative_index.json`、`coverage_summary.csv`、`unmapped_stations.csv.gz`、`validation_summary.json`。未完成/失败任务会阻止发布。覆盖表同时给出 MATCHED/OUTSIDE_DOMAIN/TOO_FAR/NO_SOURCE_PATCH 和有效/缺测样本总数；处理完成不意味着所有位置都有有效 CF。

## 4. 惰性读取与容量快照

```python
from station_cf_reader import open_station_cf, iter_station_cf
from station_cf_catalog import capacity_snapshot
from grid_cf_io import read_json

index = '/path/to/cf_stations/run1/index/authoritative_index.json'
ds = open_station_cf(
    index, model='CANESM5', climate_scenario='ssp126',
    station_scenario='ssp585', tech='wind', years=(2030, 2039),
)
# CF 数组惰性读取；.values/.load() 才读取指定部分。
small = ds.isel(time=slice(0, 240), station=slice(0, 100)).load()

for block in iter_station_cf(
    index, model='CANESM5', climate_scenario='ssp126',
    station_scenario='ssp585', tech='wind',
    time_chunk=240, station_chunk=1024,
):
    # 每块有界，CF 缺测为 NaN；按实际分析定义处理有效分母。
    pass

prepared = read_json('/path/to/cf_stations/run1/prepared.json')
capacity_gw = capacity_snapshot(
    prepared['catalogs']['ssp585'], 'wind', 2040,
    station_ids=small.station_id.values.tolist(),
)
```

`open_station_cf` 返回 xarray Dataset，通过本地 NetCDF 按需后端读取，无需 Dask；索引、站点/时间坐标会载入内存，全球 CF 数组不会自动载入。不要对全球 Dataset 无条件调用 `.load()`。选择少量站点可传 `station_ids=[...]`，返回顺序与请求一致；未知 ID 报错，已知未覆盖 ID 返回 NaN 和空间状态。

时间坐标使用原生日历的 cftime；不同模式/技术不自动对齐。不完整快照年份（例如 2035）不能由 `capacity_snapshot` 读取，也不默认插值。2030 快照没有某位置且整个快照已核验时，其装机返回 0；CF 仍保存该位置全部源气候时间。

## 5. 本地小样本与验证范围

可通过测试夹具创建完全独立于 BCSD 的小型源 CF、场站表和参考表：

```bash
PYTHONPATH=tests .venv/bin/python - <<'PY'
from pathlib import Path
from station_cf_fixtures import build
from prepare_station_cf import prepare
from extract_station_cf import extract_unit, publish
from grid_cf_io import read_json
root = Path('tmp/station_cf_example').resolve()
config = build(root, climates=('ssp126', 'ssp245'), techs=('wind', 'solar'))
prepared = prepare(config, root / 'station_output', sample=True)
for task in read_json(prepared)['tasks']:
    extract_unit(prepared, model=task['model'],
                 climate_scenario=task['climate_scenario'], station_scenario=task['station_scenario'],
                 tech=task['tech'], patch=task['source_patch'], processes=1)
print(publish(prepared))
PY
```

夹具创建命令用于新的空目录；重复验证可读取已有 prepared，按提取器恢复规则执行。

```bash
.venv/bin/python -m pytest \
  tests/test_station_cf_catalog.py tests/test_station_cf_mapping.py \
  tests/test_station_cf_extract.py tests/test_station_cf_reader.py \
  tests/test_grid_cf.py tests/test_station_gather_crop.py -q
```

本轮完成代码及合成数据验证：新增 23 项测试和已有 46 项相关回归均通过；语法编译和 diff 空白检查通过。文中双 Climate、双 Station、双技术、双 patch 示例已实际完成 16 个 unit、32 份年份 NC 的发布和非对角情景读取，产物位于本地忽略目录 `tmp/station_cf_example/`。这些是测试数据，不是全球场站生产结果。

真实远程 CSV 与优化导出批次比对、真实网格 CF 的密集/稀疏 patch 性能测试和全量场站生产尚未执行；不能用合成样本耗时承诺真实吞吐。本轮不创建超算运行目录或调度脚本。

## prepare 并行与时间轴复用

`prepare_station_cf.py --processes 16` 使用 spawn 进程池并行扫描源 CF 的 patch；允许1–16个进程。Python API 的 `prepare(..., processes=1)` 保留串行默认值。SCNet 生成器默认 `--prepare-processes 16 --prepare-cpus 16`，且不能超过 `--prepare-cpus`。

每份源 NC 仍读取并校验实际时间轴、坐标、掩膜、元数据和 sidecar。时间轴按完整内容（含 dtype/shape）、units、calendar 作为缓存键；只有这些全部一致才复用已验证的日历转换结果。目录和最近邻 mapping 的算法、任务集合、完成判据不变。日志按 source_index/catalogs/mappings 输出完成数量、已用时间和当前阶段剩余时间估计；阶段估计不包含后续阶段。

修改运行中进程尚未加载的代码或共享输出并不能安全加速已有任务。并行版本应使用独立固定 SHA checkout 和新输出根；旧作业可继续完成，两者的 prepared 身份分别保留。代码版本变更后，后续提取账号也必须部署与所选 prepared 一致的 SHA。
