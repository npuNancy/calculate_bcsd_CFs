# 从全网格 CF 提取全球场站 CF：Climate × Station 全组合实施方案

编制日期：2026-10-01。状态：核心代码与本地合成数据验证已完成；真实样本性能验证及场站 CF 生产尚未执行。实际接口见 [使用说明](场站CF提取使用说明.md)。

## 1. 目标与范围

从乌镇 1872 已完成的全球网格容量因子（capacity factor，下文 CF）中，以最近邻方式提取全球场站位置的风电、光伏 CF，完整保留 2015–2060 年的源时间序列。气候情景与场站装机情景作为两个独立维度，支持对应情景和交叉情景的反事实分析。

本阶段的数据链路为：

```text
已完成的网格 CF + 三套场站位置/装机表
                  ↓
冻结 CF 文件索引与场站目录
                  ↓
全局最近邻映射 → 按来源 patch 分组
                  ↓
逐年份段多进程提取 CF
                  ↓
场站 CF 分片 + 场站目录/容量记录 + 映射 + 全局索引/覆盖报告
```

本方案实施范围为科学数据处理代码、文件契约、本地验证和性能测试；不创建或完善 `infos/scnet_patchify_stations/`，不涉及账号分工、Slurm、提交监控或部署。输入从已计算好的网格 CF 开始，不重新执行 BCSD，也不从气象变量重新计算 CF。

主要设计选择：

| 问题 | 本方案选择 |
|---|---|
| 情景维度 | `climate_scenario`、`station_scenario` 独立且均为必填 |
| 处理单元 | GCM × Climate × Station × tech × source_patch，共 3,384 个逻辑 unit |
| 空间方法 | 规则经纬网轴最近邻 `nearest_regular`，先全局选格点，再确定来源 patch |
| CF 数值 | 直接复制源 `wind_cf` / `solar_cf`，保留 float32 和缺测 |
| 时间 | 原样保留源 CF 的时间、单位、日历及八个年份分片 |
| 并行 | 一层 `spawn` 进程池，默认最多 8 个进程，各处理独立年份文件 |
| 保存 | NetCDF4 `(time,station)` 分片，显式 Climate/Station 目录层级 |
| 容量 | 独立目录保存，不乘入 CF，不按投运时间将 CF 置零 |
| 全球访问 | 权威 JSON 索引 + 按需读取接口；不要求合并成单个全球 NC |
| 完成 | 3,384 个 unit 均有明确终态，且完整场站目录中的每个位置都有覆盖解释 |

## 2. 当前输入与已有代码依据

### 2.1 已完成网格 CF

```text
源运行根：/work/home/acjpoxgsdu/cf_grid/cf_grid_v2_20261001_013720
CF 目录： <源运行根>/outputs
完成记录：<源运行根>/runtime/completion.json
冻结配置：<源运行根>/inputs/release.json
patch 定义：<源运行根>/inputs/patch_manifest.json
```

本地完成记录表明：1,128 个网格 unit、9,024 份年份 NC，最终验收时间为 2026-10-01 03:22:54（Asia/Shanghai），年份 NC 大小合计 9,117,147,416,414 bytes；没有要求合并 final。输入源科学代码 SHA 为 `df01ef4a9dccfce7758541cec0f652238c42bb52`。

现有文件布局由 `patchify_grid_cf.py` 定义：

```text
outputs/<model>/<climate_scenario>/<patch>/<tech>/
  manifest.json
  blocks/<grid_unit_identity>/
    cf_2015-2020.nc
    cf_2015-2020.nc.json
    ...
    cf_2056-2060.nc
    cf_2056-2060.nc.json
    manifest.json
```

读取入口以技术目录的 `manifest.json` 为准，其 `blocks[]` 记录实际文件路径、年份、序号和状态；不得通过查找最新 identity 目录选择输入。`provenance.scenario` 和 NC 的 `scenario` 在新流程中均解释为 **Climate**。

当前八段为：2015–2020、2021–2026、2027–2032、2033–2038、2039–2044、2045–2050、2051–2055、2056–2060。实现应读取 manifest 中的实际合同，再核对这八段覆盖 2015–2060，不自行重新切分。

源 NC 主要字段：

| 字段 | 维度 / 类型 | 用途 |
|---|---|---|
| `wind_cf` 或 `solar_cf` | `(time,lat,lon)` / float32 | 唯一科学数值输入，单位 `1` |
| `time` | `(time)` / 源数值类型 | 保留数值、`units`、`calendar` |
| `lat`, `lon` | 一维坐标轴 | 最近邻几何与源行列定位 |
| `domain_mask` | `(lat,lon)` / int8 | 1 为计算域，0 为域外 |
| 全局属性 | model、scenario、patch_id、tech、identity、grid_fingerprint 等 | 身份校验 |

源缺测值由 `grid_cf_io.FILL = np.float32(9.96921e36)` 定义；CF 有效范围为 `[0,1]`。当前源写出配置为时间 chunk 240、空间 tile 64×64、压缩等级 2，读取时仍从每个文件读取实际 chunk 信息。

上游网格 CF 的 `input_mode=blocks/final` 只作为来源 provenance 保留。本提取器统一读取上述 **CF 年份 NC**，不根据这一字段重新选择 BCSD blocks/final。

### 2.2 场站 CSV

用户指定目录：

```text
/work/home/acbw9wpn5k/project_climate_patchify/shared_run_20260913/inputs/stations
```

按参考项目约定显式配置：

| 规范 Station 名称 | 场站 CSV | 原始数据标签 |
|---|---|---|
| `ssp126` | `stations_SSP1-2.6.csv` | SSP1-2.6 |
| `ssp245` | `stations_SSP2-4.5.csv` | SSP2-4.5 |
| `ssp585` | `stations_SSP5-6.0.csv` | SSP5-6.0 / 用户所称 ssp560 |

`ssp560` 只允许作为 **Station 输入别名**，解析后统一为 `station_scenario=ssp585`，同时保存 `station_source_label=SSP5-6.0` 和原 CSV 身份。它不构成第四种情景，也不表示存在 Climate ssp560 数据。输出不同时创建 `ssp560` 与 `ssp585` 两份任务。

字段依据参考项目的实际 catalog 实现：`year,type,lon,lat,capacity_gw`。本次编写方案没有重新全量读取远程 CSV；正式实现的准备阶段需要核验文件名、表头、内容 hash、字段范围及年度含义，不能把参考文档中的数量直接当成本批实际数量。

### 2.3 可复用与需要新增的部分

| 现有内容 | 处理方式 |
|---|---|
| `patchify_grid_cf.py` 的 manifest、年份合同和 `spawn` 组织方式 | 用作 CF 输入和并行组织依据，不调用它的气象读取/计算入口 |
| `grid_cf_io.py` 的 digest、file_identity、原子 JSON、输出锁 | 可直接复用这些小工具；新写站点 schema 校验与 writer |
| `patchify_station_cf.py` | 保持现有用途；其 loader 丢弃重复位置时只保留一条容量，不适合作为新 catalog 的完整实现 |
| `cf_physics.py` | 新提取链路不调用物理计算函数 |
| 参考项目的 `station_catalog.py`、`station_match.py` | 对齐 ID、经度规范化、年度容量记录及覆盖语义；实现移植到本仓库，不依赖另一个仓库的绝对路径 import |

参考方案中的事件变量、int8 三态和十个时间段属于极端事件产品。本 CF 产品使用 float32 CF 和现有八段；不直接套用那些 schema 或文件数。

## 3. Climate × Station 组合及反事实语义

### 3.1 笛卡尔积

```text
models  = CANESM5, MPI-ESM1-2-HR, MRI-ESM2-0, BCC-CSM2-MR
Climate = ssp126, ssp245, ssp585
Station = ssp126, ssp245, ssp585
tech    = wind, solar
patch   = 冻结 patch_manifest 中的 47 个 patch
```

每个模型、技术、patch 均处理以下九种组合：

| Climate ＼ Station | ssp126 | ssp245 | ssp585（源 SSP5-6.0） |
|---|---|---|---|
| ssp126 | 对应情景 | 交叉情景 | 交叉情景 |
| ssp245 | 交叉情景 | 对应情景 | 交叉情景 |
| ssp585 | 交叉情景 | 交叉情景 | 对应情景 |

合计：`4 × 3 × 3 × 2 × 47 = 3,384` 个逻辑 unit，其中对应情景 1,128 个、交叉情景 2,256 个。若所有 unit 都有来源场站，则有 `3,384 × 8 = 27,072` 份输出年份 NC；实际文件数需扣除没有场站的 unit，逻辑 unit 数保持 3,384。

`patch` 在这里明确为 **最近网格的来源 patch（source_patch）**。它仍取同一套 47 个 ID，但不等同于按场站坐标先做 bbox 分类的 `station_patch`。

### 3.2 计算定义

设气候情景为 C，场站情景为 S，技术为 T：

```text
station_cf[M,C,S,T,t,s] = grid_cf[M,C,T,t,iy(s),ix(s)]
s 属于 catalog[S,T]
```

Station 决定位置集合、站点身份及下游装机表；Climate 决定取哪一套源 CF。例：`Climate=ssp126, Station=ssp585` 表示在 SSP5-6.0 场站布局上抽取 ssp126 气候对应的 CF。

同一模型、Climate、技术、时间及同一源格点的 CF 不随 Station 改变；不同 Station 的站点数量和布局可能不同。CF 是单位装机的无量纲系数，抽取时不乘容量、不做容量加权，不改变网格 CF 使用的风机/光伏参数。光伏结果也是已经按源格点位置计算好的 CF，不在场站经纬度重新运行太阳几何。

任务 ID 建议：

```text
station-cf-v1/<model>/climate-<C>/station-<S>/<tech>/<source_patch>
```

所有任务记录、manifest、输出路径、索引筛选、日志摘要和读取 API 必须同时携带两个情景字段；不提供含义不明的单一 `scenario` 参数，也不在缺省时暗中令 Station 等于 Climate。

## 4. 场站目录、ID 与容量记录

### 4.1 一次准备，跨 Climate/GCM 复用

按 Station × tech 构建六套位置 catalog。位置粒度沿用参考项目：同一 Station/tech 下跨年份同址记录归为一个“场站位置”，并非保证对应一个具有实体电厂编号的独立电厂。

处理步骤：

1. 校验必需列、数值有限性、纬度范围、整数年份、wind/solar 类型和非负容量；保留 `source_row` 以便回溯。未知类型、非数值容量等输出问题清单并阻止该目录发布，不静默填 0 或删除。
2. 经度规范化到 `[-180,180)`，明确原始经度允许区间，拒绝明显非法经度后再归一化；坐标使用 float64 保存。
3. 检查 `(year,type,lon,lat)` 重复键；重复默认报错，不能自行双计容量。
4. 按规范化后的原始位置精确去重，保留全部年度容量行，生成稳定 ID。
5. 各技术 catalog 按 `station_id` 稳定排序；写入原始 CSV SHA256、实现版本、验证统计和 catalog hash。

沿用参考项目身份格式，原来的 scenario 在这里明确解释为 Station：

```python
station_id = sha1(
    f"{station_scenario}|{tech}|{lon:.4f}|{lat:.4f}".encode()
).hexdigest()[:20]
```

ID 不含 GCM、Climate、patch、年份。固定 Station/tech 的同一位置在所有 Climate 和 GCM 中保持相同 ID。不同 Station 即便同坐标，ID 也不同；跨 Station 对比用技术和规范化坐标建立显式位置对照表，不能误用 station_id 等值连接。不同原始坐标如四位格式化后发生碰撞，应报告并停止，不能据此合并。

### 4.2 年度装机快照

**容量语义确定为 `capacity_semantics=snapshot_total`：每行 `capacity_gw` 表示该位置在 `year` 对应规划年份的累计装机存量，不是该年的新增建设量。** 年份为 2030、2040、2050，catalog 汇集三个快照的位置并集；容量表保留每个快照的原始记录。

保存独立的 `capacity_rows.csv.gz`，至少包括：`station_scenario,tech,station_id,year,capacity_gw,source_row`。其中 `year` 是 `snapshot_year` 的原始列名，元数据明确其含义。CF 文件不保存一个会被误解为全期有效的固定 `capacity_mw`。

读取规则：

```text
capacity(S,T,s,Y) = capacity_gw[Station=S, tech=T, station_id=s, year=Y]
```

同一快照内每个键必须唯一；2030、2040、2050 的值不能跨年求和。例如 2030 年 1 GW、2040 年 1.5 GW，则 2040 年装机为 1.5 GW，该十年净增为 0.5 GW。需要增量时才对两个完整快照按位置对齐后做差；首个 2030 快照无法推断 2030 年当年新增量。

`first_snapshot_year=min(year)` 仅表示位置首次出现在规划快照中的年份，不表示实际投运年份。CF/catalog 可保存该字段，不能将其命名为 `activation_year` 或用于推断 2015–2029 的装机为零。位置未出现在某个**完整且核验通过的快照**时，依照上游仅导出正容量位置的契约解释为该快照未配置容量；整个快照/来源文件缺失则是输入错误，不能当作零装机。CSV 容量四位小数舍入后可能为零，不能以 CSV 值为零删除其已有位置记录。

**主 CF 保留全部源时间，不按快照年份置零或删除。** 下游以独立 `capacity_snapshot_year` 选择 2030、2040 或 2050 布局，可以将该快照用于明确指定的气候时间窗口；该参数不增加本次 CF 的 3,384 个 unit，也不要求把同位置 CF 按三个快照重复保存。

现有三快照不能唯一确定 2031–2039、2041–2049 或范围外各年的实际装机。下游若需逐年装机，必须显式选择并记录插值/阶梯保持/外推规则；本接口默认仅接受已有快照年份，其他年份报错，不隐式累计或前向填充。

依据来自三个 `Optimization_10km_ssp*` 的保存逻辑和 CSV 导出代码，详见第 13 节。用户指定远程 CSV 的具体导出批次尚需在输入冻结阶段与源快照比对；该校验用于确认文件来源，没有理由继续沿用参考事件项目的建设年份累计解释。

## 5. 全局最近邻映射

### 5.1 方法的精确定义

采用与参考方案相同的规则经纬轴最近邻，方法名为 `nearest_regular`：

```text
iy = argmin(abs(grid_lat - station_lat))
ix = argmin(abs(((grid_lon - station_lon + 180) % 360) - 180))
match_dist_deg = max(abs(delta_lat), abs(periodic_delta_lon))
```

这是规则经纬坐标的最近邻定义，`match_dist_deg` 是角度 L∞ 距离；不称其为球面最短距离最近邻。可另记 haversine `match_dist_km` 作为空间代表性诊断，但不以该值改变默认映射。

默认最大角距离 `0.15°`，沿用参考项目兼容门槛，并写入 mapping identity。实际网格分辨率、原点和轴方向必须从源 CF 坐标验证。正常 0.1° 完整格网最近轴差约不超过 0.05°，因此额外报告精确命中（≤1e-6°）、≤0.05°、0.05–0.15°和超门槛的数量，调查异常边界匹配。

等距容差取 `1e-8°`；候选按规范化坐标 `(lat,lon)` 升序选定，重复坐标按经过验证的 core owner 归属处理。保留原文件的实际零基行列索引，兼容坐标升序和降序。

实现采用排序轴 `searchsorted`、邻点候选和经度首尾周期候选，分批处理场站；不构造“全球轴长 × 全部场站数”的距离矩阵。首尾、纬度极点、日期变更线及恰好半格位置须有确定性测试。

### 5.2 先选格点，再决定 patch

只从已验证的 CF 文件读取一维经纬轴和二维 domain mask，结合冻结 patch core bbox 建立全局网格元数据。网格几何不需要重新读取上游气象或 land plan。

1. 对每个空间组核对规则轴的间距、原点、方向和 core 边界，建立理论规则网格与坐标到 `source_patch/local_iy/local_ix` 的归属表。
2. 理论几何由已核验的规则轴和 patch 定义确定；稀疏覆盖区的空缺仍保留。如果无法确定格网原点/分辨率，停止该组准备，不猜测。
3. 在完整理论几何上选择最近格点，再查它是否有本轮源 patch。缺失位置记为 `NO_SOURCE_PATCH`，不吸附到空白区域另一侧的可用格点。
4. 选中的格点可以位于场站 bbox 的相邻 patch；根据 `source_patch` 安排提取。另存几何归属 `station_patch`，供质量报告使用，不能用它限制候选搜索。
5. core bbox 西闭东开、南闭北开，北极端点特殊纳入；跨日期变更线按周期处理。重复 core 坐标或不一致的网格归属属于输入冲突，停止发布映射。

每个合法 catalog 位置在每个空间组只有一条 mapping 记录和一个空间状态。没有 source_patch 的位置保留在全局 mapping 中，不为它虚构第 48 个 patch。

### 5.3 domain 与缺测

最近邻先按几何选择，不寻找“最近有效陆地格点”。选中格点固定后，不随时间缺测更换位置。

| mapping_status | 含义 | 保存方式 |
|---|---|---|
| `MATCHED` | 有来源、距离合格、domain=1 | 进入来源 patch NC，逐值复制 CF；局部缺测保留 |
| `OUTSIDE_DOMAIN` | 最近来源格点 domain=0 | 保留在来源 patch NC，整列 fill |
| `TOO_FAR` | 最近几何格点超过门槛 | 全局 mapping 登记；reader 补 fill，不计入可提取站点 |
| `NO_SOURCE_PATCH` | 最近格点无来源 patch | 全局 mapping 登记；reader 补 fill |

非法 CSV 行不进入合法 catalog，写单独错误报告并在准备阶段解决；不得以 `NO_SOURCE_PATCH` 掩盖数据错误。

CF 的 0 是有效值（例如夜间光伏或无发电风况），不能等同缺测。domain=1 但某些时刻源 CF 缺测时，仅相应时刻 fill。所有场站均 fill 的 unit 与没有场站的 unit 区分处理。

### 5.4 缓存与复用

mapping 的键包含：catalog hash、Station、tech、源坐标指纹、mask 指纹、patch manifest hash、最近邻方法/门槛/等距规则、实现版本。空间一致的 GCM/Climate 可共用 mapping，但输入 CF 文件身份仍按 GCM/Climate 单独保留。

不能仅因名义分辨率相同就认为四模型、三 Climate 几何相同。准备阶段对所有源分片核对坐标/mask 一致性；这些是空间元数据校验，不扫描全部 CF 时间数组。指纹一致后才合并空间组。mask 不同则独立记录覆盖状态。

同一格点映射到多个站点时，只读取一次唯一 `(iy,ix)`，再广播到各 `station_id`。不同 Station 的共享位置可复用几何计算，首版仍按独立 unit 保存 CF，不增加跨 unit 的共享写入机制。

## 6. 代码新增与接口设计

保持当前仓库顶层脚本/模块风格，已新增六个职责明确的文件：

| 文件 | 职责 / 主要接口 |
|---|---|
| `station_cf_catalog.py` | 情景规范化、CSV 校验、稳定 ID、目录/原始容量表；`build_catalog(...)` |
| `station_cf_mapping.py` | 空间组、全局最近邻、边界归属、domain 状态与覆盖报告；`build_mapping(...)` |
| `station_cf_io.py` | 只读源 CF manifest/sidecar、输入计划、分块 gather、站点 NC writer、身份/恢复检查 |
| `prepare_station_cf.py` | 顺序构建冻结输入索引、六套 catalog、mapping 和 3,384 个任务描述；发布 `prepared.json` |
| `extract_station_cf.py` | 单 unit 入口、八年份任务的 spawn 并行、unit 锁、manifest、提取后审核；`extract_unit(...)` |
| `station_cf_reader.py` | 双情景索引查询、站点/时间子集读取与缺覆盖补 fill；`open_station_cf(...)`、`iter_station_cf(...)` |

相应新增 `tests/test_station_cf_catalog.py`、`tests/test_station_cf_mapping.py`、`tests/test_station_cf_extract.py`、`tests/test_station_cf_reader.py`。只有实际出现共用逻辑时才提取通用工具，不先引入额外框架。

准备入口的配置使用 JSON，示意核心字段：

```json
{
  "grid_cf_root": "/work/home/acjpoxgsdu/cf_grid/cf_grid_v2_20261001_013720",
  "stations_root": "/work/home/acbw9wpn5k/project_climate_patchify/shared_run_20260913/inputs/stations",
  "station_files": {
    "ssp126": "stations_SSP1-2.6.csv",
    "ssp245": "stations_SSP2-4.5.csv",
    "ssp585": "stations_SSP5-6.0.csv"
  },
  "station_reference_files": {
    "ssp126": "/path/to/verified-upstream-export/ssp126/stations_10km.csv",
    "ssp245": "/path/to/verified-upstream-export/ssp245/stations_10km.csv",
    "ssp585": "/path/to/verified-upstream-export/ssp560/stations_10km.csv"
  },
  "models": ["CANESM5", "MPI-ESM1-2-HR", "MRI-ESM2-0", "BCC-CSM2-MR"],
  "climate_scenarios": ["ssp126", "ssp245", "ssp585"],
  "station_scenarios": ["ssp126", "ssp245", "ssp585"],
  "techs": ["wind", "solar"],
  "years": "2015-2060",
  "spatial_method": "nearest_regular",
  "max_distance_deg": 0.15
}
```

`station_reference_files` 必填，指向对应优化批次独立导出的 CSV；准备器逐记录核对快照，不把输入 CSV 自身当来源证明。当前接口比较上游导出的 CSV，NPZ 通过上游已有工具先导出。

patch 列表从冻结源 manifest 读取，生产必须恰为 47 个；仅合成验证可显式使用 `--sample` 缩小范围。`prepared.json` 指向已发布的配置、输入索引、catalog、mapping 和任务表及其 hash。准备过程串行发布共享资产，提取子进程全部只读这些资产。

已实现 CLI（真实数据需先完成路径与快照来源配置）：

```bash
.venv/bin/python prepare_station_cf.py \
  --config <station_cf_config.json> --output-root <station_cf_run_root>

.venv/bin/python extract_station_cf.py \
  --prepared <station_cf_run_root>/prepared.json \
  --model CANESM5 \
  --climate-scenario ssp126 --station-scenario ssp585 \
  --tech wind --patch R02C09 \
  --processes 8 --time-chunk 240 --station-chunk 1024 \
  --compress-level 2
```

提取输出根从 `prepared.json` 读取，避免命令行另传路径导致同一 prepared 身份写到不同集合。`--processes` 为正整数，最大实际 worker 数为 `min(processes,未完成年份段数,8)`；`--processes 1` 提供串行测试/性能基线。

首版使用已有 NumPy、pandas、netCDF4、xarray，以及已有环境中的 cftime 和 Python 标准库，不要求新增 scipy、Dask、pyarrow 或 Zarr。接口不以安装风机库作为提取 CF 的科学前置条件。

## 7. 输入准备和身份冻结

准备阶段按顺序完成：

1. 校验源完成记录，建立 1,128 个源 CF unit 的精确 manifest 列表。核对四模型、三 Climate、两技术、47 patch 无遗漏/重复。
2. 对每个源 unit 检查 `COMPLETED`、八份 NC/sidecar、block identity、路径/stat、年份合同及源变量 schema；冻结源 manifest/sidecar hash、CF 路径/size/mtime_ns 和源 run/release 身份。
3. 读取必要的时间坐标、空间轴和 mask，建立空间组/时间合同。这个阶段不读取完整 CF 数值，也不对 TB 级 NC 做全文件 hash。
4. 构建三 Station 的六套 catalog 和原始容量表，校验原 CSV 内容 hash。
5. 建立 mapping，报告所有场站去向。按 source_patch 分组，写 3,384 条任务描述，标出每个 unit 的站点数、唯一格点数和源八段路径。
6. 验证逻辑任务数、九种情景配对、所有来源和 catalog/mapping hash 后，最后原子发布 `prepared.json`。

源清单的 `source_scenario` 适配到 `climate_scenario`，不能将 Station 传给源路径定位函数。提取开始和结束重新核对当前 unit 的输入身份/stat；发现已冻结文件变化则报错，不自动切换到另一个源 identity。

身份分层：

- `catalog_identity`：CSV 内容、Station 命名及 ID/容量语义版本。
- `mapping_identity`：catalog、网格/掩膜/patch 指纹及匹配参数。
- `unit_identity`：model、Climate、Station、tech、source_patch、源 CF 身份、mapping/catalog 身份、年份、提取实现及输出编码配置。
- `block_identity`：unit identity + 源 block identity + 年份段/时间指纹。

并行进程数不改变科学身份；chunk/压缩作为输出编码配置冻结。路径、元数据和身份均包括 Climate/Station，避免交叉情景互相覆盖。

## 8. 多进程与流式提取算法

### 8.1 一层年份进程池

每个逻辑 unit 内以现有八个年份 CF 文件为八个任务：

```text
父进程：输入/映射校验，输出锁，列出未完成分片
    ├─ worker：2015–2020 CF → station CF 文件
    ├─ worker：2021–2026 CF → station CF 文件
    ├─ ...
    └─ worker：2056–2060 CF → station CF 文件
父进程：收集返回摘要，逐次更新 manifest，审核八段后完成
```

使用 `ProcessPoolExecutor(..., mp_context=multiprocessing.get_context("spawn"))`。父进程不持有跨 spawn 传递的 NetCDF/HDF5 句柄；只传小型参数及 mapping 路径，worker 内独立打开只读源和自己的输出。进程间不传 CF 大数组，不嵌套进程池。

每个进程独占一个输出 NC；父进程独占 unit manifest 写入。限制底层数值库线程为 1，避免进程与线程相乘。默认 8 是年份并行度上限，不是性能承诺；通过 1/2/4/8 进程对比决定实际推荐值。

一个年份失败时，父进程记录原因，未启动任务可取消，已完成分片保留；所有已运行 worker 正常关闭文件后结束。下一次运行只处理未完成且身份一致的分片。

### 8.2 单年份任务

1. 打开源 NC 一次，核验变量、维度、时间、grid identity、domain mask 和 mapping 索引范围。
2. 根据源 `chunking()` 将唯一 `(iy,ix)` 按空间 chunk 分组。对源连续布局也提供有界 tile 读取路径。
3. 将站点按 `(source_chunk_y,source_chunk_x,iy,ix,station_id)` 排序，以便连续写入，保存这一固定输出顺序。catalog 的全局 ID 顺序可以不同；跨文件连接一律用 ID。
4. 沿时间读取有限 chunk；对包含选中格点的源空间 chunk 读取连续切片 `var[t0:t1,y0:y1,x0:x1]`，在 NumPy 内存数组上用配对索引 `tile[:,local_y,local_x]` gather。
5. 不能对磁盘 netCDF4 变量直接同时传两个离散索引数组，避免正交索引产生笛卡尔积。一个源格点的多个站点从已读结果广播。
6. 按连续 station block 写输出，使用有界缓冲；domain 外列直接填缺测，不必读取其 CF tile。没有任何被选中格点的空间 chunk 不读取。
7. 同步统计 valid/missing 数和有效最小最大值，记录源读取、gather、写出、总耗时与进程峰值 RSS。
8. 关闭文件，审核 header/统计及固定种子抽样值，原子发布 NC/sidecar；向父进程返回小型完成摘要。

对于单个源 chunk 映射到大量场站的情况，在内存中分 station block 广播，不一次创建该 chunk 对应的所有站点时间矩阵。输出缓存大小显式有界，不能先构造 `(46年全部时间,全球全部场站)`。

### 8.3 数值与时间保持

源 CF 是无 scale/offset 的 float32。推荐关闭自动 mask/scale 后读取 raw float32，显式识别 `_FillValue`、`missing_value` 和非有限值：正常有效值逐位复制，缺测统一编码为输出 fill；若发现本批预期外的 scale/offset 或 dtype，则报 schema 不一致，不能默默改变精度。

有效范围检查不通过时报告异常，不通过 clip 将错误值压回 `[0,1]`。源 fill 不能进入统计，更不能变成 0。正常来源的有限值应在无损压缩写出后与源选中值精确相等。

每个分片原样复制源 `time` 数值、dtype、units、calendar。禁止通过 pandas 强制转为公历，不补齐源本来不存在的时刻，不重采样或修复不同技术的时间偏移。八段间按日历解码后检查时间单调、无重叠和源一致性；不能仅用文件名判断连续。

### 8.4 内存和 I/O 预算

单进程主要数值缓冲近似为：

```text
4 × T_chunk × Y_tile × X_tile
+ 4 × T_chunk × S_buffer
+ 缺测掩膜/索引 + HDF5 cache + 元数据
```

以 `T=240,Y=X=64,S=4096` 为例，两项 float32 数值缓冲分别约 3.75 MiB；实际 RSS 还包含库缓存、mapping 与临时数组，应测量后乘以进程数并加父进程开销。避免无意生成 float64 全量副本，HDF5 chunk cache 设置明确上限。

三套 Station 会对同一个 Climate 源进行最多三次逻辑读取，这是首版独立 unit 的可解释代价。先用唯一格点去重和 chunk 分组减少重复解压；若实测确为瓶颈，后续可以在相同 model/Climate/tech/patch/年份下对三套 Station 合并读取再分发，但仍保留独立 3,384 个逻辑身份。该优化不作为首版完成条件。

## 9. 结果保存、分片契约和读取

### 9.1 输出目录

建议在 1872 使用与源 CF 独立的新运行根，例如以下占位路径；具体目录名在实现使用时确定：

```text
/work/home/acjpoxgsdu/cf_stations/<RUN_ID>/
  prepared.json
  inputs/
    input_release.json
    grid_cf_index.json
    station_cf_config.json
  catalogs/<station_scenario>/<catalog_identity>/
    wind/stations.csv.gz
    solar/stations.csv.gz
    capacity_rows.csv.gz
    catalog.json
  mappings/<mapping_identity>/
    mapping.nc
    coverage.json
  outputs/<model>/climate_<C>/station_<S>/<source_patch>/<tech>/
    manifest.json
    blocks/<unit_identity>/
      cf_2015-2020.nc
      cf_2015-2020.nc.json
      ...
      cf_2056-2060.nc
      cf_2056-2060.nc.json
    audit.json
  index/
    units.json
    authoritative_index.json
    coverage_summary.csv
    unmapped_stations.csv.gz
    validation_summary.json
```

该目录仅定义数据产品，不表示本次已创建远程目录。只读源网格 CF，站点结果不写入源运行目录。相对路径以对应索引的运行根解析，同时保留原始源绝对路径及身份用于追溯。

不默认生成合并的 46 年单文件，也不按每个站点生成一个小文件。目录中 Climate/Station 两层缺一不可。

### 9.2 年份 NC schema

| 字段 | 维度 / 类型 | 含义 |
|---|---|---|
| `time` | time / 源 dtype | 源值及 units/calendar |
| `station` | station / int32 | 文件内序号 |
| `station_id` | station / 字符串 | 由 Station/tech/位置生成的稳定 ID |
| `lon`,`lat` | station / float64 | 场站位置 |
| `source_grid_lon`,`source_grid_lat` | station / float64 | 被选择格点位置 |
| `source_iy`,`source_ix` | station / int32 | 本源 patch 的零基行列索引 |
| `mapping_status` | station / int8 | MATCHED 或 OUTSIDE_DOMAIN；全状态表在 mapping 中 |
| `domain_mask` | station / int8 | 选中源格点的 domain 状态 |
| `match_dist_deg`,`match_dist_km` | station / float32 | 角距离和距离诊断 |
| `wind_cf` 或 `solar_cf` | time,station / float32 | 单位 `1`，有效范围 `[0,1]`，源 CF 的缺测语义 |

`first_snapshot_year` 作为辅助坐标记录首次出现的快照年份，不代表投运年份，也不用于屏蔽 CF。完整容量记录放 catalog，避免在每个气候组合/年份文件重复写入容量历史。

全局属性至少包含：`schema_version=station-cf-v1`、model、climate_scenario、station_scenario、station_source_label、tech、source_patch、source_run_id、source_cf_identity、source_block_identity、catalog_identity、mapping_identity、unit_identity、block_identity、match_method、max_distance_deg、capacity_semantics=snapshot_total、capacity_time_mask=off、years、代码 SHA 和编码配置。源物理参数及 BCSD 来源通过 source provenance 引用保留。

CF 变量保留 `wind_cf`/`solar_cf` 名称，避免下游误将容量系数当发电量；设置 units、valid_range、坐标属性和 fill。文件必须能被 netCDF4 和 xarray 正确读取，并验证字符串 ID round-trip。

### 9.3 Chunk 与压缩

初始输出候选：`chunksizes=(min(240,n_time),min(1024,n_station))`、float32、zlib level 2、shuffle 开启。一个完整 240×1024 CF chunk 约 960 KiB。不使用量化或有损压缩。

性能测试比较时间 chunk 240/960、station chunk 256/1024，以及压缩 1/2/4；同时测按时间批量和按少量站点读取。先确定吞吐与空间的平衡，再冻结一个编码配置，不在处理中自动改变已发布文件编码。

### 9.4 空 unit、缺覆盖与完整性

- 有来源分配的站点数为 0：生成 `EMPTY_NO_STATIONS` manifest/audit，不创建八份空 NC；仍计为 3,384 个逻辑 unit 中的一个已处理项。
- 有站点但全为 OUTSIDE_DOMAIN：写八份全 fill CF，状态为 COMPLETED，覆盖报告标明没有有效域内站点。
- TOO_FAR/NO_SOURCE_PATCH：保留全局 catalog/mapping 和原因，由 reader 在相应站点序列补 fill；不伪造 source_patch。
- 源 CF 文件应存在却缺失/损坏：是输入失败，不能归类为 EMPTY_NO_STATIONS 或用全 fill 文件替代。

全局完成同时要求所有逻辑 unit 的终态和所有 catalog 站点的覆盖说明。处理完成不等于 100% 场站具有有效 CF。

### 9.5 原子发布和恢复

每个分片使用唯一 `.partial.<pid>.<uuid>`，写完关闭并校验后同目录原子 rename 为 NC，再原子写 sidecar，最后父进程更新 manifest。sidecar 是该分片完成证据的一部分；仅 NC 存在不能跳过。

unit 入口持独占锁；分片也持独占输出锁。恢复仅复用身份、时间、station 顺序、文件 stat 和 sidecar 完整匹配的已完成分片。孤立 NC/partial 记为不完整，不覆盖已发布权威结果；需重做时采用新的明确尝试路径，manifest 只选择一个有效产物。科学/映射/编码身份变化使用新的运行根或版本，不自动加 overwrite。

所有数组工作结束后，父进程逐份复核八段并写 COMPLETED manifest。全局索引由一个串行汇总步骤在所有 unit 证据齐备后发布，子进程不并发改全局 JSON。

### 9.6 读取接口

```python
open_station_cf(
    index_path,
    model="CANESM5",
    climate_scenario="ssp126",
    station_scenario="ssp585",
    tech="wind",
    station_ids=[...],
    years=(2030, 2039),
)
```

两个情景参数显式必填。先根据索引找到 source_patch 和相交年份，再读取所需分片。按用户给定 station_id 顺序返回，并携带 mapping_status 与缺测。全球目录中合法但未覆盖的站点补 fill；不在指定 Station catalog 的 ID 应报错，不视为已知缺覆盖。

同组合跨 patch 时间轴必须一致才能拼接；缺覆盖站点使用相同 model/Climate/tech 的已验证参考时间合同恢复时间。若空间分片时间不一致则明确失败，不用外连接自动填出一个看似完整的时间轴。

提供 `iter_station_cf(...)` 按 patch/年份/时间块迭代，供全球统计和损失计算使用；全局默认访问不一次性 load 全部数据。跨模型、跨技术保留原生日历与时间差异，不强制合并为统一 datetime64 轴。

## 10. 与极端事件、装机和损失的接口

这次交付是场站 CF，不计算事件或损失，但身份设计需让后续可以可靠连接：

| 数据 | 情景决定方式 | 连接依据 |
|---|---|---|
| CF | 气候数值来自 Climate，位置来自 Station | model、Climate、Station、tech、station_id、time |
| 场站极端事件 | 源事件来自 Climate，位置也必须来自同一 Station | 同上，并核对源网格/映射和事件有效性 |
| 装机容量 | 按 Station 和独立 capacity_snapshot_year 选取累计装机快照 | Station、tech、station_id、snapshot_year |
| 损失/暴露 | 三者连接后的分析结果 | 双情景及具体损失定义/容量规则 |

固定 Station、改变 Climate，可以分析固定布局和容量下的气候差异；固定 Climate、改变 Station，可以分析布局/容量改变下的结果。后者若只比较同站 CF，需显式使用共同位置集合；全球总量比较可以包含布局差异，但必须报告站点/容量覆盖和分母。

参考极端事件项目现有对应情景输出不能直接充当全部交叉情景事件输入。例如 Climate ssp126 × Station ssp585 所需的是 ssp126 网格事件在 ssp585 场站位置的提取，不能把已有 ssp126 场站事件改标签，也不能读取 ssp585 气候事件替代。

后续事件提取器应使用同样的双情景和 Station ID 契约，并核对其最近邻 tie-break、源网格指纹和 domain。相同方法名并不足以证明映射完全相同；发现不同格点来源应先解决再计算联合指标。

示意发电量换算：`CF × installed_capacity_gw × Δt_hours = GWh`。实际损失还需事件、反事实基准和缺测策略，不能仅从 CF 直接得出。时间间隔从源时间合同确定；有效 CF/事件与缺测分别处理，不能在连接时将缺测填 0 或 False。不同日历只在有明确分析规则后再对齐。

本产品的 CF 不应用 RQ3 中针对错误放大风电发电量/损失字段的临时 `×0.1` 修正；容量因子复制的数值以本轮网格 CF 为准。

## 11. 本地测试、性能验证与验收

### 11.1 必要测试

| 类别 | 重点案例 |
|---|---|
| 双情景 | 九种配对齐全，3,384 ID 唯一；off-diagonal 选正确气候文件/Station CSV；两情景均进入路径和身份 |
| Station 别名 | ssp560 只归一化 Station；不增加任务数；Climate 拒绝 ssp560 |
| catalog | 跨年同址、零容量、坏行、重复键、四位 ID 冲突、年度原始行守恒、快照不累加、缺快照报错、first_snapshot_year 非投运年 |
| 最近邻 | 精确/半格等距、升降序、跨 patch、经度首尾、北极、缺 patch、门槛边界 |
| 覆盖 | OUTSIDE_DOMAIN 不向内陆吸附，TOO_FAR/NO_SOURCE_PATCH 不消失，空 unit 与全 fill 区分 |
| 数值 | CF 0/1、fill/NaN、多站同格、配对索引非笛卡尔积、原有效值精确相等、不 clip |
| 时间 | noleap/365_day、公历、闰日、源 wind/solar 首点差异、八段边界、无隐式重采样 |
| 并行 | 真实 CLI 的 1/2/4/8 进程科学数据一致；不同完成顺序不改变 manifest block 顺序 |
| 恢复 | 已完成跳过、partial、缺 sidecar、源 stat 改变、identity 冲突、子进程失败、输出锁竞争 |
| reader | 双情景筛选、跨 patch ID 重排、跨年子集、未覆盖补 fill、未知 ID 拒绝 |
| 反事实不变量 | 固定 model/Climate/tech/源格点，不同 Station 取值完全相同；固定 Station 跨 Climate ID 不变 |

测试使用小型合成 CF NC 和场站 CSV，不需要真实 BCSD 数据。比较不同进程输出时比较数值/坐标/身份，不要求包含耗时等属性的整个文件二进制一致。

验证命令：

```bash
.venv/bin/python -m pytest \
  tests/test_station_cf_catalog.py tests/test_station_cf_mapping.py \
  tests/test_station_cf_extract.py tests/test_station_cf_reader.py -q
```

再运行与实际复用工具相关的现有 `tests/test_grid_cf.py` 和 `tests/test_station_gather_crop.py` 回归测试。上述新测试文件已实现；当前已完成合成 CF 端到端测试和现有网格 CF 回归测试。

### 11.2 真实 CF 样本性能验证

实现后选择风/光、两类日历、最密/稀疏 patch、边界/缺测样本及一个非对角情景组合：先单年份段，再完整八段。测量：

- 1/2/4/8 进程的墙钟时间、总吞吐与 RSS；记录源读、gather、输出写分别耗时。
- 实际站点数、唯一源格点数、选中空间 chunk 数、每 chunk 复用量。
- 输出大小、压缩比、按小站点集/连续时间批次的读性能。
- 小样本全部逐元素验证；真实大样本使用固定种子独立重读，覆盖每年份段和边界/缺测站点。

只做独立抽样的结果标明为抽样，不宣称对所有 TB 级文件二次全量重读。正式提取时对每个流块检查有效值域和缺测统计，独立审核补充检验源位置/时间取值。

### 11.3 完成标准

单个非空 unit：八段 NC/sidecar 完整、身份正确、站点顺序一致、时间与源相等、源有效值提取一致、所有文件可读，manifest/audit 为 COMPLETED。空 unit：catalog/mapping 身份验证通过，来源分配数为 0，manifest/audit 为 EMPTY_NO_STATIONS。

全局审核：

1. 恰好 3,384 个任务身份，无遗漏/重复；其中 `COMPLETED + EMPTY_NO_STATIONS = 3,384`，失败/待处理为 0。
2. 真实年份 NC 数为 `8 × 非空已完成 unit 数`，上限 27,072；八段总 time 数与相应源一致。
3. 对每个 model/Climate/Station/tech，有 `catalog_count = MATCHED + OUTSIDE_DOMAIN + TOO_FAR + NO_SOURCE_PATCH`。前两类在各 source_patch 文件中恰出现一次；后两类在 mapping/未覆盖表中有完整记录。
4. 每站 `valid_count + missing_count = 源 time_count`，有效 CF 值域正确，domain 外整列 fill。
5. 每个源格点多站结果相同，且反事实不变量通过。
6. 输出覆盖报告同时展示数据处理完成率、有效映射率、时间缺测率，不能仅报文件数。
7. 全局 reader 能按两情景读取样例并连接相应 Station 容量记录；不依赖文件行序或旧单情景路径。

## 12. 存储预算与实施顺序

### 12.1 存储估算

不能按源网格 9.12 TB 简单乘三得到输出空间。站点输出取决于目录规模、覆盖率、同格多站重复、原生 time 长度及压缩比。

准确的 float32 逻辑体积为：

```text
B_raw = 4 × Σ(model,Climate,Station,tech,patch)
                [N_output_station × Σ(year_block) N_time]
```

输出站点数包括 OUTSIDE_DOMAIN；没有来源的站点用 catalog/mapping 表示，不占 CF 矩阵实体列。CSV、索引、映射、临时分片和失败尝试空间另计。

参考文档称三 Station 的 wind/solar 唯一位置合计约 1,606,397，这仅作为算例。假设每模型/Climate/tech 均约 134,400 个三小时时点且全部位置写入，则未压缩 CF 约为：

```text
4 bytes × 4 models × 3 Climate × 134400 × 1,606,397
≈ 10.36 TB（十进制）
```

三 Station 已包含在 1,606,397 中，不能再乘一次 3。该数不是实际文件预测；准备阶段以本批目录与源时间长度重算，性能样本测压缩后 bytes/元素，并按密集/稀疏、风/光和缺测比例分组外推。输出数据空间建议另外预留临时/恢复开销及至少 20% 余量。本方案不据此承诺具体运行时长。

### 12.2 实施步骤和交付物

| 顺序 | 开发内容 | 阶段验收 |
|---|---|---|
| 1 | catalog、双情景配置/别名、输入 CF 索引、任务身份 | 六套 catalog；九配对；3,384 唯一 unit；容量原始记录完整 |
| 2 | 全局最近邻与覆盖映射 | 边界/缺 patch/domain 测试通过；每站去向可解释 |
| 3 | 单进程单年份 CF gather/writer | 数值/缺测/time 与源一致，原子产物闭环 |
| 4 | 八段 spawn、manifest、锁和恢复 | 1/2/4/8 进程一致，失败后复用正确 |
| 5 | 索引、全球 reader 和双情景下游连接 | 任意配对/时间/站点可读取；缺覆盖与未知 ID 区分 |
| 6 | 真实 CF 样本、性能/存储评估、独立审核 | 冻结参数建议、报告限制与覆盖，不仅验证程序退出码 |

最终代码交付：上述模块/入口、必要测试、小型样本配置及读取示例。数据交付：场站 CF 年份分片、catalog 与原始容量行、mapping、逐 unit manifest/audit、唯一权威索引、覆盖和验证报告。超算运行材料属于后续单独工作，不包含在本实施文档对应的改动中。

## 13. 容量语义证据与实现验证项

### 13.1 上游代码证据（2026-10-01 核查）

上游仓库：`/data6/yanxiaokai/project_climate/globally-interconnected-solar-wind-system-addresses-future-electricity-demands`。本次实际核对了 ssp126、ssp245、ssp560 三份实现，而非只依据函数注释或成本模型中的 incremental 名称判断。

| 证据位置（相对于上游仓库） | 实际操作 | 含义 |
|---|---|---|
| `Optimization_10km_ssp126/downscale_stage.py:48–54` | `dC = C_star - prev_y.sum()`；`y_cum = prev_y + dy` | 求解新增量后回加到累计状态 |
| `Optimization_10km_ssp126/downscale_stage.py:64–94` | `save_year` 从 `cum[tech]` 取 `y`，保存 `y[sel]` 到 `*_capacity_twp` | 年度 NPZ 保存累计装机，不是 dy |
| `Optimization_10km_ssp245/downscale_stage.py:49–55,65–96` | 相同的差分求解、累计更新和累计输出 | ssp245 同口径 |
| `Optimization_10km_ssp560/downscale_stage.py:49–55,65–96` | 相同的差分求解、累计更新和累计输出 | ssp560 同口径 |
| `plot_10km_optimal_stations.py:85–109` | 按像元聚合年度 NPZ 的 `*_capacity_twp`，乘 1000 转 GW | 仅空间聚合及单位换算，没有年度差分 |
| `plot_10km_optimal_stations.py:190–205,225–240` | 每年加载后直接输出 `year,type,lon,lat,capacity_gw` | CSV 仍是年度累计装机快照 |

核查时 `downscale_stage.py` SHA256：

- ssp126：`5380f88eeee83c9c4390fddbb74d575ffb0670033e8f151e54af39a0ad438f29`。
- ssp245 与 ssp560：`e82584672e2d7960fd11110c4a862b0915784bb7296df131d64e4ce9ee59c94b`（两份文件一致）。

参考事件方案的“按 year≤目标年累计 capacity_gw”与上述上游快照语义不一致，新 CF 方案使用本节核实的快照契约。现有事件 CF 数值提取若未应用容量/投运 mask，不因这一语义差异自动失效；任何按该规则计算过的装机权重、暴露或损失，后续接入前需单独核查。本次只修订本仓库方案，不修改其他项目结果。

### 13.2 输入冻结时的来源核对

上游导出脚本默认文件名为 `stations_10km.csv`，用户指定远程文件为 `stations_SSP*.csv`；仅凭名称不能证明复制/重命名过程未变更容量。本次已确认上游生产/导出口径，尚未全量比对这三份远程 CSV。实现准备阶段需记录实际来源批次，并按 `(year,type,lon,lat)` 核对记录数、坐标集合及逐位置容量，验证远程 CSV 与对应年度 NPZ/导出 CSV 相符；若仅改名且内容相同，可用文件 hash 确认。CSV 保留四位小数，NPZ 对照按该导出舍入规则比较；总量误差容差需计入各行舍入累计。

核对不通过就报告来源或语义冲突，不自动把它当增量重解释。通过后将 `capacity_semantics=snapshot_total`、快照年份列表和来源 hash 写入 catalog/prepared 身份。该项是可执行的数据验证，不再请求用户凭经验判断增量或快照。

### 13.3 其他实现验证项

以下属于实现验证项，不需用户预先选择：远程 CSV 当前 hash/数量、空间组实际数量、真实日历/时间长度、输出 chunk 和最优进程数。根据代码验证和样本实测冻结，不以经验数替代证据。

本方案依据：

- 本仓库 [网格 CF 主入口](../patchify_grid_cf.py)、[网格读写](../grid_cf_io.py)、[原场站入口](../patchify_station_cf.py)。
- 本地已验收的 `infos/scnet_patchify_grid/completion_status/completion.json`、`units.json`、`final_audit.json`（运行记录被 Git ignore，代码检出不自带这些文件）。
- 用户指定参考方案：`/data6/yanxiaokai/project_climate/extreme_event_definitions/infos/scnet_patchify_stations/实施方案.md`。
- 参考项目实际 `grid_extreme_signals/station_catalog.py` 与 `grid_extreme_signals/station_match.py`，用于核对场站字段和 ID 约定；容量语义以上游 Optimization 代码为准。

本次已实现准备、提取、恢复、审核发布和惰性读取。合成数据验证包含真实 CLI 的 1/2/4/8 进程、双情景数值不变量、光伏时间偏移、快照容量、空 unit 和异常处理。未执行新的远程数据处理；真实输入核对、性能结果和全量生产状态须另行记录，不能由本地测试推定。
