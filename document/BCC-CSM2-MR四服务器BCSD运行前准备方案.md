# BCC-CSM2-MR 四服务器运行前准备方案

> 调查与实施日期：2026-08-28。
>
> 本文同步记录运行前准备实施结果：四台模型软链和全部作业脚本已经建立，监控器和本地模式
> 资料已经就绪。没有提交或取消作业，没有修改远端代码、ACL 或共享缓存，也没有全量读取
> NetCDF 数据数组。
>
> 目标模式：`BCC-CSM2-MR`；训练期 `1979-2014`；应用期 `2015-2060`；情景为
> `ssp126/ssp245/ssp585`；变量为 `rsds/pr/tas/uas/vas/hurs`。

## 1. 结论

本轮任务规模为：

```text
26 regions × 3 SSP × 6 variables = 468 个最终输出
```

运行前调查与实施结论如下：

1. BCC-CSM2-MR 的 GCM 实体数据已经存在于 185 账号统一目录中，六个目标变量共有 `410`
   个文件、约 `649.358 GiB`，无零字节文件。文件名时间范围覆盖训练期 `1979-2014` 和
   应用期 `2015-2060`；抽样元信息中的变量名、维度、单位、日历、时间戳和 `160 × 320`
   网格符合当前流程要求。
2. 四个生产账号均可读取中央目录；四台的 `~/bcsd/data/CMIP6/BCC-CSM2-MR` 已建立为指向
   185 统一目录的软链，默认 `--gcm-base "$HOME/bcsd/data/CMIP6"` 现可找到该模式。
3. 四服务器区域分配复用
   `document/MPI-ESM1-2-HR四服务器运行与共享缓存方案.md` 第 4.2 节，不调整区域归属。
4. 共享中间数据 ACL 方案已经实现，本轮只复用，不再设计、不迁移、不改权限。
5. Step1 的 156 个“区域 × 变量”组合经主流程同款区域候选逻辑与元数据级 resolver 复核后为
   `87 HIT + 69 HIT COMPAT + 0 MISS`，全部可复用。
6. 时间偏移口径：BCC 与 MRI 一致，`rsds/pr=90` 分钟，`tas/uas/vas/hurs=0` 分钟。
   `tas/uas/vas` 不能复用 CANESM5 或 MPI-ESM1-2-HR 的 180 分钟 Step1；当前可命中
   MRI 已生成的 0 分钟共享缓存。
7. 已生成 226 个单作业脚本和 16 个 manual-rolling 包（48 个三阶段物理脚本），合计 274
   个物理脚本。静态验证覆盖 468/468 个最终组合且无重复，所有脚本均通过 `bash -n`。
8. 生产顺序：先运行全部区域的 `ssp126`，其中 India、México、Australia、Brazil 三阶段区域
   在各自服务器优先；任一区域 `ssp126` 六变量全部 `VALID` 后，再解锁该区域 `ssp245/ssp585`。

## 2. 固定任务范围与账号

### 2.1 区域全集

```text
Australia, Austria, Brazil, Chile, Denmark, Egypt, France, Germany,
Greece, India, Ireland, Italy, Japan, México, Netherlands, Poland,
Portugal, Romania, South Africa, South Korea, Spain, Sweden, Turkey,
Ukraine, United Kingdom, Vietnam
```

每个区域均运行：

```text
SSP：ssp126、ssp245、ssp585
变量：rsds、pr、tas、uas、vas、hurs
训练期：1979-2014
应用期：2015-2060
```

同一区域的三个 SSP 和六个变量必须留在同一账号，不允许跨服务器拆分或重复运行。

### 2.2 四个生产账号

| 名称 | SSH 主机 | 账号 |
|---|---|---|
| 乌镇151 | `scnet-wuzhen-151` | `acw4hq2non` |
| 乌镇173 | `scnet-wuzhen-173` | `acrxpqssp4` |
| 乌镇9259 | `scnet-wuzhen-9259` | `acnqt6koeq` |
| 乌镇9359 | `scnet-wuzhen-9359` | `ac6xcutdaj` |

## 3. BCC-CSM2-MR GCM 输入调查

### 3.1 实体数据与生产入口

实体数据已经存在：

```text
/work/home/acbw9wpn5k/bcsd/data/CMIP6/BCC-CSM2-MR
└── 3hr
    ├── historical
    ├── ssp126
    ├── ssp245
    └── ssp585
```

四个生产账号当前状态均为：

```text
~/bcsd/data/CMIP6/BCC-CSM2-MR
-> /work/home/acbw9wpn5k/bcsd/data/CMIP6/BCC-CSM2-MR
```

### 3.2 六变量文件清单与时间覆盖

| 实验 | 变量 | 文件数 | 文件名起止时间 | 本轮年份覆盖 | 零字节 |
|---|---|---:|---|---|---:|
| historical | `rsds/pr` | 每变量 22 | 1950-01-01 01:30 ～ 2014-12-31 22:30 | 覆盖 1979-2014 | 0 |
| historical | `tas/uas/vas` | 每变量 22 | 1950-01-01 00:00 ～ 2014-12-31 21:00 | 覆盖 1979-2014 | 0 |
| historical | `hurs` | 12 | 1979-01-01 00:00 ～ 2014-12-31 21:00 | 精确覆盖 1979-2014 | 0 |
| ssp126 | `rsds/pr` | 每变量 16 | 2015-01-01 01:30 ～ 2062-12-31 22:30 | 覆盖 2015-2060 | 0 |
| ssp126 | `tas/uas/vas/hurs` | 每变量 16 | 2015-01-01 00:00 ～ 2062-12-31 21:00 | 覆盖 2015-2060 | 0 |
| ssp245 | `rsds/pr` | 每变量 16 | 2015-01-01 01:30 ～ 2062-12-31 22:30 | 覆盖 2015-2060 | 0 |
| ssp245 | `tas/uas/vas/hurs` | 每变量 16 | 2015-01-01 00:00 ～ 2062-12-31 21:00 | 覆盖 2015-2060 | 0 |
| ssp585 | `rsds/pr` | 每变量 16 | 2015-01-01 01:30 ～ 2062-12-31 22:30 | 覆盖 2015-2060 | 0 |
| ssp585 | `tas/uas/vas/hurs` | 每变量 16 | 2015-01-01 00:00 ～ 2062-12-31 21:00 | 覆盖 2015-2060 | 0 |

按实验汇总：

| 实验 | 六变量文件数 | 大小 |
|---|---:|---:|
| historical | 122 | 194.497 GiB |
| ssp126 | 96 | 151.629 GiB |
| ssp245 | 96 | 151.615 GiB |
| ssp585 | 96 | 151.617 GiB |
| **合计** | **410** | **649.358 GiB** |

### 3.3 抽样 NetCDF 元信息

本次只读取了 historical 和 ssp126 每个目标变量各一个文件的头部元信息，以及时间、纬度、
经度坐标的首尾标量，没有读取数据数组。

| 变量 | 数据维度 | 单位 | 日历 | GCM 时间槽 |
|---|---|---|---|---|
| `rsds` | `(time, lat, lon)` | `W m-2` | `365_day` | 01:30、04:30、…、22:30 |
| `pr` | `(time, lat, lon)` | `kg m-2 s-1` | `365_day` | 01:30、04:30、…、22:30 |
| `tas` | `(time, lat, lon)` | `K` | `365_day` | 00:00、03:00、…、21:00 |
| `uas` | `(time, lat, lon)` | `m s-1` | `365_day` | 00:00、03:00、…、21:00 |
| `vas` | `(time, lat, lon)` | `m s-1` | `365_day` | 00:00、03:00、…、21:00 |
| `hurs` | `(time, lat, lon)` | `%` | `365_day` | 00:00、03:00、…、21:00 |

六变量抽样网格一致：

```text
lat = 160，范围约 -89.141519 ～ 89.141519
lon = 320，范围 0.000000 ～ 358.875000
```

### 3.4 manifest 状态

`historical/ssp126/ssp245/ssp585` 四个目录当前均没有 `manifest.tsv`。当前 `step02/step03`
输入发现使用文件 glob，不读取该清单，所以不阻塞运行；若未来把 SHA256 manifest 作为强制
审计依据，应重新生成完整清单，不能手填。

## 4. 四服务器区域分配

本表直接复用 MPI-ESM1-2-HR 方案第 4.2 节。区域归属、区域数量、网格总数和执行形态均不变。

| 服务器 | 单作业区域 | 三阶段区域 | 区域数 | 网格总数 | 负载占比 |
|---|---|---|---:|---:|---:|
| `scnet-wuzhen-151` | Greece、Romania、Germany、Spain、Chile、South Africa、Japan、Sweden、France | — | 9 | 158890 | 22.16% |
| `scnet-wuzhen-173` | Austria、Italy、Turkey | India、México | 5 | 185573 | 25.88% |
| `scnet-wuzhen-9259` | Denmark、Ukraine、United Kingdom、Poland、Ireland | Australia | 6 | 184495 | 25.73% |
| `scnet-wuzhen-9359` | Portugal、Netherlands、South Korea、Vietnam、Egypt | Brazil | 6 | 188155 | 26.24% |
| **合计** | **22 个** | **4 个** | **26** | **717113** | **100%** |

四台当前均未发现 BCC-CSM2-MR 的最终 NetCDF；468 个最终组合全部需要生成 BCC 自身最终结果。

## 5. 中间结果复用与时间偏移

### 5.1 Step1 复用判定

Step1 只处理训练期 ERA5-Land。当前语义缓存键比较：

```text
region、bbox、ERA5 输入变量、输出变量、训练年月、3 小时步长、
target_offset_minutes、ERA5 源文件 identity、转换算法版本
```

模式名本身不决定新语义键，实际时间偏移决定能否跨模式复用。旧 artifact 由兼容 resolver
转换为相同语义后再判断。

### 5.2 时间偏移对照

| 输出变量 | ERA5-Land 变量 | CANESM5 | MPI-ESM1-2-HR | MRI-ESM2-0 | BCC-CSM2-MR |
|---|---|---:|---:|---:|---:|
| `rsds` | `ssrd` | 90 分钟 | 90 分钟 | 90 分钟 | 90 分钟 |
| `pr` | `tp` | 90 分钟 | 90 分钟 | 90 分钟 | 90 分钟 |
| `tas` | `t2m` | 180 分钟 | 180 分钟 | 0 分钟 | 0 分钟 |
| `uas` | `u10` | 180 分钟 | 180 分钟 | 0 分钟 | 0 分钟 |
| `vas` | `v10` | 180 分钟 | 180 分钟 | 0 分钟 | 0 分钟 |
| `hurs` | `rh` | 0 分钟 | 0 分钟 | 0 分钟 | 0 分钟 |

结论：

- `tas/uas/vas` 不能复用 CANESM5 或 MPI-ESM1-2-HR 的 180 分钟 Step1；
- BCC 与 MRI 时间偏移完全一致，可复用 MRI 已生成且语义校验通过的 Step1；
- 不能通过重命名缓存目录、修改 manifest 或平移时间坐标伪造命中。

### 5.3 元数据级 resolver 实测

复核使用 `Step1CacheResolver(full_validation=False)`，只读取 manifest，并检查 artifact 路径、
变量、大小、mtime 和 schema；没有读取 NetCDF 数据数组。

| 服务器 | Step1 总数 | HIT | HIT COMPAT | MISS | 备注 |
|---|---:|---:|---:|---:|---|
| `scnet-wuzhen-151` | 54 | 34 | 20 | 0 | South Africa 使用 `South-Africa` ERA5 目录 |
| `scnet-wuzhen-173` | 30 | 15 | 15 | 0 | — |
| `scnet-wuzhen-9259` | 36 | 20 | 16 | 0 | United Kingdom 使用 `United-Kingdom` ERA5 目录 |
| `scnet-wuzhen-9359` | 36 | 18 | 18 | 0 | South Korea 使用 `South-Korea` ERA5 目录 |
| **合计** | **156** | **87** | **69** | **0** | **全部可复用** |

三区域的复核结论：

```text
South Africa   payload region = "South Africa"，ERA5 目录 = South-Africa
United Kingdom payload region = "United Kingdom"，ERA5 目录 = United-Kingdom
South Korea    payload region = "South Korea"，ERA5 目录 = South-Korea
```

作业脚本中的 `--region` 继续使用带空格的区域名；主流程会按 `region_raw` 和 `region_fs`
候选路径查找 ERA5-Land 输入。

### 5.4 Step2 以后不能跨模式复用

| 阶段或产物 | 旧模式/MRI → BCC | 原因 | BCC 内部跨 SSP |
|---|---|---|---|
| Step1 ERA5 高频观测 | 可复用 | 由显式时间偏移和 ERA5 源 identity 决定 | 可以 |
| Step2 观测到 GCM 网格 | 不可复用 | payload 包含模式名、BCC GCM 网格样例和 Step1 identity | 可以 |
| historical Step3 | 不可复用 | payload 包含模式名及 BCC historical GCM 文件 identity | 可以 |
| apply Step3 | 不可复用 | 模式和 SSP 均不同，且产物位于本次运行目录 | 不跨 SSP |
| Step4 粗网格订正 | 不可复用 | 依赖当前模式和当前 SSP 的 GCM | 不跨 SSP |
| Step5 SD map | 不可复用 | 依赖 BCC 自身的 Step1/Step2 key 和文件 identity | 可以 |
| Step6 remap 网格/权重 | 不可跨模式复用 | 缓存路径按 `<model>/<region>` 隔离 | BCC 内部可以 |
| Step6 分块与最终输出 | 不可复用 | 每个模式、SSP、变量均需新输出 | 不跨 SSP |

## 6. 作业脚本与监控准备状态

- 单作业脚本位于四台的 `~/jobs/bcc_csm2_mr/`，共 226 个；
- 三阶段包位于 173、9259、9359 的 `~/jobs/bcc_csm2_mr_manual_rolling/`，共 16 个包、48 个
  物理脚本、72 条变量链；
- 三阶段区域拆为 `ssp126_time0`、`ssp126_shared`、`ssp245`、`ssp585`；
- 274 个物理脚本全部通过 `bash -n`；
- 脚本/plan 参数覆盖 468/468 个最终组合，无重复、无遗漏；
- `infos/scnet/monitor_bcsd_jobs.py` 已登记 BCC 四服务器矩阵，支持普通作业、manual-rolling、
  最终输出校验和固定优先顺序；
- 本地 `infos/scnet/models/BCC-CSM2-MR/` 已创建，包含 `config.md`、`slurm_scripts.md` 和
  留空具体进度的 `completion_status/progress_summary.md`；
- 本轮未启动常驻监控，也未提交任何作业。

## 7. 强制生产顺序

### 7.1 第一阶段：全部区域 ssp126

第一阶段运行 `26 × 6 = 156` 个 `ssp126` 最终任务。India、México、Australia、Brazil 的
三阶段区域在各自服务器优先，分别运行 `ssp126_time0` 和 `ssp126_shared`。

普通区域提交对应的六个单变量脚本：

```text
job_<REGION>_BCC-CSM2-MR_ssp126_<VAR>.sh
```

### 7.2 第二阶段：后续 SSP

某区域 `ssp126` 六变量全部通过后，再提交该区域的 `ssp245/ssp585`。预期缓存行为为：

```text
Step1：全部命中 BCC 已可用 artifact
Step2：命中 BCC ssp126 生成的缓存
historical Step3：命中
Step5：命中
Step3 apply、Step4、Step6 最终输出：按 SSP 新计算
```

## 8. 正式提交前验收清单

### 8.1 GCM 输入

- [x] 四台 `~/bcsd/data/CMIP6/BCC-CSM2-MR` 均解析到 185 统一目录；
- [x] 四台均能读取中央 GCM 目录；
- [x] 410 个目标 GCM 文件不可读和零字节数量均为 0；
- [x] historical 六变量覆盖 1979-2014；三个 SSP 六变量覆盖 2015-2060；
- [x] 不把缺失的 `manifest.tsv` 当作 glob 运行阻塞项。

### 8.2 缓存与时间偏移

- [x] 四台共享缓存软链均解析到同一 ACL 中央目录；
- [x] BCC 偏移配置为 `rsds/pr=90`、`tas/uas/vas/hurs=0`；
- [x] 元数据预检得到 156 个可复用 Step1 组合；
- [x] 生成脚本明确传入 `--cmip6-model BCC-CSM2-MR`；
- [ ] 首批实际日志确认 `ssp126` 六变量 Step1 均命中语义校验通过的 artifact；
- [ ] 后续 SSP 日志确认六变量 Step1 继续命中 BCC 已有 artifact。

### 8.3 区域与作业

- [x] 四服务器区域表与 MPI 第 4.2 节完全一致，26 个区域无重复、无遗漏；
- [x] 同一区域的三个 SSP 和六变量没有跨账号；
- [x] 22 个区域使用单作业，India、México、Australia、Brazil 使用 manual-rolling 三阶段；
- [x] 作业脚本和 plan 中的模式、区域、情景、变量、年份和输出路径均正确；
- [x] 四台没有 BCC 在队作业或既有最终输出；
- [x] 274 个物理脚本通过语法检查，468 个最终组合无重复、无遗漏；
- [x] 监控器已经支持 BCC，并使用独立的 BCC completion/usage 状态目录；
- [ ] 任一区域 `ssp126` 六变量全部 `VALID` 后才解锁该区域后续 SSP；
- [ ] 提交前复核余额、活动作业数、脚本语法和实际核数。

## 9. 本轮边界与证据口径

本轮远端操作包括：

- 目录、软链、文件名、数量、大小、权限和 manifest 的只读检查；
- historical 与 ssp126 每变量一个文件的 NetCDF 元信息抽样，只读取头部和坐标首尾标量；
- Step1 resolver 的 156 组 `full_validation=False` 元数据预检；
- 四台 BCC 既有输出、脚本和作业数量检查；
- 四台建立 BCC GCM 软链；
- 四台生成 226 个单作业脚本和 16 个 manual-rolling 包；
- 静态解析全部脚本与 plan，并对 274 个物理脚本执行 `bash -n`。

本轮没有执行：

- 全量读取任一 GCM 或 ERA5-Land NetCDF 数据数组；
- 修改 ACL、共享缓存或远端代码；
- 执行 `sbatch`、启动常驻监控、取消作业或删除任何文件；
- 写入或改动任何 GCM、ERA5-Land、缓存、运行产物或最终输出。

下一步仅在再次完成余额、重复作业和输出检查后，按第 7.1 节提交 `ssp126` 任务。
