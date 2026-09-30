# 全格点 CF 计算实施方案

编写日期：2026-09-30。分支：`develop-patch-grid`；代码基线：`c845d15`。

实现已落地，当前入口和验证结果见 [实现与本地验证记录](全格点CF实现与本地验证记录.md)。本文保留设计依据，具体CLI以 README 为准。本阶段范围为全格点容量因子（CF）计算、NetCDF 存储、分块与多进程、本地验证；不建设 `infos/scnet_patchify_grid/`，不规划该目录内部文件，也不涉及新超算作业的生成、提交和监控。

## 1. 推荐方案

1. 新增 `patchify_grid_cf.py`，在 BCSD 原生格点上计算风电、光伏 CF，复用 `cf_physics.py` 的公式和参数。新入口不需要场站 CSV、空间匹配、装机容量或场站筛选。
2. 最终数据为 `wind_cf(time, lat, lon)` 或 `solar_cf(time, lat, lon)`，float32，取值范围为 `[0,1]`；域外与气象输入无效的位置保存缺测，不能写成有效的零。
3. 一个 unit 为 `1 GCM × 1 SSP × 1 tech × 1 patch`，对应一个作业内的完整计算入口：父进程调度八个年份块任务，默认发布八个年份 CF NetCDF，作为正式产物。通过 argparse 开关 `--merge-final` 控制是否在同一作业内额外流式合并为 2015–2060 年的完整 CF NetCDF，默认不合并。
4. 生产输入固定为 BCSD `production-V2` 的八个年份 point block，由该批次 manifest 解析 `blocks/` 下的实际路径，同时读取同批次 land plan 和 sidecar。完整 BCSD final 仅用于测试参考，不作为生产计算输入。
5. 年份段之间使用一层 `spawn` 进程池，`--processes` 默认值为 8，最多八个 worker；每个 worker 内按时间 chunk 和空间 tile 流式读取、计算、写出。八个计算任务可由较少 worker 分批执行。启用 `--merge-final` 时，父进程在 worker 全部结束后独占 final 写入，不拆成年份作业分次提交。这里只定义计算单元，不实现超算脚本。

“所有格点”指上游 BCSD 计算域的所有格点，不再受场站位置限制。输出保留 patch 完整矩形经纬度，海洋及上游未计算的位置为缺测。现有 BCSD 没有提供的海洋、南极洲或被过滤区域，不能仅靠 CF 改造产生有效结果。

## 2. 当前实现及可复用部分

### 2.1 计算链路

当前生产入口为 `patchify_station_cf.py`，输出实际为 `(time, station)`。BCSD 年份块输入为 `(time, point)`，两者的 point/station 含义不同。

```text
场站 CSV → 技术筛选、经纬度去重、patch core 筛选
    ├─ final：读取规则网格 → nearest/bilinear 抽取场站天气
    └─ blocks：读取年份合同及 land plan → point 索引映射 → 抽取场站天气
        ↓
按技术参考时间轴对齐 → 单位转换 → 风光物理核
        ↓
流式写 (time, station)；blocks 模式最后合并完整时间段
        ↓
<output-root>/<model>/<scenario>/<patch>/<tech>.nc + .json
```

| 代码位置 | 当前行为 | 网格改造方式 |
|---|---|---|
| `find_final`、`_var`、`_coord` | 定位 final、解析变量名和坐标别名，检查 sidecar 的 patch/variable | 保留为测试参考；生产 reader 从 production-V2 的 block manifest 定位输入 |
| `_load_block_plan`、`_block_identity` | 读取 `time_block_contract.blocks`、`final_tasks[].block_files` 和 sidecar | 提取不依赖 stations 的输入计划，沿用真实 manifest 字段 |
| `load_stations`、`_select_patch_stations`、`_match`、`_point_map` | 场站筛选、空间抽取和距离阈值 | 新入口改为原生网格描述与 point→grid 索引映射 |
| `_read_block_chunk` | 按时间邻点插值后 gather 到 station | 复用时间邻点思想，重写连续读取和网格填充部分 |
| `_doy_hour` | 保留日序和含分钟、秒的 UTC 小时 | 复用；支持的日历需显式验证 |
| `_solar_points_block` | 逐站调用三维光伏核 | 改为每个矩形 tile 一次调用光伏核 |
| `cf_physics.py` | NumPy 风光计算，支持网格形状 | 保留公式，只在读写层处理有效性 |
| `_StreamingStationWriter`、`_compute_blocks` | station writer、年份 worker、最终合并 | 新建网格 writer、年份调度及流式 final 合并 |

新入口以 `--years 2015-2060` 定义一个完整 unit，年份块范围从 manifest 获取。断点恢复自动跳过身份匹配且校验通过的年份块，不通过拆分年份提交来完成正常运行。

### 2.2 物理定义保持不变

| 技术 | 输入和参考轴 | 当前公式要点 |
|---|---|---|
| wind | `uas, vas`；参考 `uas.time` | 10 m 风速模长，100 m、1/7 幂律外推，`GE120/2500` 功率曲线，25 m/s 切出，CF 裁剪至 `[0,1]` |
| solar | `rsds, tas, uas, vas`；参考 `rsds.time` | 辐射转 kW/m²、温度转 °C，Erbs 分解、双轴跟踪、温度修正，系统系数 0.8056，CF 裁剪至 `[0,1]` |

CF 是逐时、逐格点计算，时间插值以外没有跨时刻依赖。网格版直接将 `(T,Y,X)` 传入物理核；光伏坐标参数必须是 `lat(Y)`、`lon(X)`，不能把两个长度为 K 的点坐标当作两条网格轴，否则会广播成 `K×K`。

经纬度原样写出；太阳几何若为数值等价需要将经度映射到 `[-180,180)`，只对计算参数转换，不改变输出坐标。用跨 180° 的测试验证与原场站公式一致。

场站代码先插值天气再计算非线性 CF。新网格 CF 在格点中心可与旧入口对照；在任意站点插值网格 CF，通常不等于先插值天气再算 CF。现有场站 Loss 读取器后续需另行适配，本次不将其列为验收依赖。

### 2.3 当前超算运行方式：仅作背景

`scnet/create_cf_patch_jobs.py` 是转发入口，实际生成逻辑在 `infos/scnet_patchify/create_cf_patch_jobs.py`。每个脚本处理一个 `model × scenario × patch × tech`；patch 列表必须显式传入，生成器不调用 `sbatch`。

当前生成器默认 `wzhctest`、10 CPU，计算参数默认 `input-mode=final`、`processes=8`。其中年份多进程只在 blocks 路径生效，不能把参数默认值当作 final 已并行的证据。脚本激活共享 `climate` 环境，并把 BLAS/OpenMP 线程限制为 1；作业、日志和输出目录通过参数指定。

现有 CF 年份块文档记录了 2026-09-14 的八站小样本：4 worker 快于 8 worker，串并行 CF 最大差为零。该结果说明 I/O 和进程开销需要实测，不能作为全网格吞吐保证。运行方式依据本地源码和历史记录；另已只读确认 production_v2 的目录结构，未核验实时队列或生产完成状态。

## 3. BCSD 网格和输入合同

已检查本地上游 `../bcsd/global_bcsd/`，版本为 `de26ecb`：

- `land_plan.py` 的生产网格为 1801×3600、0.1°，纬度升序，经度 `[0,360)`；有效计算点来自 land plan，不从首时次有限值推断。
- `patches.py` 当前有 47 个活动 patch，通常为 30°×30°，排除部分纯海区域和南极洲。
- `finalize.py` 将 `(time, point)` block 按 `y_index/x_index` 放回 patch core 网格，final 为 `(time, lat, lon)`，海洋保留 fill。默认 chunk 是 `(240,64,64)`，压缩级别为 2；实际文件仍需查 header。
- 普通 patch 常见 300×300，北极 patch 可能为 301×300；从真实坐标读取维度，不写死。

生产数据位置（已通过 SSH 只读确认目录存在）：

```text
Host: scnet-wuzhen-1866（乌镇1866）
生产根目录: /work/share/aczlvkl1ac/bcsd_runs/production_v2/
完整 BCSD 输出: /work/share/aczlvkl1ac/bcsd_runs/production_v2/outputs/
年份块目录: /work/share/aczlvkl1ac/bcsd_runs/production_v2/blocks/
清单目录: /work/share/aczlvkl1ac/bcsd_runs/production_v2/manifests/
```

用户提供的 `outputs/` 是完整 BCSD 输出目录；`blocks/` 和 `manifests/` 与它同级。因此 `--bcsd-root` 传入其父目录 `production_v2/`，年份块及 land plan 的具体路径仍由该批次 manifest 解析。以上路径用于该远程文件系统上的运行，不是本地文件路径；本地验证使用小样本输入。

### 3.1 坐标来源和校验

生产网格：从同批次 land plan 的完整坐标轴和 patch core 范围取得输出轴；检查其与对应 final（如存在）逐元素相同。保留 land plan 的坐标 dtype 和顺序，不能用 `arange` 重建，也不能只取 block 中实际出现过的纬度/经度，否则会丢失整行或整列海洋。

使用 `y_index/x_index` 构建 point→patch 索引，校验索引范围、唯一性、`flat_index=y_index*N_global_lon+x_index`、`global_point` 对应的 land plan 点身份，以及 point 的实际 lat/lon。同一运行内各变量和各年份的点集合及顺序必须一致；首版遇到差异报错，不猜测重排。掩膜与该批次 land plan/patch 点集合交叉验证。

只写 core 网格，不将上游工作区 buffer 带入输出；保留 core 边界的半开区间及北极端点规则。

### 3.2 年份文件

从 `time_block_contract.blocks` 获取年份区间，当前合同为：

| 段 | 年份 |
|---|---|
| 00 | 2015–2020 |
| 01 | 2021–2026 |
| 02 | 2027–2032 |
| 03 | 2033–2038 |
| 04 | 2039–2044 |
| 05 | 2045–2050 |
| 06 | 2051–2055 |
| 07 | 2056–2060 |

所有变量合同须相同、连续且无重叠，源文件路径取自 `production-V2` 的 manifest，并核验为同一生产批次的 blocks。输入根目录通过参数传入，不能硬编码或沿用历史 production_v1 路径。八块及其 sidecar、land plan 不完整时报告依赖缺失；身份冲突时报错。

一个 unit 处理全部八段。父进程在启动计算前验证完整输入合同，并构建每块所需的本块及相邻块路径。已远程核对 production_v2 的 manifest/header，并用新 reader 通过 BCC-CSM2-MR/ssp126/R02C09/solar 的四变量八块预检。实测 manifest 保留旧根路径，按目录身份重定位至当前 blocks 并校验 sidecar；气象块缺单位，新增 `--input-units` 显式声明。详见验证记录。

### 3.3 时间对齐和单位

1. wind 使用 `uas.time`，solar 使用 `rsds.time`，输出保留参考轴的原始数值编码、units、calendar 和分钟偏移。
2. 按真实源时间搜索插值左右邻点，不能用目标的时间位置去切其他变量。当前 final 路径按目标位置加一格 padding 的实现不作为新 reader 的通用保证。
3. block 边界按需读取前后年份块的少量邻点；内部 time chunk 同样读取必要邻点。物理核不需要额外时间窗口或空间 halo。
4. 不外推。完整输入首尾没有插值邻点的位置为缺测；内部数据断档报错，不能跨缺口插值制造数据。输出时次数量以验证过的参考轴为准，不机械要求所有变量都是 `46×365×8`。
5. 支持 Gregorian/proleptic Gregorian 与 noleap/365_day，保留 cftime；跨变量日历不兼容或其他未验证日历明确报错。验证单调、无重复及 3 小时间隔，时间 units 不同可统一换算后比较。
6. 将 rsds 的 W/m²、kW/m²和 tas 的 K、°C 按元数据确定性转换；风速检查 m/s。未知或缺失单位报错，不按某个 tile 的数值中位数猜测，防止分块大小改变结果。历史文件确有缺单位时，先明确其合同再增加显式配置。

## 4. NetCDF 输出约定

建议使用独立的输出根 `outputs/cf_grid/<run-version>/`：

```text
<root>/<model>/<scenario>/<patch>/
    wind.nc                       # 仅 --merge-final 时生成的 wind 46 年文件
    wind.nc.json
    solar.nc                      # 仅 --merge-final 时生成的 solar 46 年文件
    solar.nc.json
<parts-root>/<model>/<scenario>/<patch>/<tech>/<identity>/
    cf_2015-2020.nc[.json]
    cf_2021-2026.nc[.json]
    ...
    cf_2056-2060.nc[.json]
    manifest.json
```

每个 unit 默认发布自身技术的八个年份文件及清单；启用 `--merge-final` 时额外发布一个 final。`parts-root` 默认设为 `<root>/<model>/<scenario>/<patch>/<tech>/blocks/`，允许显式指定独立路径，其下按 identity 隔离。年份文件是持久产物，不按临时文件自动清理。

每份 NC：

```text
time(time)                 参考变量原始 CF 时间编码
lat(lat), lon(lon)          BCSD patch 原生轴，原样保留
domain_mask(lat, lon)       int8，1=计算域，0=域外
wind_cf(time, lat, lon)     float32，或 solar_cf
```

CF 变量的 `units="1"`，`valid_range=[0,1]`，显式设置 float32 `_FillValue=9.96921e36`；读入后可解码为 NaN。零值表示有效输入下的零出力，包括夜间、低风速和切出风速，不表示缺测。

wind 有效性为域内且 uas/vas 有限；solar 还要求 rsds/tas 有限。读取时先正确解码源 fill，单位转换和时间插值后再形成有效性掩膜，物理核计算后还原缺测。现有物理核会将部分非有限输入转为零，必须在写出层覆盖这些位置。光伏夜间输入缺失也写缺测。

首版保留无损 NetCDF4 zlib/shuffle，不量化为整数。chunk 初始候选 `(240,64,64)`，裁剪至实际维度，并与计算 tile 对齐；压缩级别测试 1、2、4 后选择。纯域外 tile 保持未写 fill chunk，避免把海洋显式写零。

属性和 sidecar 至少包含：schema、model/scenario/patch/tech、code SHA及源码摘要、物理参数及风机曲线身份、输入模式、源 manifest/sidecar 身份、land plan 身份、grid/mask/time 指纹、年份范围、实际 time count、单位转换、缺测策略、关键依赖版本和完成状态。风机功率曲线来自 windpowerlib，需记录版本和实际曲线摘要。

年份文件的 `manifest.json` 按合同顺序记录八段身份、路径、时间范围及完成状态，同时记录本次 `merge_final` 参数和独立的 final 状态。默认模式下，八块及其 sidecar 全部验证通过即为 `COMPLETED`，final 状态为 `not_requested`。启用 `--merge-final` 时，八块完成后为 `merge_pending`，必须等 final 和 sidecar 验证通过才将本次 unit 标记为 `COMPLETED`。final sidecar 记录八个源 CF 块身份、总时间数、坐标指纹和源 BCSD 批次。

## 5. 分块计算与性能

### 5.1 作业内完整调用链

```text
一个 unit / 一个作业 / 一次计算入口
    父进程：校验 production-V2 的八段合同、输入身份与网格
        ↓ spawn 进程池，最多八个 worker
    block 00 … block 07：每块读取相应 BCSD block
        → 时间 chunk × 空间 tile
        → 时间对齐、单位转换、有效性掩膜
        → 风光物理核 → 独立 CF 年份文件及 sidecar
        ↓ 全部块完成并验证，释放 worker
    默认：验证八个年份文件及 sidecar → 更新清单 → unit 完成
    --merge-final：父进程按时间及空间小块流式合并
        → wind.nc 或 solar.nc（2015–2060）
        → 校验并发布 final sidecar → 更新清单 → unit 完成
```

八个 block 是计算任务数量，`--processes` 默认值固定为 **8**，可同时处理八块；允许用户显式调小，例如 `--processes 4` 时同一进程池分批处理，仍只有一个作业。恢复时实际 worker 数不超过待计算块数。

风机曲线每 worker 初始化一次，太阳几何的时间项可预计算；不得为优化改变公式或浮点计算顺序而不验证。完整 BCSD final 的小范围读取只放在测试参考中，新入口不增加 final 生产模式。

### 5.2 production-V2 point block 读取

每个年份 worker 打开自身年份的各变量 block，邻年仅用于时间插值。预先建立空间 tile 到 point 位置的映射，将该 tile 气象值填入局部 `(T,Y,X)` 缓冲后调用同一物理核；海洋位置保持无效。输出轴仍是完整 patch 轴。

不直接把 `_read_block_chunk` 的 station gather 改成所有点列表：需要根据实际 point chunk 布局合并连续 point 索引区间，限制缓存容量，避免每个 tile 反复解压同一压缩块。局部散射可在内存中进行，磁盘读取尽量连续；记录实测读放大。若 point 顺序导致重复读取严重，应按输入连续 point 批次读取并在有限 tile 缓冲中组装，不能以提速为由加载整年份天气。

生产直接读取 point blocks，测试中用同源 final 小样本作为空间重建和数值对照。性能测试用于选择 chunk、tile 和 worker 数，不改变输入批次及 blocks 来源。

### 5.3 进程与发布

第一版仅按年份段多进程，worker 数不超过待处理段数。每个 worker 自行打开输入并独占其年份临时 NC；不嵌套空间进程池。一个年份段仍可通过时间与空间 chunk 控制内存，单段重试也可单进程运行。

发布顺序：独有临时 NC → 关闭并验证 → 原子替换正式 NC → 原子发布 sidecar。使用每输出文件锁阻止两个运行同时覆盖；父进程统一更新该 unit 的 manifest；unit 锁防止重复入口同时更新清单。NC 有而 sidecar 缺失的情况视为未完成，不凭文件存在直接跳过。

恢复粒度为年份文件；匹配身份且验证通过的年份复用，不完整年份重算。输入、网格、科学参数或实现变化不得复用旧结果；已有不兼容产物要求新输出根或显式覆盖。计算完成前复查输入身份，避免中途输入变化仍发布成功。

仅启用 `--merge-final` 时执行 final 合并，由父进程按合同顺序逐个打开 CF 年份文件，按 `time_chunk × tile_y × tile_x` 读写，禁止整体加载46年数据或使用会整体物化数组的 concat。合并前核验技术、坐标、掩膜、dtype/fill、时间编码和计算身份；核验时间无重复、无缺口且与八段参考轴拼接完全一致。保持域外 fill 和 chunk 压缩设置，合并阶段不重新计算 CF。

final 写入独有临时文件，关闭后验证 header、时间数和分块数据一致性，再原子发布 NC 与 sidecar。合并中断时保留已完成年份块，下次只重做合并。两种模式都保留八个年份文件，不自动删除。`merge_final` 只控制产物发布，不加入年份块的科学计算身份；后续使用相同参数和输出路径增加 `--merge-final`，应复用已验证的八块，仅补做合并。已有 final 不能替代对八块的验证，默认模式也不要求 final 存在。

### 5.4 体量与内存估计

以下按 noleap、3 小时、300×300、float32 估算，单位为十进制 GB/TB：

| 范围 | CF 数组未压缩逻辑大小 |
|---|---:|
| 单 patch、单技术、46 年 | 48.36 GB |
| 单 patch、单技术、6 年 | 6.31 GB |
| 单 patch、两技术、46 年 | 96.71 GB |
| 47 patch × 4 model × 3 scenario × 两技术、46 年 | 约 54.54 TB |

矩阵估算将每个 patch 都近似为 300×300，未计北极额外行、闰年、元数据和压缩开销。实际磁盘量受域外未写 chunk 和压缩影响，不能沿用二值事件信号的压缩率。默认仅保留八块，单 unit 未压缩 CF 逻辑大小约 48.36 GB、全矩阵约 54.54 TB；启用合并后八个 CF 年份块与 final 同时保留时，单 unit 的未压缩 CF 逻辑大小约为 96.71 GB；全矩阵约 109.09 TB。实际占用以压缩文件实测为准；覆盖已有 final 时还需额外容纳旧 final 和新临时 final。合并增加一次完整 CF 读取和写出，但通过小块复制不需要约48 GB的内存。

240×64×64 的一个 float32 数组约 3.75 MiB，光伏四个原始气象数组合计 15 MiB；这只是下界，还需计入插值邻点、多个浮点中间数组、NetCDF cache 和进程运行时。总内存约为 worker 峰值之和加父进程，不允许直接沿用场站数驱动的 `tb` 公式。

本地先测试 32×32、64×64 tile 与 120/240/480 时间 chunk；在合适组合上比较 1/2/4 worker，必要时再测 8。记录读取、插值、物理核、压缩写出、总耗时、峰值内存和输出大小，不只报告核函数速度。

## 6. 文件级实施顺序

| 阶段 | 建议修改 | 验收点 |
|---|---|---|
| 1：输入与网格描述 | 新建 `grid_cf_io.py`，集中管理坐标、掩膜、时间对齐、单位、输入身份；从场站入口提取确有共用需求的小函数 | 校验 production-V2 manifest、八块及 land plan；测试参考坐标一致 |
| 2：串行网格输出 | 新建 `patchify_grid_cf.py`，加入 point block reader、网格 writer和完整 unit 参数；直接复用 `cf_physics.py` | 小样本从天气到 `(time,lat,lon)` 端到端；科学与缺测测试通过 |
| 3：分块读取优化 | 在已验证 reader 上优化连续 point I/O 和有限 tile 缓冲，补齐邻年边界测试 | blocks 与同源 final 测试参考、跨年边界结果等价 |
| 4：多进程与恢复 | 在新入口增加一层年份 `spawn` 调度、身份校验、输出锁、manifest 和可选流式 final 合并 | 1/2/4/8 worker 一致；两种发布模式、块失败及 merge 失败可恢复 |
| 5：文档与本地基准 | 更新 README 的网格 CLI/格式说明和 `CLAUDE.md` 的入口约定，新增网格测试及验证记录 | 本地真实 CLI 成功；保留既有场站回归；记录实际性能和局限 |

模块先保持上述两个新增文件，不预建通用计算框架。现有场站入口可继续作为回归参考，旧结果保持独立路径。`cf_physics.py` 不为此次空间扩展更换风机、参数或太阳模型；文件中“仅供场站”的说明可在实施时同步更新。

现有 NumPy、pandas、xarray、netCDF4 和 windpowerlib 可完成方案，首版不新增依赖。测试前检查环境是否包含所需测试工具，运行使用 `.venv/bin/python`。

CLI（除生产根目录外，路径及任务参数按实际填写；缺单位的 production_v2 blocks 还需如下显式单位参数）：

```bash
.venv/bin/python patchify_grid_cf.py \
  --bcsd-root /work/share/aczlvkl1ac/bcsd_runs/production_v2 --model <MODEL> --scenario <SSP> \
  --patch <PATCH> --patch-manifest <PATCH_MANIFEST> \
  --land-plan <LAND_PLAN> --tech solar \
  --input-units tas=K rsds=W/m2 uas=m/s vas=m/s \
  --years 2015-2060 \
  --tile-shape 64 64 --time-chunk 240 --processes 8 \
  --compress-level 2 --output-root <GRID_CF_OUTPUT_ROOT>
```

manifest 从 `production-V2` 根目录按验证后的命名规则发现，缺失则报错。`--land-plan` 与 manifest 指向的批次须一致。示例为一个 unit 的一次入口调用，在内部处理八个 block，默认不合并 final；`processes=8` 是默认值，资源不足时可显式调小。

argparse 默认值：

```python
parser.add_argument("--processes", type=int, default=8, help="年份块计算进程数，默认8")
parser.add_argument(
    "--merge-final",
    action="store_true",
    default=False,
    help="将八个年份 CF 文件流式合并为 2015–2060 年完整 NetCDF",
)
```

上述 CLI 不传此参数时只生成八个年份文件；需要完整46年文件时，在同一命令末尾追加 `--merge-final`。两个模式均在一个 unit 的同一作业内完成，参数不改变 BCSD blocks 输入或计算进程结构。

## 7. 测试与验收

新增 `tests/test_grid_cf.py`，复用现有合成数据构建思路，覆盖以下行为：

1. **科学等价**：在每个小网格中心放置测试站，以 nearest 精确命中；有效输入下网格与场站 CF 在 `rtol=1e-6, atol=1e-6` 内一致。涵盖风机切出前后、零风、夜间、日出日落、高纬和跨 180° 经度。旧入口吞掉缺测的位置单独核验，不以旧零值为正确缺测结果。
2. **网格身份**：lat/lon 的值、dtype、顺序与 BCSD 完全相同；海岸、全域外 tile、尾部不足 tile、北极额外行和相邻 patch core 覆盖正确。
3. **时间等价**：整段参考与多 chunk/多年份拼接一致；包括 90 分钟偏移、源变量正负偏移、源首尾缺邻点、跨年和不同变量时次数量。缺口、重复时刻及日历冲突须报错。
4. **日历**：Gregorian 闰年和 noleap 均保留原始数值时间和 units/calendar；太阳几何使用现有 `_doy_hour` 语义。
5. **缺测与单位**：NaN/Inf/fill、海洋、缺单位和错误单位；0 CF 可正常解码为零，solar 夜间无效输入仍缺测，所有有效 CF 均有限且在 `[0,1]`。
6. **输入与参考**：生产 point block reader 与测试用同源 final 在有效域逐时逐点一致；错误 point 身份、重复索引、不同点顺序或缺失 sidecar 不得通过。
7. **并行与恢复**：不同 worker 数、tile 和 time chunk 结果一致；不完整 NC、坏 sidecar、输入变化、并发写同一路径、块重试和 final 合并失败均遵循发布规则；final 分块读取结果与八块按合同拼接完全一致。
8. **实际 CLI**：以小型但时间合同完整的合成输入运行新入口，测试真实 spawn 子进程和参数解析；验证不传 `--processes` 时为8，非正数明确报错；验证不传 `--merge-final` 时为 False、只发布八块且可完成，传入时为 True 并发布完整 final；验证默认运行后追加开关复用八块，以及合并失败后重试。文件 header 与 sidecar、manifest 相互一致。小样本合同可用短年份范围，不依赖 46 年生产数据。

执行 `.venv/bin/python -m pytest tests/ -q` 保留原有测试。当前既有测试仅覆盖场站 gather 裁剪及 final 流式结果，不代表 blocks、网格缺测或恢复已验证。

完成标准是：直接读取 production-V2 八个年份 blocks，一个 unit 内多进程计算后默认发布八个原生网格年份文件，`--merge-final` 可额外发布46年 final，计算与有效性正确，计算及合并内存由 chunk/tile 控制，块和合并均可恢复，端到端测试通过。全球单文件拼接、下游 Loss 改造和超算部署留待后续独立任务。

## 8. 本次核查与参考依据

- 已核对当前 `patchify_station_cf.py`、`cf_physics.py`、两个既有测试、SCNet 生成器及 CF 年份块文档。
- 已核对本地 BCSD `land_plan.py`、`patches.py`、`finalize.py` 和年份常量；实现阶段已补充远程真实 header 和一个完整 unit 的元数据校验，见验证记录；尚未运行全量生产计算。
- 本次任务中运行现有测试：**4 passed，2 warnings**。warning 分别为 NumPy 二进制兼容性提示及分钟时间戳的浮点编码提示；未调整环境或依赖。
- 本地用 `(16,3,4)` 合成气象、包含 `[0,90,180,270]` 经度，比较光伏三维核与逐点调用：最大绝对差为 **0**。这仅验证物理核已支持规则网格，不替代未来 reader/writer 端到端验收。
- 已用本地 Python 核算本文数组体量及原始输入缓冲大小；未声称已验证生产规模性能。

参考文档：

1. 本仓库 [CF 年份块多进程方案](patchify_CF_复用BCSD年份块的多进程并行方案.md)。
2. 相邻仓库 [全网格极端天气事件识别实施方案](../../extreme_event_definitions/document/全网格极端天气事件识别实施方案.md)。
3. 相邻仓库 [全网格实现与本地验证记录](../../extreme_event_definitions/document/全网格实现与本地验证记录.md)。
4. 相邻仓库 [极端事件并行化分析](../../extreme_event_definitions/document/extreme_event_parallelization_qa.md)。

参考方案用于借鉴原生坐标、显式缺测、分块、身份和本地验收；CF 的具体计算和分文件设计以上述本仓库代码及 BCSD 年份合同为依据。相邻仓库链接依赖当前目录布局，不随本仓库 clone 自动提供。
