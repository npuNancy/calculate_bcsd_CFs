# Patchify CF 计算：复用 BCSD 年份块的多进程并行方案

## 1. 结论

CF 的时间维度可以复用 BCSD Step6 的年份块，并在一个 CF unit 内启动 8 个进程：

```text
model × SSP × technology × patch
└── 8 个年份块进程
    ├── block-00
    ├── block-01
    ├── ...
    └── block-07
```

每个进程只读取自己负责的 BCSD 年份块 NC 文件，写出一个 CF 年份块文件；全部子进程完成后，再由父进程执行一次合并，生成现有下游所使用的完整 `wind.nc` 或 `solar.nc` 及其 sidecar。

这里的“合并”是 CF 作业内的本地合并步骤，不是 Slurm 年份分片方案。8 个进程仍属于同一个 CF job，避免把一个逻辑 unit 拆成 8 个独立调度任务。

这种方案比让 8 个进程共同打开完整 2015–2060 NC 更合理：每个进程只访问约八分之一的时间数据，峰值内存和读请求更小，且不同进程读取的文件路径天然互斥。

但是，这不是只增加 `multiprocessing.Pool(8)` 就能完成的改动。必须同时定义：年份块合同、块文件读取器、CF 块输出格式、原子合并规则和最终 sidecar 规则。还需要修改SLURM作业脚本的核数（主要是为了内存分配，3.5GB/核；其次至少需要8核用于8进程）。

## 2. BCSD 年份块合同

当前 `bcsd/global_bcsd` 中的 `APPLY_YEAR_BLOCKS` 为 8 个连续日历年份块：

| block | 年份范围 | 年份数 |
|---|---:|---:|
| 00 | 2015–2020 | 6 |
| 01 | 2021–2026 | 6 |
| 02 | 2027–2032 | 6 |
| 03 | 2033–2038 | 6 |
| 04 | 2039–2044 | 6 |
| 05 | 2045–2050 | 6 |
| 06 | 2051–2055 | 5 |
| 07 | 2056–2060 | 5 |

这 8 块连续覆盖 2015–2060，不重叠、不缺口。CF 不应重新按“每 4/5 年”猜测边界，而应从 BCSD 的三阶段 manifest 读取 `apply_year_blocks` 和每个变量的 block 路径。

### 2.1 远程结果实测

已通过 SSH 在 `scnet-wuzhen-185` 的共享目录做只读检查。当前生产根目录为：

```text
/work/home/acbw9wpn5k/bcsd/global_bcsd/production_v1
```

实测的一个文件为：

```text
blocks/CANESM5/ssp126/rsds/R01C09/rsds_0000000_0017520.nc
```

该文件的 NetCDF 结构为：

```text
dimensions: time=17520, point=36215
data variable: rsds(time, point)
point variables: global_point, y_index, x_index, flat_index, lat, lon
time units: days since 1850-01-01 0:0:0.0
calendar: 365_day
```

因此，当前 `global_bcsd` 远程结果确认是**陆地点 block**，不是 `(time, lat, lon)` 规则网格 block。实测路径使用 `blocks/`，而不是带前导点号的 `.blocks/`；如果其它历史 BCSD 结果目录仍使用 `.blocks/`，CF 必须从对应 manifest 解析路径，不能仅按目录名推断格式。

同一 patch 的 `rsds` block 文件名依次使用以下时间索引范围：

```text
0000000–0017520
0017520–0035040
0035040–0052560
0052560–0070080
0070080–0087600
0087600–0105120
0105120–0119720
0119720–0134320
```

每个 block 旁边有 `block_<var>_<start>_<stop>.json` sidecar，包含 `kind=global-step6-block`、`time_start`、`time_stop`、输入文件 identity 和输出路径；年份范围仍应以三阶段 manifest 的 `apply_year_blocks` 为准。实测的 `rsds` block 单文件约 1.1 GB，说明“每个进程只读一个 block”确实能显著缩小单进程的时间读取范围，但 8 个进程仍会并发读取约 8 个 GB 级文件。

同一年份块的变量时间轴并不一定相同。对同一远程 patch 的第一个块实测为：

| 变量 | 文件名 | `time` 数量 | 首时次 | 末时次 |
|---|---|---:|---:|---:|
| `rsds` | `rsds_0000000_0017520.nc` | 17,520 | 01:30 | 22:30 |
| `tas` | `tas_0000000_0017519.nc` | 17,519 | 03:00 | 21:00 |
| `uas` | `uas_0000000_0017519.nc` | 17,519 | 03:00 | 21:00 |
| `vas` | `vas_0000000_0017519.nc` | 17,519 | 03:00 | 21:00 |

四个文件均为 `calendar=365_day`、`days since 1850-01-01`。因此，`start/stop` 是每个变量自己的时间索引，不能用 `rsds` 的索引位置直接读取 `uas/tas/vas`；CF 必须按实际时间戳对齐。

每个 CF 子进程的输入身份应至少包括：

```text
model
SSP
technology
patch
block_index
start_year
end_year
每个 BCSD 变量的 block 文件路径及 sidecar/hash
```

对于风电，需要 `uas`、`vas`；对于光伏，需要 `rsds`、`tas`、`uas`、`vas`。这些变量共享日历年份块边界，但必须分别验证自己的时间轴、calendar 和时次数量。


每个 CF unit 执行前应该先确定唯一的输入模式：`block` 或 `final_fallback`。两种模式不能在同一个 unit 内按变量或年份混用。

### block 模式的检查顺序

1. 读取 BCSD 三阶段 manifest，解析本 `model/SSP/patch` 对应的 8 个年份块；
2. 对风电检查 `uas/vas`，对光伏检查 `rsds/tas/uas/vas` 的每个 block 文件和 sidecar；
3. 检查 sidecar 的 `kind`、patch、variable、`time_start/time_stop`、输入 identity 和输出路径；
4. 读取每个文件的实际 `time`，检查 `365_day`/calendar、单位、单调性、3 小时步长和年份覆盖；不同变量允许时次数量不同，但必须能按真实时间戳对齐；
5. 检查 block 的 `point` 元数据（或规则网格元数据）与 `land_plan`、空间匹配合同一致；
6. 只有 8 个年份块及其所需变量全部通过检查，才选择 `block` 模式并启动 8 个进程。

### 回退到完整最终文件

如果 block manifest、某个 block 文件或其 sidecar 缺失，但对应的完整最终 BCSD 文件和最终 sidecar 已通过现有完整性校验，则可以选择 `final_fallback`：

- 使用当前完整最终文件读取逻辑完成该 CF unit；
- 在 unit 状态、控制日志、parts manifest 和 CF 最终 sidecar 中记录 `block_fallback=final`；
- 同时记录回退原因、最终文件路径、最终 sidecar hash、代码 SHA 和实际输入模式；
- 该 unit 不再启动 8 个年份块进程，而按完整最终文件的串行/现有流式路径运行。

以下情况不能静默回退：block 文件已经存在，但 `patch_id`、变量名、calendar、时间轴、`global_point/flat_index`、输入 identity 或空间匹配指纹不一致。这些属于确定性输入错误，应停止 unit 并报告；否则会把损坏或串错的 block 隐藏成“成功的完整文件”。如果完整最终文件也不存在或校验失败，则记为依赖缺失/输入不完整，不能开始 CF 计算。

## 3. 推荐的进程结构

### 3.1 父进程职责

父进程不直接计算全部年份，而负责：

1. 读取并验证 BCSD block manifest；
2. 读取场站表、patch manifest 和 station identity；
3. 建立稳定的场站顺序及空间匹配信息；
4. 启动 8 个子进程；
5. 等待 8 个 block 结果；
6. 校验并合并 CF block 文件；
7. 原子发布完整 CF 文件和最终 sidecar。

建议使用 `multiprocessing` 的 `spawn` 上下文。不要在父进程已经打开 xarray/NetCDF 数据集后使用 `fork` 让子进程继承文件句柄；每个子进程应自行打开和关闭自己的 BCSD 文件。

### 3.2 子进程职责

每个子进程只处理一个年份块：

```text
读取 block-XX 的 BCSD 变量
→ 校验该块时间轴和变量身份
→ 进行网格/陆地点到场站的提取
→ 执行风电或光伏 CF 物理核
→ 写出唯一的 CF block 文件
→ 写出 block sidecar
```

子进程不应写最终的 `wind.nc`/`solar.nc`，也不应修改其它 block 的文件。

建议的外部临时布局如下：

```text
cf_parts/<model>/<ssp>/<patch>/<technology>/
├── block_00.nc
├── block_00.nc.json
├── block_01.nc
├── block_01.nc.json
├── ...
├── block_07.nc
├── block_07.nc.json
└── parts_manifest.json
```

这些文件应使用唯一 partial 文件写入，校验成功后原子 rename。某个 block 失败时，只重试该 block，不重算其它 7 个 block。

## 4. 如何复用 `.blocks` 中的 NC 文件

### 4.1 先区分两种 BCSD block 格式

当前仓库中存在两套相关实现，不能仅凭目录名判断格式：

1. `bcsd/step06_apply_sd_3h.py` 的旧/兼容路径使用类似：

   ```text
   <output>/.blocks/<output_stem>/<output_stem>_<start>-<end>.nc
   ```

   通常是 `(time, lat, lon)` 的规则网格 block。

2. `bcsd/global_bcsd` 的三阶段实现使用 manifest 指定的 `blocks/...` 路径，Step6 block 通常是 `(time, point)` 的陆地点文件，并带有 `global_point`、`y_index`、`x_index`、`flat_index`、`lat`、`lon` 等坐标。

因此，CF 方案必须以实际 BCSD manifest 和 block 文件 schema 为准。不能假定所有远程结果都叫 `.blocks`，也不能假定所有 block 都是 `(time, lat, lon)`。

### 4.2 如果 block 是规则网格 `(time, lat, lon)`

这是对现有 `patchify_station_cf.py` 最接近的情况。子进程可以：

1. 读取对应年份 block 的实际 `time` 和规则网格坐标；
2. 用同一个 patch manifest 做场站筛选；
3. 使用固定的空间匹配结果提取场站天气；
4. 对当前年份块计算 CF；
5. 输出 `(time, station)` 的 CF block。

如果 BCSD 变量之间时间轴存在偏移，子进程必须以技术对应的参考时间轴为准：当前实现风电以 `uas` 为参考，光伏以 `rsds` 为参考；其它变量按真实时间戳插值到参考轴。子进程在 block 边界可以只读取相邻 block 的少量边界时次作为插值邻点，然后只保留本 block 的核心时间；不能按文件索引直接配对，也不能把边界插值结果错误地截断成缺测。

### 4.3 如果 block 是陆地点 `(time, point)`

不能直接把它当作现有 CF 入口所需的规则网格。当前 CF 入口的 `find_final()` 和 `_match()` 逻辑面向带 `lat/lon` 网格轴的文件，而 Step6 陆地点 block 使用 `point` 维度。

推荐增加一个“Step6 point-block reader”，而不是把 point block 重新展开成完整全球网格：

1. 从 `land_plan` 和 block 的 `global_point/flat_index` 恢复陆地点身份；
2. 预先建立稳定的 `grid_match_map`，并固定其算法版本；
3. 对每个场站记录其 nearest 或 bilinear 所需的全局网格点索引和权重；
4. 子进程从 block 的 `point` 维度直接 gather 所需的 `global_point`；
5. 对 block 中不存在的 ocean/无效点保留缺测，并沿用当前 finite-weight renormalization 规则；不能为了“找到最近陆地点”而静默改用另一个空间点；
6. 在所有年份 block 中强制使用同一 station 顺序、point 顺序和权重指纹。

当前 CF 默认空间方法为 `nearest`。对该模式，最稳妥的第一版是让 `grid_match_map` 保存规则目标网格的全局索引，再通过 `global_point` 查找该索引是否存在于陆地点 block；如果不存在，应按原规则网格语义得到缺测，而不是改成“最近的可用陆地点”。bilinear 模式还必须保存四个全局点及其权重，并明确 ocean corner 缺失时的有限权重重归一化语义。

这样可以避免读取完整最终网格，也避免每个进程重复构造全球网格。但需要特别验证：point-block 的空间匹配语义必须与原规则网格 CF 结果一致，尤其是海岸附近的 nearest/bilinear 匹配和 ocean 缺测。若无法证明等价，不应直接把 point block 作为生产输入；应先用同一 patch 同时读取最终网格和 point block 做逐站、逐时次对照。

### 4.4 变量时间轴对齐规则

年份块复用的关键不是“八个文件的索引相同”，而是“八个文件覆盖相同的日历年份合同”。每个 CF worker 应执行以下步骤：

1. 从 manifest 得到本 block 的 `start_year/end_year`，分别打开所需变量的对应 block 文件；
2. 解码每个变量自己的 `time`，检查 calendar、单位、单调性和 3 小时规律；
3. 按技术选择参考时间轴：风电用 `uas`，光伏用 `rsds`；
4. 对其它变量按真实时间戳插值或精确对齐，而不是按 `time` 数组位置对齐；
5. 在年份块首尾按需打开相邻 block 的边界记录，提供线性插值邻点；
6. 计算完成后只写参考时间轴在本年份块范围内的核心时次。

相邻 block 的边界读取应限制为必要的 1–2 个时间步，不能为了插值重新读取完整相邻年份块。对于第一个 block 的物理边界，如果没有前置 source 时次，应沿用串行完整文件实现的边界语义，并在 sidecar 中记录产生的缺测/边界处理，而不是擅自外推数据。

## 5. CF block 输出和合并

每个 block 输出建议仍然是 station-only：

```text
dimensions: (time, station)
variables:  wind_cf 或 solar_cf
coordinates: time, station_id, lon, lat, capacity_mw
```

8 个 block 必须共享并验证：

- `station_id` 顺序完全一致；
- `lon`、`lat`、容量和场站数量完全一致；
- `model`、`scenario`、`patch_id`、`tech` 一致；
- calendar、time units、时间步长一致；
- block 年份范围与 manifest 一致；
- 相邻 block 的时间连续且没有重复时刻；
- 每个 block 的输出变量名和 dtype 一致。

合并时按 manifest 顺序连接时间维度，而不是按文件名的字典序盲拼。合并后必须重新生成现有格式的最终文件：

```text
outputs/cf/<model>/<ssp>/<patch>/wind.nc
outputs/cf/<model>/<ssp>/<patch>/wind.nc.json
```

最终 sidecar 应记录：

```json
{
  "status": "COMPLETED",
  "block_count": 8,
  "blocks": [
    {"index": 0, "start_year": 2015, "end_year": 2020, "sha256": "..."}
  ],
  "source_bcsd_manifest_sha256": "...",
  "station_identity_sha256": "...",
  "output_time_count": "sum(block time_count; must match the campaign calendar contract)"
}
```

只有最终合并文件和 sidecar 同时通过校验后，CF unit 才能标记为 `succeeded`。单独存在 8 个 block 文件不能视为 CF 完成。

## 6. 内存、I/O 和进程数建议

### 6.1 I/O 优势

每个进程只读取一个年份块，理论上可以避免 8 个进程共同扫描完整 2015–2060 文件。相较于让所有进程打开同一个完整 NC，优势包括：

- 每个进程的时间范围明确；
- 不需要对完整时间轴做一次性索引；
- NetCDF chunk 读取更容易与年份块边界对齐；
- 失败重试只涉及一个年份块；
- block 文件可复用 BCSD Step6 已完成的中间结果。

但 8 个进程仍可能同时读取同一组变量的 8 个文件。如果远程 `.blocks` 位于共享存储，瓶颈可能从单进程内存转移到共享文件系统总带宽。因此要以实测读吞吐和 `AveCPU` 判断，而不能假定 8 倍加速。

### 6.2 内存边界

现有 CF 流式逻辑按时间块限制临时数组。8 个进程后，总内存约为：

```text
8 × 单进程 block 临时数组
+ 8 × xarray/NetCDF 解码开销
+ 父进程和 merge 缓冲
```

每个子进程都应保持时间块流式处理，不得因为已经按年份分片就把整个年份 block `.values` 一次性加载。尤其是密集 patch 和光伏的四变量输入，必须继续保留内部时间 block。

性能和等价性测试固定使用以下测试账号：

```text
账号名称：乌镇199
SSH：scnet-wuzhen-199
用户名：aclym5felp
```

测试必须选用 `BCC-CSM2-MR` 的一个任务 unit，不使用 CANESM5、MPI-ESM1-2-HR 或 MRI-ESM2-0，以免影响当前正在运行或优先调度的作业。推荐选择一个场站数量较多、8 个 BCSD 年份 block 均已通过检查的 BCC unit，例如：

```text
model=BCC-CSM2-MR
SSP=ssp126
technology=solar
patch=R03C09
```

`R03C09` 只是密集 patch 的候选；开始测试前必须通过 `scnet-wuzhen-199` 检查用户 `aclym5felp` 的 account-wide `squeue`，并检查 campaign `task_state` 和提交 manifest，确认该 unit 没有 `RUNNING/PENDING/CONFIGURING/COMPLETING` 或即将提交的 active job。测试提交时必须使用乌镇199对应的有效计费账号参数，不能沿用10个正式 worker 的用户名或计费账号。如果候选 unit 已被控制器占用，应在 BCC-CSM2-MR 内选择另一个满足同样条件的 unit，不能抢占、取消或覆盖已有作业。

测试必须使用乌镇199上 checkout 外的版本化 pilot 输出、日志和运行状态目录，目录所有者为 `aclym5felp`；不能写入正式 `outputs/cf`，也不能覆盖现有最终文件。BCSD block 和 manifest 只读使用，开始前应确认 `aclym5felp` 对共享输入有读取权限、对 pilot 目录有创建文件和原子 rename 权限。block 加速测试应强制要求 `input_mode=block`；如果 block 不完整，不得自动回退到完整最终文件后把该结果当成 block 加速数据。回退路径可以另做独立的兼容性测试，但不能混入性能对照。

在这个未占用的 BCC unit 上测试 `n=1、2、4、8`。只有在 `n=8` 满足以下条件时才进入正式方案：

- walltime 相对 `n=1` 有稳定下降；
- `MaxRSS` 保留至少 20% 余量；
- CPU 利用率确实提高；
- 没有共享存储错误或 NetCDF 锁冲突；
- 合并结果与串行结果逐点一致或在预先定义的浮点容差内一致。

子进程内部应设置：

```text
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
```

避免每个 Python 子进程再启动一组 BLAS 线程而发生过度订阅。

## 7. 失败、重试和一致性规则

建议将 CF unit 状态细化为：

```text
block_pending
block_active
block_succeeded
block_retryable
block_deterministic_failure
merge_pending
merge_failed
succeeded
```

重试规则：

- 某个 block 的临时 I/O、节点故障或超时：只重试该 block；
- OOM：先根据该 block 的 `MaxRSS` 调整资源，再重试；
- 输入 identity、calendar、时间轴、站场匹配错误：停止整个 CF unit；
- merge 发现缺块、重复时间或 station identity 不一致：不得发布最终 CF 文件；
- 旧 block 文件和旧 merge 文件保留用于审计，不覆盖成新的不兼容格式。

## 8. 与当前 CF/下游依赖的关系

CF 年份块并行不会改变逻辑依赖：

```text
BCSD 8 个年份块
        ↓
CF 8 个年份块进程
        ↓
CF merge / canonical CF
        ↓
Loss 使用完整 CF
```

Extreme 不消费 CF，因此可以继续独立运行；Loss 必须等待完整 CF merge 成功，不能只看到部分 CF 年份块就提交。

当前已有的下游验证器和任务状态通常按一个 `cf/.../patch/tech` 最终文件判断成功。因此引入年份块后，必须把“8 个 block 完成 + merge 成功”作为同一个逻辑 CF unit 的完成合同，不能让监控器把某个 block 误记成完整 CF 成功。

## 9. 推荐实施顺序

1. 从 BCSD manifest 读取 8 个年份块及其输入 sidecar，不扫描或猜测目录；
2. 先实现只读的 block schema 检查；
3. 在乌镇199（`scnet-wuzhen-199`，用户 `aclym5felp`）固定选择一个未占用的 `BCC-CSM2-MR` unit，在隔离 pilot 目录中做 `n=1` 串行 block 基准；
4. 继续在乌镇199使用同一个 BCC unit、相同资源口径和同一个 block reader 测试 `n=2、4、8`；
5. 对每个 block 输出做逐点等价性检查；
6. 实现原子 merge 和最终 sidecar；
7. 通过 merge 后的 canonical CF 文件接入现有 Extreme/Loss 依赖；
8. 最后再评估是否需要将 8 进程固定为默认值。

不建议直接在当前已经运行的单 job campaign 中切换输入格式或输出合同。应先用版本化的实验目录完成 pilot，确认 `.blocks` 实际 schema、空间匹配等价性、MaxRSS 和读写吞吐后，再决定是否用于正式 campaign。

## 10. 最终判断

“复用 BCSD Step6 的 8 个年份块，每个进程只读对应 NC 文件”是一个**方向正确且比完整 NC 多进程读取更好的 CF 并行方案**，特别适合减少单进程的时间范围和重复 I/O。

但需要注意两个前置条件：

1. 远程 `.blocks` 的实际文件格式必须与 CF reader 对接；如果是 `(time, point)`，需要 point-block reader 和稳定的空间匹配合同，不能直接套用当前规则网格 reader；
2. 8 个年份结果必须先合并并验证，才能发布现有下游使用的完整 CF 文件。

因此建议把它定义为：

```text
一个 CF logical unit
= 8 个年份块 worker
  + 一个本地 merge
  + 一个 canonical CF sidecar
```

这是当前最值得优先验证的 CF 多进程方向；不需要把任务进一步拆成 8 个独立 Slurm 分片任务。
