# 全格点 CF 实现与本地验证记录

日期：2026-09-30。分支：`develop-patch-grid`。本次完成计算代码与本地验证，不生成、提交超算作业。

## 已实现

| 文件 | 职责 |
|---|---|
| `grid_cf_io.py` | BCSD manifest/sidecar、年份合同、点身份、land plan、时间与单位校验；连续 point 区间读取；网格 writer、锁、身份与原子发布 |
| `patchify_grid_cf.py` | 全格点风光 CF、时间/空间分块、年份 spawn 进程池、默认8进程、可选 `--merge-final`、状态清单与恢复 |
| `tests/test_grid_cf.py` | 科学参考、网格和时间编码、缺测、生产路径重定位、单位声明、并行、CLI、输入变更和失败恢复 |
| `README.md`、`CLAUDE.md` | 当前入口、生产输入约定、输出和恢复使用说明 |

复用 `cf_physics.py` 原有公式及参数，仅更新其模块用途说明；原场站代码和超算生成器保持不变。无新增依赖。

默认发布八个年份 NC，格式 `(time,lat,lon)`、float32，BCSD patch 原生坐标及参考时间编码保持不变。
八块和sidecar验证通过即可完成；加 `--merge-final` 才生成完整时间段 final，之后增加该开关复用已有八块。
实际worker数不超过待处理块数，进程数变化不使结果身份失效。不同科学输入或分块/编码配置使用新输出根或显式覆盖。

缺测编码与有效的零CF分开。风电要求 uas/vas 有效；光伏要求 rsds/tas/uas/vas 全部有效，包括夜间。
时间插值只使用真实源时间邻点，不外推；源时间精确命中时只使用该点，避免相邻NaN污染。
单位转换在插值后进行，保持与场站入口的浮点运算语义一致。

## production_v2 的实测适配

通过 `scnet-wuzhen-1866` 只读检查目录、manifest、sidecar 和 NetCDF 元数据；未计算生产规模 CF。

- 根目录为 `/work/share/aczlvkl1ac/bcsd_runs/production_v2/`，`blocks/`、`manifests/`、`outputs/` 同级。
- 抽查 manifest 的 `final_tasks[].block_files` 仍带旧 `production_v1` 绝对路径。reader 按严格的
  `blocks/model/scenario/variable/patch/filename` 后缀定位当前根目录下的文件，要求实际sidecar指向当前文件；
  不尝试读取旧根的气象块。源manifest摘要与实际读取文件身份均纳入输出来源。
- 同批次共享 land plan 确实仍位于
  `/work/home/acbw9wpn5k/bcsd/global_bcsd/production_v1/shared/land_plan.nc`。
  输入参数必须与manifest一致，同时检查sidecar记录的路径、文件大小与修改时间。
- 抽查 uas/vas/tas/rsds 块均没有 `units` 属性，原有块上游粗网格样本也未提供单位。
  新增 `--input-units` 显式补充缺失单位。README示例使用既有CF流程采用的 K、W/m²、m/s；
  这不是通过缺失的header独立证明单位。若文件本身提供单位，显式参数不得与之冲突。
- 抽查 point block 为 `(time,point)`，气象数据float32，chunk `(240,4096)`，
  时间为 `days since 2015-01-01`、`365_day`。rsds 相对其他变量偏移90分钟。

用新 `load_plan` 对 **BCC-CSM2-MR / ssp126 / R02C09 / solar** 运行完整元数据预检：

| 检查项 | 结果 |
|---|---|
| 每变量年份块数 | 8，四变量全部通过 |
| patch矩形网格 | 300×300 |
| land plan 有效点 | 89,210 |
| rsds各块时次数量 | 17,520 × 6 + 14,600 × 2 |
| rsds总时间数 | 134,320 |
| 校验范围 | sidecar身份、输入路径、年份覆盖、3小时间隔、跨块连续性、日历、point索引与顺序、掩膜、land plan身份 |

远程预检只在临时目录生成了一个patch边界JSON并自动删除，没有修改生产数据。其他unit未逐一检查；
例如早期抽查的一个MPI目录未发现对应块文件，不能据此声称全部模型/情景已具备完整输入。

## 本地小样本性能

使用 `tests/test_grid_cf.py::build(shape=(32,32))` 构建2024–2025两年noleap合成输入，
共5,840时次、992个有效格点，两个年份文件。含随机天气、NaN/Inf和90分钟时间偏移。
输入chunk为 `(120,128)`；所有试验使用独立输出根，默认不合并。

| 技术 | worker | tile | time chunk | 压缩级别 | 入口总耗时 / s | 两个年份NC合计 / bytes | 最大单进程RSS / MiB |
|---|---:|---|---:|---:|---:|---:|---:|
| wind | 1 | 32×32 | 240 | 2 | 1.77 | 12,978,790 | 345.3 |
| wind | 2 | 32×32 | 240 | 2 | 1.66 | 12,978,790 | 227.7 |
| wind | 2 | 16×16 | 480 | 1 | 1.82 | 13,098,235 | 228.6 |
| solar | 1 | 32×32 | 240 | 2 | 1.82 | 9,494,039 | 345.3 |
| solar | 2 | 32×32 | 240 | 2 | 1.51 | 9,494,039 | 262.1 |
| solar | 2 | 16×16 | 480 | 1 | 2.23 | 9,557,535 | 262.7 |

总耗时包括元数据校验、进程启动和发布。RSS来自进程生命周期高水位；单进程试验与夹具创建共用父进程，
其RSS包含此前创建输入的内存峰值，不能用于比较单worker实际峰值，也不能把表中RSS当作整个作业内存。
逐块读取比较表中三组配置：wind和solar的最大绝对差均为0。
本样本仅两个年份任务，不用于推断8进程性能；8块、2/4/8进程的结果一致性由测试单独验证。
本机未独占，压缩率与时间不能外推到生产数据。用户指定的默认8进程保持不变。

原始结果位于本地被忽略的 `tmp/grid_cf_benchmark/results.json`，运行脚本为 `tmp/benchmark_grid_cf.py`。
上表及参数作为长期记录，不依赖临时产物永久存在。

每个worker的sidecar另记录读取/插值、物理核、写切片及总耗时。写切片时间不包含关闭文件时的全部flush开销，
该开销包括在总耗时内。读取数组字节数不是磁盘实际流量。输入变量有64 MiB有界chunk cache，
当前及邻年句柄的cache、进程数、临时数组需要一起计入生产内存预算。

## 验证范围与限制

本地测试包括完整noleap/Gregorian年份、反向纬度、北极行、0/180/270°经度、尾部tile、域外区域、
正负90分钟偏移、跨年份边界、缺测与零值、旧场站物理核及同源final参考、2/4/8进程一致性。
恢复测试包括新增合并开关、坏sidecar、配置/输入变化、计算期间输入变化、输出锁及合并中断。
真实CLI验证默认8进程、不合并及追加合并路径。

结果身份使用输入stat、元数据文件摘要、源码摘要及运行库版本。恢复检查stat及NetCDF头，
不在每次恢复时完整扫描TB级数据。无法检测被刻意保持大小和mtime的外部数据修改；应将生产输入视为不可变批次。
年份文件保留，计算失败重算对应年份；不提供tile级断点续写。合并按小块复制并对写回数据逐块比较，
不整体加载46年数据。

完整测试：`.venv/bin/python -m pytest tests/ -q`，**28 passed**（最终回归20.43秒）。
7条warning均来自旧场站参考：6条为含Inf的SciPy插值提示，1条为分钟时间戳浮点编码提示；新入口测试无失败。
Python编译检查及 `git diff --check` 通过。
