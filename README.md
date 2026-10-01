# BCSD 全格点容量因子

`patchify_grid_cf.py` 从 BCSD production_v2 的年份 point blocks 或完整年份气象 NC 计算全格点风电或光伏容量因子，
输出 `wind_cf(time, lat, lon)` 或 `solar_cf(time, lat, lon)` NetCDF。
坐标来自同批次 land plan 的 patch core，保留完整矩形网格；海洋和必要气象输入缺失的位置写缺测。
有效零出力仍为 `0`，CF 为 float32、范围 `[0,1]`。

一个 unit 是 `model × scenario × patch × tech`。一次入口调用在内部用 spawn 进程池处理年份块，
默认 **8进程**；每个 worker 按时间和空间分块，独占自己的年份文件。
默认只发布八个年份文件；显式传入 **`--merge-final`** 才由父进程流式合并完整2015–2060年文件。
本阶段不生成或提交超算作业。

## 运行

使用现有 `.venv/bin/python`，无新增依赖。以下数据路径属于乌镇1866远程文件系统，
本地运行时需替换为本地输入路径。

```bash
.venv/bin/python patchify_grid_cf.py \
  --bcsd-root /work/share/aczlvkl1ac/bcsd_runs/production_v2 \
  --model BCC-CSM2-MR --scenario ssp126 --patch R02C09 \
  --patch-manifest /path/to/patch_manifest.json \
  --land-plan /work/home/acbw9wpn5k/bcsd/global_bcsd/production_v1/shared/land_plan.nc \
  --tech solar --years 2015-2060 \
  --input-units tas=K rsds=W/m2 uas=m/s vas=m/s \
  --input-mode auto --processes 8 --tile-shape 64 64 --time-chunk 240 \
  --compress-level 2 --output-root /path/to/cf_grid/run1
```

风电使用 `--tech wind --input-units uas=m/s vas=m/s`。风光物理公式沿用 `cf_physics.py`：
GE120/2500、100 m 幂律外推、25 m/s切出；Erbs辐射分解、双轴跟踪和温度修正。

`--patch-manifest` 保持既有格式：`{"patches":{"R02C09":{"core_bbox_360":[60,90,30,60]}}}`。
`--land-plan` 必须与该批次 manifest、sidecar 声明一致。production_v2 当前引用共享的旧路径 land plan，
不代表读取旧批次气象 blocks。

输入约定：

- 从 `manifests/<model>__<scenario>__<variable>__<patch>.json` 解析年份合同和文件名；
  `--input-mode auto` 默认检查整个 unit 的所需变量：全部 blocks 及 sidecar 存在则从
  `--bcsd-root/blocks/<model>/<scenario>/<variable>/<patch>/` 读取；任一缺失则整个 unit
  从 `outputs/<model>/<scenario>/<variable>/` 的最终 NC 读取。`blocks` / `final` 可显式固定来源。
  实测 manifest 存在旧根目录绝对路径，reader 保留并验证其目录身份，将根定位到当前批次，
  再要求对应 sidecar 的输出路径指向实际文件。完整 NC 按 manifest 八个年份合同划分时间索引，
  每个 spawn worker 独立只读打开同一输入文件，按时间和空间小块读取，分别写自己的 CF 年份文件。
- 四类气象输入的部分 header 缺少单位。`--input-units VAR=UNIT` 是显式单位声明，
  仅补充缺失元数据；若与文件已有单位冲突则报错。上例沿用既有 CF 的 K、W/m²、m/s 合同，
  若换用其他生产批次，应确认后填写。缺单位且未指定时停止，不按数据值猜测。
- wind 参考 `uas.time`，solar 参考 `rsds.time`；按真实时间邻点线性插值，必要时读取相邻年份块边界。
  不外推，不跨数据缺口插值。支持 Gregorian/proleptic Gregorian 与 noleap/365_day，3小时时间步。
- `--years` 必须等于 manifest 全部合同段的连续范围，生产默认2015–2060；本地夹具可使用较短完整年份合同。
- blocks 的点身份、顺序和有效域必须一致；final 的网格坐标、shape、时间轴及 sidecar 必须匹配。
  已有文件损坏、权限失败或身份冲突不会触发回退；选定来源缺文件则报错。

## 输出与恢复

```text
<output-root>/<model>/<scenario>/<patch>/
    <tech>/manifest.json
    <tech>/blocks/<identity>/
        cf_2015-2020.nc[.json]
        ...
        cf_2056-2060.nc[.json]
        manifest.json
    wind.nc[.json]      # 仅 wind --merge-final
    solar.nc[.json]     # 仅 solar --merge-final
```

`--parts-root` 可将年份文件改存到 `<parts-root>/<model>/<scenario>/<patch>/<tech>/<identity>/`。
年份文件是持久结果，不自动清理。默认八块及sidecar通过验证即完成；请求合并时必须等合并文件验证完成。
在原命令上追加 `--merge-final`，可复用八个已完成文件，仅补做合并。

输入来源（blocks/final）、输入身份、源码摘要、关键依赖、物理参数、网格、单位和分块/编码配置绑定结果身份；
进程数与合并开关不影响年份文件身份。参数或输入变更须换输出根，或显式 `--overwrite` 重建。
文件缺失、sidecar损坏、文件stat/header不匹配时重算相应年份；合并中断保留年份文件供恢复。
使用 unit 和输出文件锁、独有临时文件及原子发布，防止重复运行同时写同一路径。

sidecar记录各阶段耗时、读取数组字节数及进程峰值RSS。读取字节数是返回数组体量，
不是物理磁盘流量；NetCDF压缩解码可能发生读取放大。每个打开的输入气象变量设置64 MiB有界chunk cache，
内存预算需计入各进程、当前及相邻年份句柄和物理核中间量。

## 验证和场站入口

全格点超算部署、ACL共享、14账号分工及15分钟监控见
[SCNet运行包](infos/scnet_patchify_grid/README.md)。生成器默认覆盖四模型的全部1,128个unit，
每个计算账号都预生成完整作业包，结果直接写乌镇1872的home共享目录。

```bash
.venv/bin/python -m pytest tests/ -q
```

测试涵盖公式对照、原生坐标、缺测、时间边界、日历、串并行、恢复和真实CLI。
实现范围和本地测量见 [验证记录](document/全格点CF实现与本地验证记录.md)，设计见
[实施方案](document/全格点CF计算实施方案.md)。

原场站入口 `patchify_station_cf.py` 继续保留，输出 `(time, station)`，现有命令及 SCNet 生成器不变。
网格CF表示格点上的技术出力能力，不使用装机容量。现有场站Loss读取器不能直接读取新格式；
非线性物理核意味着“插值网格CF”通常不等于“插值天气后计算场站CF”。

## 从已完成网格 CF 提取场站 CF

`prepare_station_cf.py` 准备场站快照目录和全局最近邻映射，`extract_station_cf.py` 按
`model × climate_scenario × station_scenario × tech × source_patch` 提取已有 CF，
生产范围为 **3,384 个 unit**。内部最多八个 spawn 年份进程，输出压缩 `(time,station)` NC，
保留源时间/日历和缺测；`station_cf_reader.py` 提供惰性读取及有界迭代。
装机记录按 2030/2040/2050 **累计快照**独立保存，不跨年累加、不乘入 CF。

配置需要独立上游场站 CSV 导出作为快照来源核对。完整 CLI、发布和读取示例见
[场站 CF 提取使用说明](document/场站CF提取使用说明.md)；
科学契约见[实施方案](document/全网格CF提取全球场站CF实施方案.md)。
