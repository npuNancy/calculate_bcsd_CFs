# 全网格 CF → 全球场站 CF：SCNet 运行入口

本目录服务一个固定生产范围：**4 GCM × 3 Climate × 3 Station × 2 tech × 47 source_patch = 3,384 个提取 unit**，2015–2060 年八个年份段。生成器无需指定 GCM 或情景；始终纳入 CANESM5、MPI-ESM1-2-HR、MRI-ESM2-0、BCC-CSM2-MR 和全部九种 Climate × Station。

输入是乌镇 1872 已完成的网格 CF：

```text
/work/home/acjpoxgsdu/cf_grid/cf_grid_v2_20261001_013720/
```

BCSD 已抽检通过；不执行 BCSD 全量检查、不读取天气重新计算 CF。场站最近邻、快照容量和数据接口沿用[使用说明](../../document/场站CF提取使用说明.md)及[配置示例](../../document/station_cf_config.example.json)。Station ssp585 对应源 SSP5-6.0；装机为年度累计快照，不能跨年累加。

## 文件与工作流

| 文件 | 用途 |
|---|---|
| [goal.md](goal.md) | 运行、监控、补槽、重试、完成判据 |
| [运行前准备.md](运行前准备.md) | 部署、ACL、输入/容量、全账号完整脚本库存与阶段交接 |
| [accounts.csv](accounts.csv) | 1 汇总账号 + 14 worker 的固定角色表 |
| [作业分工/作业组合提交顺序.md](作业分工/作业组合提交顺序.md) | 初始分工和 ready 队列优先级 |
| [作业分工/patch_assignment.csv](作业分工/patch_assignment.csv) | 47 patch 的逻辑归属基线 |
| [create_jobs.py](create_jobs.py) | 标准库生成器，仅生成，不提交 |
| [run_job.py](run_job.py) | 计算节点预检、调用既有入口、写 JobID 回执 |
| [CodingAgent_goal运行提示词.md](CodingAgent_goal运行提示词.md) | 后续授权持续运行 goal 时粘贴使用 |
| `completion_status/progress.md` | 本地忽略的极简进度表；新 checkout 按 goal 模板初始化 |

| 阶段 | 科学入口 | 逻辑作业数 | 依赖 | 完成证据 |
|---|---|---:|---|---|
| prepare | `prepare_station_cf.prepare(config, shared_root, processes=16)` | 1 | 输入配置、ACL、代码、全包及启动空间检查通过 | 生产 prepared、catalog/mapping/任务表、prepare receipt + Slurm 成功 |
| extract | `extract_station_cf.extract_unit(...)` | 3,384 | prepare 验收、资产只读冻结、空间复核通过 | 每 unit 八段和 audit，或合法 EMPTY_NO_STATIONS；receipt + Slurm 成功 |
| publish | `station_cf_publish.publish(prepared, ledger_path)` | 1 | 所有 extract 已验收成功或合法空任务 | 全局索引、覆盖/审核报告、publish receipt + Slurm 成功 |

prepare 和 publish 都在 **worker 的 Slurm 计算节点**执行。1872/acjpoxgsdu 不运行这两个全局阶段，也不运行提取作业；它创建共享根、处理小配置/ACL并汇总轻量状态。

全部结果直接进入新运行根 `/work/home/acjpoxgsdu/cf_stations/<RUN_ID>/`，不分散写各账号 share，也不采用软链接汇总。各 worker 仅保存代码、脚本、外部环境配置和日志。

## 作业包与最小命令

每个 worker 生成 **3,386 份**：prepare 1 + extract 3,384 + publish 1。14 账号共 **47,404 份脚本副本**，其中提取脚本 47,376 份。逻辑正常运行只提交 3,386 个作业，副本不是重复提交任务。

```bash
python3 infos/scnet_patchify_stations/create_jobs.py \
  --jobs-dir /path/outside/checkout/jobs/station_cf_v1 \
  --code-sha <固定40位SHA> \
  --config-sha256 <inputs/campaign_config.json的64位SHA256> \
  --dry-run
```

实际生成去掉 `--dry-run`，在全部 14 worker 各执行一次。默认 wzhctest、prepare 16 CPU，extract/publish 10 CPU，均为24小时；extract 内 8 个 spawn worker、time chunk 240、station chunk 1024、压缩 2。prepare 默认 16 个 spawn 进程建立源索引（`--prepare-processes` 可调），目录和 mapping 顺序生成；publish 默认 8 个线程读取已验收回执（`--publish-processes` 可调，不能超过 `--publish-cpus`），结合 prepared 和冻结映射生成索引与覆盖报告；不遍历、打开或 stat 提取输出文件，不重新计算有效值/缺测值总量，审核摘要记录该验收口径；实际分区计费/内存规则仍需核实。

脚本通过 `SCF_ENV_FILE` 读取账号环境，使用 `source "$SCF_CLIMATE_ACTIVATE" climate` 激活环境。路径、计费账号均不写入 `#SBATCH`；日志是相对 `logs/%x-%j.*`，控制器在提交前创建 logs，并显式传 `--account=<实际worker用户名>`、`--chdir=<worker work root>` 和 `--export=ALL`。

## 进度口径

表头 ssp126/ssp245/ssp585 表示 **Climate**。每格累积该 GCM/Climate 下的 3 Station × 2 tech × 47 patch，即 **282**；合法空 unit 经回执与记账确认后也计入。八段不是八个 unit。

Markdown 仅保存 Last checked 和 4×3 表。环境/ACL/库存等准备步骤写 `runtime/preparation.json`，prepare/extract/publish 阶段写 `runtime/stages.json`，逐任务和调度详情写 `runtime/units.json`、JSONL；本地镜像放被忽略的 completion_status 中。提取表满格但 publish 未完成时，goal 仍未完成。

本目录的创建、生成器测试和本地样例不代表远程已经部署、生成或提交。后续只有用户明确授权运行时，才执行提示词中的提交监控循环。
