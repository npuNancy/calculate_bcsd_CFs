# 全球场站 CF：运行契约

本文件针对固定全量任务。用户授权正式运行时创建并持续执行 goal，不设 token budget；仅建设或阅读本目录不触发提交。

## 1. 目标、身份与计算合同

完成 CANESM5、MPI-ESM1-2-HR、MRI-ESM2-0、BCC-CSM2-MR；Climate 与 Station 各 ssp126/245/585 的全部九配对；wind/solar；47 个来源 patch，共 **3,384** 个 extract unit，2015–2060 年八段。实际入口由 `run_job.py` 调用既有 prepare/extract/publish API，不运行旧气象到场站 CF 入口。

```text
station-cf-v1/prepare
station-cf-v1/<model>/climate-<C>/station-<S>/<tech>/<source_patch>
station-cf-v1/publish
```

Climate 选择源 CF，Station 选择位置/装机目录，两个情景不得互相替代。Station ssp585 是 SSP5-6.0 场站的规范标签；容量 `snapshot_total`，不累加快照，不乘入 CF。保留全期 CF、原生日历、有效 0 与缺测。最近邻先选全局格点再确定 source_patch；默认不合并年份 NC。

BCSD 已抽检，接受此前提，不新增全量 BCSD 检查、hash 或重新降尺度。计算节点现有 prepare 对 CF 的时间/坐标/domain/sidecar 做必要索引校验，extract 仅读对应 CF，publish 做既有审核；不是再检查气象输入。所有 CSV 全量比较、NC 数组/空间元数据读取、科学审核均在计算节点，不放登录节点。

## 2. 角色、路径与部署

账号以本目录 accounts.csv 为准：14 worker、各最多 20 个 account-wide active job；1872/acjpoxgsdu 只汇总。prepare/publish 也必须由 worker 运行，计入该账号 20 槽。

```text
汇总根 /work/home/acjpoxgsdu/cf_stations/<RUN_ID>/
  inputs/campaign_config.json、release.json（用户输入，先冻结）
  prepared.json、catalogs/、mappings/（prepare 完成后冻结）
  inputs/station_cf_config.json、grid_cf_index.json、input_release.json（prepare 生成后冻结）
  outputs/<model>/climate_<C>/station_<S>/<source_patch>/<tech>/...
  index/units.json（prepare 的不可变任务表，不是调度台账）
  index/authoritative_index.json、coverage_summary.csv、unmapped_stations.csv.gz、validation_summary.json
  runtime/preparation.json、stages.json、units.json
  runtime/submissions.jsonl、reassignments.jsonl、observations.jsonl、cycles.jsonl
  runtime/receipts/<unit_id>/<JobID>.json、.submit.lock
worker $HOME/cf_station_runs/<RUN_ID>/campaign.<profile>.env、jobs/<profile>/、logs/
本地 infos/scnet_patchify_stations/completion_status/progress.md 与 JSON 镜像（忽略）
```

共享根必须物理位于 1872 home，直接写入，不改 share、不开各账号输出副本。不覆盖已完成网格 CF 运行。源码统一完整 SHA、干净 develop-patch-grid；先本地测试/commit/push，再 HTTPS clone/fast-forward。远程执行 Git 前 `module load apps/git/2.30.2` 并 `git --version`。脏目录不 reset；活动任务 checkout 不 pull。HTTPS 凭据缺失报告缺项，不复制 SSH 密钥或把 token 写 URL。

作业环境激活路径在外部配置：`SCF_CLIMATE_ACTIVATE=/work/home/acbpgywfpz/miniconda3/bin/activate`；生成脚本使用变量展开，满足各账号脚本字节一致的要求。生成器无 GCM/情景筛选参数，始终全量生成。

## 3. 阶段依赖与准备关卡

1. 基础准备：按《运行前准备》验证账号、部署、共享环境、源 CF/站点表代表路径、共享 ACL/跨账号读写/flock、个人容量依据、全部完整作业包。每 worker 3,386 份，14 份 manifest 与逐脚本 SHA 相同，所有脚本 bash -n 通过。
2. prepare：选择一个有槽 worker，唯一活动 prepare job。运行快照 CSV 来源比对、六目录、全局 mapping 和 CF 索引。成功后必须同时有 Slurm COMPLETED/0:0、prepare receipt、生产 prepared（sample=False、3,384 唯一身份）与固定 SHA/config。
3. 中间冻结：检查 prepared 的 `storage_estimate`，结合个人可用空间、压缩实测/保守逻辑量和 ≥20% 余量；冻结生成资产，完成跨账号读取检查。将 prepared 文件 SHA256 写入所有后续 profile 环境，确定 `stages.prepare=verified`。这些关卡没通过，extract 仍 dependency_blocked。
4. extract：3,384 个独立 unit，依赖同一 verified prepare。按 ready 优先级补槽，不设模型/Climate/Station 全组完成屏障。每作业默认 8 年份进程，不拆成 27,072 个调度作业。
5. publish：只有全部 extract 的 Slurm、对应回执和产物验收通过才释放唯一 publish。必须运行在 worker 计算节点。publish 验收后才能结束 goal。

准备/全局审核阶段只有一个 ready 作业时，无法填满 280 槽是正常 DAG 限制，记录 `dependency_blocked`/`ready_exhausted`，不得重复提交共享阶段来填槽。

## 4. 固定每 15 分钟监控

首次立即检查，之后严格对齐 Asia/Shanghai `hh:00、hh:15、hh:30、hh:45`，每轮按 **检查—思考—汇总—执行** 顺序。按真实时钟计算 next_check_at；越界立即补检查并记录迟到原因；只有一个持锁控制循环。长提交批次到下个边界先交还监控，不连续提交到错过多个轮次。

- **检查**：分别用每个 worker 自己的 SSH 身份查询 account-wide `squeue`，以及该账号 `sacct` 主 Job/batch/extern 的状态、ExitCode、Elapsed、MaxRSS、取消/失败原因；读取小段日志、JobID receipt、JSON manifest/sidecar/stat。1872 跨用户 squeue 可能因可见性返回空，不能用它代替各自查询。记录查询时间和错误。
- **思考**：分类未知/活动/成功/空成功/失败/阻塞；计算各账号空槽与依赖可用性，选重分配/重试候选，分析 I/O、队列和长尾。查询失败或证据矛盾不视作空槽或成功。
- **汇总**：原子更新阶段与任务 JSON、极简 progress；每轮打印“预计剩余耗时约……，预计完成时间……（Asia/Shanghai）”，写明估计依据、上下界、样本数量与不确定性；无样本则明确无法可靠估计，不编造精确值。cycles.jsonl 保存检查耗时、决策、计划和下轮时间。
- **执行**：每次在两级锁内重查、领取、提交或允许的重试/重分配，马上保存 JobID。对每个 worker 持续补到 **20 个账号全部 active**（包括其他项目）。未满时逐账号写明其他项目占槽、无安全 ready、阶段依赖、不可达或具体阻塞；不能重复 unit 填槽。结束后写实际行动和更新后的台账。

ETA 用近期成功非空 unit 的墙钟耗时，按 wind/solar、model、Climate/Station 的站点数/唯一格点数/选中 chunk 等规模分层，结合运行剩余、未提交量、实际可用并发和排队等待。空 unit 单独建轻任务样本，不能用其短耗时推断密集 patch。共享 I/O 随并发变化需修正吞吐，长尾/重试和 publish 时间加入区间；prepare 未完成时说明缺少何种样本。提取结束但 publish 待运行时剩余耗时不为 0。整个目标成功后为 0。

## 5. 提交、并发与重分配

任一账号最多 20 个占配额非终态，包含 PENDING/RUNNING/CONFIGURING/COMPLETING 等及其他项目；14 worker 上限 280，prepare/publish 也占其中槽位。不要用 CPU 数或八个 Python 进程当作八个 Slurm 作业。

统一锁顺序：1872 `runtime/.submit.lock` → 目标 worker `$HOME/.bcsd_submit.lock`，有界等待。该账号锁与已有网格/其他项目协作；未协作项目可能产生竞态，记录并核对在途提交，不能宣称锁约束了所有外部程序。

锁内每次：

```text
重读 runtime 台账 → 验证 unit 依赖/当前 assignment
→ 排除同 unit active/submitting/unknown → 刷新目标账号全部 squeue
→ 保存 submitting（账号、尝试、脚本/profile/SHA、时间、固定输出路径）
→ sbatch --parsable --account=<实际worker用户名> --chdir=<WORK_ROOT> --export=ALL <预生成脚本>
→ 立即持久化 JobID、提交时间及原始回执
```

先 export 对应不可变 `SCF_ENV_FILE` 并创建 logs；`#SBATCH` 不展开环境变量，脚本内也不写账号特定 --account。Job 名足够区分 Climate/Station；记账查询应使用足够宽的 JobName 字段，不因截断匹配错误。

sbatch 成功响应丢失时记 unknown，通过实际账号、完整 jobname、提交窗口、squeue/sacct 查证，未查清不重投。对所有阶段保证同 unit 最多一个 active Job，文件锁只是输出保护，不代替提交去重。

**重分配规则：**

1. 每轮补槽时，本账号 ready 不足则必须检查其他账号 `not_submitted` 且 ready 的任务；有候选就领取，不只等原 owner 排队。
2. 保留 `logical_owner`、输入、输出、prepared identity 和编码配置，只改 `submit_username` 等调度字段，`assignment_version += 1`。
3. 选目标账号已有同 SHA/config/profile 完整包内的对应脚本；常规重分配不重跑 generator。
4. 在两级锁内再检查任务状态和目标账号配额，记录旧/新账号、原因、时间、assignment_version、原/新脚本路径/hash，随后补 JobID 到 reassignments.jsonl。
5. 不迁移、不取消、不复制 active 任务；submitting/unknown/失联账号在途作业不抢占。确认终态且符合规则的失败任务可另行重试，attempt 单独递增。
6. prepare/publish 的未提交任务也可选择任意 worker，但各只有一个逻辑身份。不因副本数为 14 而提交 14 次。

## 6. 资源、容量和恢复

默认wzhctest，prepare 16 CPU/16进程，extract 10 CPU/8进程，publish 10 CPU/单进程，均为24小时；科学编码固定 time_chunk=240、station_chunk=1024、compress_level=2。CPU 获取内存的具体计费规则准备时确认，不把历史每核额度当永久保证。BLAS/OpenMP 线程固定 1。

容量必须按 **1872 账号/目录实际可用额度**评估。2026-10-01 曾查得共享文件系统约 16.09 PB 空闲，但 quota 返回空，不代表个人有这些空间。已有适用的用户/平台配额确认可继续引用并注明日期与适用账号；有缺项才补查，不重复请求已给出的授权或信息。抽取前结合 prepared 的 raw_cf_bytes、实际压缩样本、临时/失败文件开销与至少 20% 余量，运行中用成功分片字节更新预测。

分类：`not_submitted, submitting, active, succeeded, succeeded_empty, retryable, resource_failure, deterministic_failure, incomplete_output, dependency_blocked, unknown`。

TIMEOUT、节点故障、确认短暂 I/O 最多自动重试两次，总 attempt≤3；CANCELLED 先查原因。OOM 或内存证据不足先调整 CPU/内存或降低 processes，不能原参数盲试；多进程 batch MaxRSS 不必然是整个作业峰值，结合 cgroup/节点证据。确定性参数、代码、ACL、输入/身份错误暂停受影响分支并报告，其他安全 ready 继续。

仅 CPU、walltime、processes 变化可复用已完成科学分片：在 14 worker 全部生成新 profile 的完整包，保存新 manifest/hash，使用新的不可变 profile 环境文件，不改变旧环境/脚本/活动 checkout。改变 time_chunk/station_chunk/compression、映射、输入或代码会改变身份，需新科学运行版本，不能在旧产物上加 overwrite。

程序复用身份一致的完整分片；孤立/不完整文件保留，新尝试目录由提取器创建。失败任务跨账号恢复时仍访问共享原路径，禁止复制年份文件、递归 chmod 旧结果或删除运行历史。

## 7. 成功与终止判据

每阶段成功均要求 **Slurm COMPLETED + ExitCode=0:0 + 同 JobID 的 COMPLETED receipt**，并核对 run_id/unit_id/真实 username/JobID/code SHA/config hash/release hash/prepared hash/profile/pack hash/script hash 与提交台账。squeue 消失、文件存在或旧结果不能代替。

- prepare：生产 prepared 与 3,384 唯一任务、源索引/目录/mapping 完整，已有科学校验通过，随后冻结/ACL/容量关卡完成。
- 非空 extract：manifest COMPLETED、独立 audit 同身份、八个年份 NC 与 sidecar 的 block identity、源 identity、station/time count、stat 一致，1872 可读写；数组验证已由计算节点执行。
- 空 extract：manifest/audit EMPTY_NO_STATIONS，prepared/mapping 证明该 unit station_count=0，blocks=[]；记为 succeeded_empty。缺源、权限错误或全缺测有站点任务不能冒充空 unit。
- publish：全局索引/validation_summary 为 production COMPLETED、3,384 unit，`year_nc_files=8×nonempty_units`（最多27,072），覆盖表及未覆盖清单完备；对应 publish JobID 成功。只读监控不重新执行 publish 的 NC 抽样。

最终：prepare/extract/publish 全部 verified，3,384 个提取 unit 为 succeeded/succeeded_empty，十二格全 282/282，无 active/submitting/unknown/retryable/incomplete/blocked，汇总账号实际可访问所有权威文件，重分配和提交记录可追溯。保存 runtime/completion.json，清空 next_check_at，再将 goal 标记完成。提取表全绿但 publish 未完成不得结束。

## 8. 极简进度和压缩恢复

本地 `completion_status/progress.md` 仅保留真实 Last checked 和下表；列是 Climate，累计三个 Station。只对已满足全部成功证据的 unit 计数；八段全部完成才记一项，合法空 unit 也记一项。每格满 282 才加勾。

```markdown
Last checked: `<实际ISO时间+08:00> Asia/Shanghai`

| model | ssp126 | ssp245 | ssp585 |
|---|---|---|---|
| CANESM5 | 0/282 | 0/282 | 0/282 |
| MPI-ESM1-2-HR | 0/282 | 0/282 | 0/282 |
| MRI-ESM2-0 | 0/282 | 0/282 | 0/282 |
| BCC-CSM2-MR | 0/282 | 0/282 | 0/282 |
```

每轮均原子更新，即使没变化；阶段进度/ETA/JobID/原因不追加到此表。详细状态存在 runtime/preparation.json、stages.json、units.json 和 cycles.jsonl，并镜像到本地忽略目录。

基础准备步骤至少记录：code、environment、inputs_visibility、snapshot_references、ACL、flock、capacity、packs（逐账号）、prepare_job、assets_freeze、extract、publish。每项有 checked_at/status/evidence/error/next_action；阶段台账另含依赖状态。任务台账含 stable unit_id、stage、双情景、logical_owner、submit_username、attempt、assignment_version、profile、脚本与各身份 hash、JobID、状态、证据和原因。

上下文压缩/恢复后先读这些持久记录，确认旧控制器是否存活、锁持有者和全部 in-flight；不以对话遗漏或命令观察超时判定作业失败，不另开并发循环或重建 RUN_ID。仅在实际完成上述终点时关闭 goal。
