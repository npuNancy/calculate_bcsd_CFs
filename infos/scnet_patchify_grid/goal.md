# production_v2 全格点 CF：运行契约

本文件供正式运行时的CodingAgent执行。单独阅读它不触发提交；用户粘贴正式运行提示词后，创建并持续执行goal，不设token budget。

## 1. 终端目标与计算合同

完成CANESM5、MPI-ESM1-2-HR、MRI-ESM2-0、BCC-CSM2-MR全部四模型、ssp126/ssp245/ssp585、wind/solar、47个patch，共1,128 unit。
稳定身份：`cf-grid-v2/<model>/<scenario>/<tech>/<patch>`；年份固定2015–2060。

| 内容 | 契约 |
|---|---|
| 计算入口 | `run_job.py` → `patchify_grid_cf.compute` |
| 一个Slurm作业 | 一个完整unit；默认8个spawn worker处理八个年份块 |
| 上游 | production_v2 的 blocks 或 final、sidecar、manifest及共享land plan；按unit冻结来源 |
| 默认结果 | 八个 `(time,lat,lon)` CF年份NC及sidecar、unit manifest、Job receipt |
| 可选合并 | `--merge-final` 默认False；开启后父进程在同一作业内合并 |
| 完成计数 | 每model/SSP有47×2=94 unit，每格完成数/94 |
| 恢复 | 同身份年份块复用；坏块重算；合并失败只补合并；不自动删除产物 |

BCSD已抽检通过，直接接受此前提。准备和监控不得另做全量BCSD扫描、全文件hash或逐块数组审核。
计算入口在计算节点核验当前unit输入属于必要合同检查，不是重复全campaign预检。
当前manifest中存在旧根路径，reader会按目录身份定位当前production_v2/blocks或outputs；land plan的旧共享路径可合法复用。

输入选择：wind 检查 uas/vas；solar 检查 rsds/tas/uas/vas。全部年份NC及sidecar存在才选blocks，
否则整个unit选final，不混合两种布局。计算入口默认auto；生产脚本显式固定分类的blocks/final。
完整NC只按时间/空间切片读取；默认8个独立只读进程，每个写一个年份结果。详见[输入来源与并行](输入来源与并行.md)。

本次输入实现、来源或分类发生变化时使用新RUN_ID、release及完整脚本包。旧运行已停止，
保留台账、日志和产物，不恢复旧控制循环，不原地覆盖旧结果。重新启动生产须用户另行授权。

## 2. 角色与持久路径

`accounts.csv`是唯一账号角色表：14 worker；乌镇1872/acjpoxgsdu只汇总。所有结果直接保存其home共享目录，不能默默换成各账号share或软链接汇总。

```text
1872：/work/home/acjpoxgsdu/cf_grid/<RUN_ID>/
    inputs/patch_manifest.json、input_inventory.csv、release.json
    outputs/<model>/<scenario>/<patch>/<tech>/manifest.json
    outputs/<model>/<scenario>/<patch>/<tech>/blocks/<identity>/cf_<years>.nc[.json]
    runtime/units.json、observations.jsonl、submissions.jsonl、reassignments.jsonl
    runtime/preparation.json、cycles.jsonl、receipts/<unit_id>/<JobID>.json
    runtime/.submit.lock
每个worker：$HOME/cf_grid_runs/<RUN_ID>/
    campaign.env、jobs/<profile>/、logs/
每个worker代码：$HOME/project_climate/repos/calculate_bcsd_-cfs
本地：infos/scnet_patchify_grid/completion_status/progress.md（忽略）
```

超算使用Git前执行 `module load apps/git/2.30.2`，再用 `git --version` 检查；生成的作业也执行此步骤。
部署固定同一完整SHA、相同科学环境；远程checkout只读。修改代码先本地测试、commit/push，再远程HTTPS clone或fast-forward更新，活动任务期间不改变它们使用的checkout。
不依赖远程SSH Git密钥，不把token放URL；HTTPS认证缺失时报告具体缺项。不得reset脏目录或覆盖用户改动。

## 3. 准备与资源

逐项执行《运行前准备》并记录preparation.json：代码、配置、ACL探测、容量、跨账号flock、每账号完整作业库存的SHA及bash语法结果。
全部14个worker各有1,128份脚本，计15,792份；正常重分配只选择目标账号预生成脚本，不临时重写脚本。
生成器不需要GCM参数；不得继承参考Extreme包的模型排除、18账号、三阶段或share存储约定。

默认10 CPU、8 worker、24小时，wzhctest，使用共享climate环境；实际内存配额准备时确认。
BLAS/OpenMP线程为1。只有CPU数不低于worker数才生成profile。
记录成功任务Elapsed及可靠峰值RSS，进程型任务的batch MaxRSS不必然等于整个作业峰值，要结合节点/cgroup证据。
空间预估以年份文件实测更新，至少预留20%资源余量；df的PB级容量不能替代1872的账号配额。
不额外建立BCSD扫描作业或独立merge作业。准备通过即按ready队列补槽，首批成功后持续修正资源估计。

## 4. 固定15分钟循环

Asia/Shanghai 的每小时 `hh:00、hh:15、hh:30、hh:45` 必须执行“检查—思考—汇总—执行”，不按耗时改成30/60分钟。
首次启动立即检查并准备，随后对齐下一个钟点；轮次超时立即补一轮并记录迟到原因，不并行运行两个控制循环。
每轮即使无变化也更新时间、进度和剩余耗时估计。

1. **检查**：14账号的account-wide squeue，终态sacct（主Job及batch步骤）、退出码、日志尾部；查询receipt、manifest、sidecar和stat。不在登录节点打开NC批量读取数组。
2. **思考**：按状态分类，更新剩余工作量、ETA、每账号空槽；识别暂态失败、资源问题、确定性错误和可重分配的待提交任务。
3. **汇总**：原子更新units.json、极简progress表；cycles.jsonl记录本轮检查耗时、预计剩余耗时、预计完成时间、依据/不确定性、补槽计划、next_check_at。每次向用户打印“预计剩余耗时约X小时，预计完成时间YYYY-MM-DD HH:MM（Asia/Shanghai）”；样本不足则给有依据区间或明确无法可靠估计，禁止编造精确ETA。
4. **执行**：锁内去重、重计槽位、提交/允许的重试/重分配，立即记JobID；结束后更新台账和本轮行动记录。每个账号持续补足至20个active jobs；全部ready耗尽、其他项目占槽、账号不可达或存在安全阻塞时写明原因，不用重复任务填槽。

ETA用近期成功unit的墙钟耗时按tech、model及patch参考规模分层估计，结合当前运行剩余时间、待提交工作量、实际可用槽位及Slurm排队。
简单下界为剩余作业小时/有效并发槽位；排队、长尾和重试加入区间。核时/秒与墙钟时间不能混用；确认全部完成后ETA为0。

## 5. 槽位、锁与重分配

每账号上限20，包含其他项目的PENDING/RUNNING/CONFIGURING/COMPLETING及所有占配额非终态；不能只按本campaign jobname计数。
查询失败不能当作0。首选空槽最多账号，同槽位轮询；初始归属用于优先调度，不是独占权。

所有提交在1872共享 `runtime/.submit.lock` 内串行完成，并与该worker既有 `$HOME/.bcsd_submit.lock` 协作。
统一加锁顺序为中心锁→账号锁，使用有界等待；参与的其他项目须遵守账号锁，否则记录竞态风险并核对在途提交。
锁内每次执行：重新读台账 → 查全局同unit是否active/submitting/unknown → 重计目标账号全部作业 → 持久化submitting → `sbatch --parsable -A <实际用户名>` → 立即持久化JobID。
提交始终带 `--chdir="$CF_WORK_ROOT" --export=ALL`，先export该账号的CF_ENV_FILE，提前创建logs。

提交前持久化unit_id、attempt、assignment_version、logical_owner、submit_username、脚本路径/hash、profile、代码SHA、release hash、时间和输出路径。
丢失sbatch回执时按jobname、账号、提交窗口与sacct查证，记unknown；查明前不得再次提交。同unit最多一个active Job。

空闲账号优先领取其他账号 **not_submitted且ready** 的unit；保留logical_owner、输入及输出路径，仅更新submit_username等调度字段。
assignment_version递增，reassignments.jsonl记录旧/新账号、脚本、原因、时间及随后JobID。
不得迁移、取消、复制active任务；失联账号的在途任务必须先查明。失败unit只能在确认终态并满足重试规则后重新排队，其attempt也递增。
共享路径与相同代码/输入身份使换账号重试可复用完整年份块；不复制parts，不改变其身份。完成结果不搬迁。

## 6. 失败及profile调整

分类：not_submitted、submitting、active、succeeded、retryable、resource_failure、deterministic_failure、incomplete_output、dependency_blocked、unknown。
本项目不使用“无场站跳过”：全网格unit均有独立完成判据。
TIMEOUT、节点故障、短暂I/O最多自动重试2次（总attempt≤3）。取消原因需查明，不能一律自动重跑。
OOM需依据资源证据增加CPU/内存或降低worker数后再试；仅变更CPU、墙钟和进程数可复用既有年份块。
改变tile/time_chunk/压缩会改变计算身份，应为受影响unit另设版本化运行根并记录权威位置；不得在共享旧路径默认加overwrite。
确定性的输入/ACL/单位/参数/代码错误暂停受影响unit并报告；不全局阻塞其他ready任务，不绕过校验。
新profile在全部14个worker新目录生成同版完整包并核对hash，旧脚本保留；不修改或取消已提交作业。

## 7. 成功与终止

成功须同时满足：Slurm COMPLETED且ExitCode=0:0；对应JobID的COMPLETED receipt；unit manifest为COMPLETED；八个年份NC和sidecar存在、非partial、身份及stat匹配；1872可访问所有产物。
请求合并时还需final及sidecar完成；默认不要求final。receipt的run/unit/user/JobID、code SHA、release hash、profile、input_mode、input_inventory_sha256和manifest身份须与台账一致。
squeue消失、文件存在或旧场站结果均不算成功。NC数组校验由计算入口在计算节点完成，监控仅用轻量证据。

只有全部1,128 unit succeeded，4×3单元格均94/94，且无active/submitting/retryable/incomplete/blocked/unknown，才能结束goal。
无法解决的外部阻塞必须如实报告，不伪造完成。上下文压缩后从台账和receipt恢复，防止重复提交。

## 8. 极简进度表

仅保留时间及下表，不追加阶段表、JobID、ETA或日志；细节统一存runtime JSON。准备期初始化为0/94而非引用示例中的成功数。
时间使用实际检查的Asia/Shanghai ISO时间。即使无变化，每轮原子重写。

```markdown
Last checked: `<实际时间+08:00> Asia/Shanghai`

| model | ssp126 | ssp245 | ssp585 |
|---|---|---|---|
| CANESM5 | 0/94 | 0/94 | 0/94 |
| MPI-ESM1-2-HR | 0/94 | 0/94 | 0/94 |
| MRI-ESM2-0 | 0/94 | 0/94 | 0/94 |
| BCC-CSM2-MR | 0/94 | 0/94 | 0/94 |
```

每格只统计该model/SSP已完成的CF unit；完成94个才标记`94/94 ✅`。准备步骤进度保存在preparation.json，计算进度按此表展示。
