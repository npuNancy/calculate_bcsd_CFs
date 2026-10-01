# CodingAgent goal 运行提示词

以下提示词用于用户决定正式部署并持续运行时粘贴；仅创建本文件不授权当前会话提交作业。

```text
请创建并持续执行一个 goal：依据 infos/scnet_patchify_stations/goal.md，完成全球网格 CF 到全球场站 CF 的正式全量提取与全局索引发布。不设置 token budget。

先读 AGENTS.md、scnet-parallelize-workflows 技能以及本目录 README.md、运行前准备.md、goal.md、作业分工/作业组合提交顺序.md、accounts.csv、patch_assignment.csv；代码接口依据 document/场站CF提取使用说明.md 和 station_cf_config.example.json。遇到上下文压缩，读取持久台账接续同一运行，不另建RUN_ID或重复控制循环。

本任务授权完成必要的测试、代码提交/推送、14 worker部署、ACL准备、全部脚本生成、Slurm提交、监控、补槽、合规重试和未提交任务重分配。已有授权和有效配额确认继续使用，不重复询问。不要覆盖已有用户修改、删除历史产物或取消其他作业。缺失的真实输入/权限/配额证据先调查，确需我提供时明确缺项，同时继续其他独立准备。

输入是1872的已完成网格CF：/work/home/acjpoxgsdu/cf_grid/cf_grid_v2_20261001_013720。BCSD已经完成抽检，不再做全量BCSD检查，不重新降尺度或从天气计算CF。场站入口：/work/home/acbw9wpn5k/project_climate_patchify/shared_run_20260913/inputs/stations。装机已经确认是snapshot_total年度累计快照；定位三份实际独立上游参考导出，不能使用占位路径或CSV自证。规范Station ssp585对应源SSP5-6.0；Climate和Station独立，处理全部九配对。

范围固定为CANESM5、MPI-ESM1-2-HR、MRI-ESM2-0、BCC-CSM2-MR × 3 Climate × 3 Station × wind/solar × 47 source_patch，共3,384提取unit，2015–2060八个年份段。不要要求我选择GCM，不缩减范围。最近邻及输出语义沿用现有科学入口。

1872/acjpoxgsdu仅汇总，所有结果直接写它的/work/home/acjpoxgsdu/cf_stations/<RUN_ID>，通过ACL让14worker和汇总账号互相访问；prepare、extract、publish全部在worker的Slurm计算节点运行。登录节点只做小JSON、stat、队列、日志检查。准备期按文档检查ACL继承、跨账号原子写入、flock、实际个人可用空间；prepare后再冻结生成资产、复核容量，并固定prepared hash。

正式代码必须是所有worker相同的干净develop-patch-grid完整SHA。远程Git前module load apps/git/2.30.2并git --version。所有14worker各生成完整3,386份脚本（prepare1+extract3384+publish1），逐文件bash -n和hash校验，跨账号manifest及脚本字节一致。脚本路径通过外部环境变量指定。常规重分配使用现有脚本，不临时生成局部任务包。环境文件不可原地修改影响已排队任务。配置用inputs/campaign_config.json，避免与prepare生成的station_cf_config.json冲突。

先运行唯一prepare，验收Slurm、JobID回执、生产prepared和中间关卡后释放extract；全部提取unit验收后才运行唯一publish。实际提取任务归属按CSV初始化，prepare/publish默认1352但可选其他空闲worker。全局阶段没有足够ready任务时如实记录，不能重复提交以填满槽位。

首次立即检查，之后严格按Asia/Shanghai每小时hh:00、hh:15、hh:30、hh:45监控。每轮必须依次“检查—思考—汇总—执行”，打印预计剩余耗时、预计完成时间、估计依据和不确定性；没有可靠样本就说明缺少依据，后续用非空任务分层耗时与实际吞吐更新，纳入排队、重试、长尾和publish时间。长批次提交也必须及时返回边界检查。

每轮通过各worker自身SSH身份查询account-wide squeue/sacct，不用1872看到的跨用户空队列判断空闲。每账号同时最多20个active，包含其他项目及prepare/publish。持续补到20；自己的ready不足就领取其他账号尚未提交且ready的任务。未补满逐账号写原因。重分配保持unit_id、logical_owner、输入输出、prepared和科学编码不变，更新submit_username/assignment_version并完整记录。active、submitting、unknown不抢占。

遵循goal两级锁顺序：共享runtime/.submit.lock，再目标账号$HOME/.bcsd_submit.lock；锁内重读台账、依赖和空槽，先写submitting，sbatch --parsable返回后立即保存JobID。响应丢失记unknown并核实，未查清不得重投。单个unit任何时候最多一个active作业。资源失败和重试按goal分类规则，不盲目重复提交。

每次验收同时核对Slurm COMPLETED/0:0、对应JobID回执、全部固定身份及产物证据；非空任务必须八段完整，合法EMPTY_NO_STATIONS可计成功，有站点但全缺测不算空任务。squeue消失或文件存在不是成功。不得在登录节点重新运行科学审核。

持续原子更新共享runtime/preparation.json、stages.json、units.json及提交/重分配/循环日志，并镜像到本地被忽略的completion_status。progress.md仅保留真实Last checked及四模型×三Climate极简表，每格累计三个Station共282，全部验收后才加✅；不往表内追加阶段说明。index/units.json是不可变科学任务表，不作动态台账。

完成判据：准备关卡和prepare verified；3,384个extract均通过验收或合法空成功；十二格282/282；唯一publish成功生成production全局索引、覆盖及审核报告，年文件数=8×非空unit数；无active/submitting/unknown/待重试或阻塞，汇总账号可访问全部权威文件，保存completion.json并清空next_check_at。只有满足这些条件才能将goal标记完成；提取表全绿但publish未完成时继续运行。
```
