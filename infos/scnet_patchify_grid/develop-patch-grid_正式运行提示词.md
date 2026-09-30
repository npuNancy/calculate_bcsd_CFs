请创建并持续执行一个goal，完成当前仓库develop-patch-grid分支的BCSD production_v2全格点CF计算，不设token budget。这条提示词授权准备、部署、生成脚本、提交、监控、允许范围内重试及任务重分配。先读AGENTS.md、infos/scnet_patchify_grid/README.md、goal.md、运行前准备.md、作业分工/作业组合提交顺序.md及两个CSV；遵循这些文件的计算与运行合同，不停在计划或首次提交。

范围固定为CANESM5、MPI-ESM1-2-HR、MRI-ESM2-0、BCC-CSM2-MR四模型，ssp126/ssp245/ssp585，wind/solar，47个patch，共1,128个unit。无需用户指定GCM。一个unit是1 GCM×1 SSP×1 tech×1 Patch，一个Slurm作业内部默认8进程计算八个年份段；最终合并由参数--merge-final控制，默认不合并，不另拆八个时间作业或独立final作业。年份固定2015–2060；输出保持BCSD原生经纬网格，八份(time,lat,lon)容量因子NC是持久结果。

上游在乌镇1866：/work/share/aczlvkl1ac/bcsd_runs/production_v2/outputs/。CF实际读取同级blocks下八个年份段文件，由manifests解析，不读取完整46年气象NC。用户已经完成BCSD抽检，直接接受此前提，不新增全量BCSD扫描、全文件hash或逐块数组审核；仅检查部署所需权限、小配置和代表路径。计算节点入口对当前unit的输入校验正常保留。manifest旧根路径由reader解析到当前production_v2；共享land plan旧路径本身不代表旧气象批次，不擅自改写上游manifest或搬动land plan。

角色以accounts.csv为准。乌镇1872/scnet-wuzhen-1872/acjpoxgsdu只汇总，不运行计算。十四个worker为乌镇1352、1731、1866、1775、1311、1331、1863、1892、1520、1870、1959、1752、1862、1359，Host及username必须从CSV读取。所有结果直接保存到1872的/work/home/acjpoxgsdu/cf_grid/<RUN_ID>/outputs/；采用ACL集中共享，不改成各账号share加软链接。运行状态和回执保存共享runtime，输入配置冻结在inputs。

按运行前准备.md执行：先module load apps/git/2.30.2并git --version；代码经本地测试、commit/push后，用HTTPS在全部worker部署同一完整SHA的干净develop-patch-grid。不覆盖脏目录。活动作业使用的checkout不得pull。确认共享climate环境可激活、依赖一致、分区计费账号和实际CPU内存规则；计算只在Slurm节点执行。生成脚本通过外部CF_ENV_FILE读取路径，环境激活用CF_CLIMATE_ACTIVATE，不在脚本嵌入任何账号home绝对路径；SBATCH日志使用相对logs，提交前创建目录并指定工作目录。

由1872建立新运行根，对必要祖先仅授予穿越，运行目录设置worker访问及默认ACL，inputs冻结后只读。十四账号实际测试创建、原子rename和跨账号读写，1872必须能读写worker产物；跨账号验证同一个flock的互斥。检查home实际配额、可用空间及inode，df文件系统总容量不代表个人配额。按年份产物实测更新空间预测，保留20%余量，默认不额外存46年final。缺权限由对应数据所有者处理，不通过复制输入绕过身份。容量不足报告具体缺项，不擅自改变汇总账号或路径。

冻结RUN_ID、代码SHA、release.json及其hash、47个patch定义、land plan stat和显式输入单位；uas/vas为m/s、tas为K、rsds为W/m2。全部worker各自调用create_jobs.py生成完整1,128份脚本，合计15,792份副本。生成时省略--patches。核对十四份manifest完全一致，每份脚本hash和bash语法均通过。记录准备证据到preparation.json；已有同版完整包校验复用，不覆盖。默认10 CPU、8进程、24小时，tile64×64、time_chunk240、压缩2，实际资源以分区规则为准。

初始化units台账，稳定unit_id为cf-grid-v2/<model>/<scenario>/<tech>/<patch>，保留logical_owner，记录submit_username、attempt、assignment_version、脚本hash/profile、SHA、release hash、JobID和产物路径。初始归属按patch_assignment.csv；按模型、SSP、tech和patch参考规模优先级排队。准备通过即提交，不再添加重复BCSD检查阶段。

首次启动立即检查，此后严格对齐Asia/Shanghai每小时hh:00、hh:15、hh:30、hh:45，依次“检查—思考—汇总—执行”。用真实时钟计算next_check_at，不把15分钟改成每次处理完再睡15分钟。超过轮次边界立即补检查并记录迟到原因，不启动并发控制循环。每轮无变化也更新记录。

检查十四账号account-wide squeue、终态sacct和ExitCode、日志尾部、receipt及manifest/sidecar/stat；登录节点不全量打开NC数组。区分待提交、提交中、运行中、成功、可重试、资源失败、确定性错误、产物不完整、依赖阻塞及unknown。思考剩余工作、瓶颈、失败原因、空槽与重分配。汇总原子更新台账及极简progress，记录本轮耗时、决策、行动、下次检查时间。执行时补槽、重试或重分配，并马上保存提交回执。

每轮必须打印预计剩余耗时和预计完成时间（Asia/Shanghai），按近期成功unit墙钟时间、tech/model/patch规模、运行剩余量、未提交量、有效并发和排队情况估算。剩余作业小时除以有效槽位仅为下界；将长尾、排队和失败计入区间。样本不足明确说明无法可靠估计及缺少何种样本，不编造精确时间；之后每轮持续修正，不混用核时与墙钟。全部成功时剩余耗时为0。

每个worker最多20个active作业，十四账号最多280槽。每轮持续补足各账号到20，计数包括其他项目的全部占配额非终态，不只统计本campaign；查询失败不能当作0。补不满时记录ready耗尽、其他项目占槽、账号不可达或具体阻塞原因，不能用重复unit填槽。每次sbatch前重新核对。

提交统一先取1872共享runtime/.submit.lock，再取目标账号既有$HOME/.bcsd_submit.lock，有界等待。锁内重读台账、查全局同unit状态、重计账号全部作业，先持久化submitting再sbatch --parsable，按实际username指定account，带正确chdir和export。马上记录JobID；回执丢失记unknown，以账号、jobname、提交窗口和sacct查证，未查清不能再提交。同unit最多一个active Job，不依靠文件锁容忍重复投递。若其他项目不协作账号锁，记录竞态并核查在途提交。

有空槽且本账号任务不足时，领取其他账号not_submitted且ready的unit，使用目标账号预生成脚本。只改变submit_username等调度字段，logical_owner及输入输出路径不变；递增assignment_version，在reassignments.jsonl写原/新账号、原因、时间、脚本与JobID。不迁移、取消或复制active任务，不抢失联账号尚未查明的在途任务；完成产物不搬迁，换账号恢复复用同身份年份块。

TIMEOUT、节点故障或短暂I/O确认终态后最多自动重试2次，总attempt不超过3；取消原因必须查清。OOM先根据资源证据调整CPU/内存或降低进程数，不能原参数无限重试。新profile在所有worker生成新版完整包并校验，普通重分配不重新生成脚本。确定性输入、权限、参数或代码错误暂停受影响unit并报告，其他ready继续。改变科学身份或编码配置不得覆盖旧结果；遵循goal中的版本化规则，禁止自动追加overwrite或删除年份文件。

成功必须同时有Slurm COMPLETED且ExitCode=0:0、对应JobID的成功receipt、完成manifest及八份年份NC/sidecar，身份和stat一致且1872可访问。请求merge时还须final成功；默认不要求final。squeue消失、文件存在、旧场站结果都不能当成功，也不存在无场站skip。progress.md须被ignore，仅写真实Last checked和四模型×三SSP表，每格统计wind+solar共94个unit，未运行从0/94开始，满94才加勾。准备和故障细节只写runtime JSON。上下文压缩后从持久台账恢复，避免重复提交；仅当1,128个unit全部成功且没有未解决状态时完成goal，不把阻塞或首次提交当作完成。
