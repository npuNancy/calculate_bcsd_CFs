# production_v2 全格点 CF 运行包

固定运行 **4 GCM × 3 SSP × 2 tech × 47 patch = 1,128 unit**。
四个GCM依次为 CANESM5、MPI-ESM1-2-HR、MRI-ESM2-0、BCC-CSM2-MR，生成时不需要指定GCM。
一个unit是一个Slurm作业，内部默认8个年份worker，默认不合并46年NC。

采用集中共享目录：所有worker直接写乌镇1872（`acjpoxgsdu`）的
`/work/home/acjpoxgsdu/cf_grid/<RUN_ID>/outputs/`，通过ACL授权。
乌镇1872不运行计算作业；14个计算账号以 [accounts.csv](accounts.csv) 为准。
BCSD production_v2 根目录为 `/work/share/aczlvkl1ac/bcsd_runs/production_v2/`，
`outputs/` 保存完整气象文件，CF实际读取同级 `blocks/`，具体文件由manifest解析。
上游已经抽检，不增加全量BCSD检查。

## 文件

| 文件 | 用途 |
|---|---|
| [goal.md](goal.md) | 固定campaign的运行契约、15分钟监控、补槽、重分配及完成规则 |
| [运行前准备.md](运行前准备.md) | HTTPS部署、共享目录ACL、配置、所有worker完整脚本库存 |
| [作业组合提交顺序](作业分工/作业组合提交顺序.md) | 优先级、初始分工及动态重分配 |
| `accounts.csv`、`作业分工/patch_assignment.csv` | 账号角色与47个patch初始归属 |
| `create_jobs.py` | 标准库作业生成器，不SSH、不提交、不读取气象数组 |
| `run_job.py` | 计算节点适配器，核对配置并调用已有CF入口，写完成receipt |
| [正式运行提示词](develop-patch-grid_正式运行提示词.md) | 3800–4000 UTF-16字符，供CodingAgent启动goal |
| `completion_status/progress.md` | 被Git忽略的极简进度表；首次准备从goal中的模板创建 |

## 作业包

```bash
python3 infos/scnet_patchify_grid/create_jobs.py \
  --jobs-dir "$HOME/cf_grid_runs/$CF_RUN_ID/jobs/cf_grid_v2_v1" \
  --code-sha <完整40位部署SHA> --dry-run
```

去掉 `--dry-run` 生成1,128份脚本。**每个worker都生成完整1,128份，共15,792份副本**，
逻辑任务仍只有1,128个；汇总账号不需要脚本。
同SHA、profile和参数在不同账号生成的脚本与manifest逐字一致，生成器拒绝覆盖。
`--patches` 仅用于本地小矩阵测试，正式准备必须省略；四个模型始终全部生成。

默认资源：`wzhctest`、10 CPU、8 worker、24小时，tile64×64、time_chunk240、压缩级别2。
10 CPU沿用既有CF资源起点，约35 GB仅是历史配额估算，实际以分区规则及sacct为准。
`--merge-final` 默认关闭，不影响1,128个unit的数量；启用时在同一作业内合并。

生成脚本完全不含账号home绝对路径，包括climate路径也通过 `CF_CLIMATE_ACTIVATE` 注入。
`#SBATCH` 日志为相对 `logs/`，提交前必须创建目录，使用
`--chdir="$CF_WORK_ROOT" --account=<实际用户名>`，不能在SBATCH指令中写shell变量。
所有代码、配置、输入和输出路径来自外部 `CF_ENV_FILE`；具体配置见准备文档。

每15分钟按上海时间整点、15、30、45分执行“检查—思考—汇总—执行”。
每账号最多20个account-wide active jobs，14账号最多280槽；每轮持续补槽，
空闲账号可领取其他账号尚未提交的unit。详细运行态存共享runtime JSON，progress.md只保留4×3计数表。

## 本地验证

```bash
python3 infos/scnet_patchify_grid/create_jobs.py --help
.venv/bin/python -m pytest tests/test_cf_grid_jobs.py -q
python3 /data6/yanxiaokai/project_climate/bcsd/utils/check_prompt_chars.py \
  infos/scnet_patchify_grid/develop-patch-grid_正式运行提示词.md
```

目录建设和本地生成测试不表示已经部署、设置远程ACL或提交作业。
生产准备及运行由正式提示词驱动，沿用 [计算入口说明](../../README.md)。
