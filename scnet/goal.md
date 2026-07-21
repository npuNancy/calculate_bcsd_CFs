# SCNet 多区域 BCSD 容量因子并行计算操作手册

## 1. 总目标

在 SCNet 超算上使用本仓库的：

- `S01E01_Simulate_Solar_CF_BCSD.py`：计算光伏容量因子；
- `S01E02_Simulate_Wind_CF_BCSD.py`：计算陆上风电容量因子；
- `scnet/create_cf_S01_jobs.py`：生成 Slurm 作业脚本；

对已经完成 BCSD 降尺度的多个 `energy × region × scenario` 组合计算容量因子。

本手册只规定通用运行方式，不固化具体超算服务器、账号、区域、CMIP6 模式、情景或年份。每次生产运行的范围必须由命令行参数和对应模式的本地状态目录明确记录。

并行单位固定为：

```text
一个 energy × region × scenario = 一个 Slurm 作业
```

光伏和风电必须拆成不同作业：solar 作业只运行 S01E01，wind 作业只运行 S01E02。不同能源、区域和情景组合之间都可以独立调度。

## 2. 固定约束

以下约束适用于所有服务器、区域、模式和情景：

1. 分区使用 `wzhctest`。
2. 账号在 `squeue` 中的全部作业总数不得超过 20；`PENDING`、`RUNNING` 和过渡状态均计入，不能只统计本项目作业。
3. 超算约按每核 `3.5 GB` 分配内存。S01E* 当前是单进程程序，多申请的核主要用于取得足够内存，不会自然获得等比例计算加速。
4. 不删除或覆盖现有 `~/bcsd/data`、`~/bcsd/outputs`、`~/bcsd/runs` 和 `~/bcsd/cache`。
5. 输入 BCSD 输出必须先校验；不能只因文件存在或上游 Slurm 状态为 `COMPLETED` 就开始计算。
6. 远程项目代码应由本地修改、测试、提交并推送后，再在远程干净工作区执行 fast-forward pull；服务器专用配置和生成结果放在仓库外的 `~/jobs`、`~/logs` 或明确的运行记录目录。
7. 提交流程与 BCSD 共用 `$HOME/.bcsd_submit.lock`，避免两个提交流程并发突破 20 个作业的硬上限。

另外必须遵守：

- 不在登录节点运行完整年份范围的 S01E* 计算，只通过 Slurm 运行；
- 不把 `region=all` 用作生产并行入口，生产范围必须显式列出；
- 不因已有输出文件存在就直接认定完成；
- 不默认使用 `--overwrite`，重算必须精确到单个组合并先记录原因；
- 不把凭据、账号密码或私钥写入仓库、作业脚本和日志。

## 3. 输入、输出与运行参数

### 3.1 输入目录

超算上的 BCSD 降尺度数据根目录固定为：

```text
$HOME/data/bcsd_outputs
```

S01E* 使用的目录结构为：

```text
$HOME/data/bcsd_outputs/<MODEL>/<REGION>/<MODEL>/*.nc
```

不同能源作业的输入依赖为：

```text
solar / S01E01：rsds、tas、uas、vas
wind  / S01E02：uas、vas
```

如果 BCSD 输出实际位于其他目录，应在 `~/data/bcsd_outputs` 下建立软链适配层；不得为适配目录结构复制整套大文件，也不得移动上游 BCSD 输出。

### 3.2 输出目录

默认输出位于远程容量因子项目中：

```text
<PROJECT_DIR>/output/
├── CFs_of_solar/<MODEL>/<REGION>/
└── CFs_of_wind/<MODEL>/<REGION>/
```

单个组合预期生成：

```text
solar_CF_<REGION>_<MODEL>_<SCENARIO>_<YEARS>_allmonths.nc
wind_CF_<REGION>_<MODEL>_<SCENARIO>_<YEARS>_allmonths.nc
```

### 3.3 运行范围

服务器由实际部署目标决定，不写入通用生成器；模式、区域、情景和年份在生成作业时显式传入：

```bash
python scnet/create_cf_S01_jobs.py \
  --cmip6-model '<MODEL>' \
  --regions '<REGION_1>' '<REGION_2>' \
  --scenarios '<SCENARIO_1>' '<SCENARIO_2>' \
  --years '<START_YEAR>-<END_YEAR>'
```

区域名包含空格时必须加引号。传入的区域名还必须与输入目录和文件名完全一致，不能擅自互换空格、连字符或非 ASCII 字符。

### 3.4 `data/maps` 与 `data/grid_of_regions` 依赖

经当前代码检查，`S01E01_Simulate_Solar_CF_BCSD.py` 和 `S01E02_Simulate_Wind_CF_BCSD.py`：

- 不读取 `data/maps`；
- 不读取 `data/grid_of_regions`；
- 只从 `--data_dir` 读取 BCSD NetCDF；
- 通过 Python 包 `global_land_mask` 在运行时构建陆地掩膜，而不是读取仓库内地图文件。

因此，按当前代码部署时，clone 本仓库后不需要额外同步 `data/maps` 或 `data/grid_of_regions`。

每次 S01E* 的数据依赖发生修改后，部署前必须重新检索确认。如果未来代码开始读取 `data/maps`，由于该目录可能未纳入 Git，必须在远程 clone 或 fast-forward pull 后，从本地项目目录同步到目标服务器的对应目录：

```bash
rsync -av --progress \
  ./data/maps/ \
  '<SSH_HOST>:<REMOTE_PROJECT_DIR>/data/maps/'
```

同步后检查远程文件数、大小和读取权限。不得在依赖 `data/maps` 的作业提交之后再补同步，也不得用宽泛 rsync 覆盖远程 `data/` 下的其他目录。

## 4. 本地监控与状态目录

状态记录保存在本地仓库，而不是散落在远程服务器。不同 CMIP6 模式必须使用相互隔离的目录：

```text
./scnet/models/<MODEL>/completion_status/
```

当前已初始化：

```text
./scnet/models/CANESM5/completion_status/
```

这里的 `CANESM5` 只是已经建立的模式级目录，不表示本手册把运行模式固定为 CANESM5。运行其他模式时应创建对应的：

```text
./scnet/models/<OTHER_MODEL>/completion_status/
```

`completion_status/` 是本地运行期状态目录，不由 Git 跟踪。仓库的 `.gitignore` 使用：

```gitignore
scnet/models/*/completion_status/
```

因此无需在目录中保留 `.gitkeep`。首次监控或初始化某个模式时，由程序或操作者用 `mkdir -p` 创建目录。需要迁移或备份运行状态时，应单独归档，不通过源码提交传播。

每个模式的状态目录建议包含：

```text
completion_status/
├── completion_<SERVER>.csv  # 组合级完成状态
├── usage_<SERVER>.csv       # 队列、耗时和核时快照
├── latest_snapshot.json     # 最近一次结构化快照
├── progress_summary.md      # 人工可读进度汇总
└── monitor.log              # 本地监控器运行日志
```

服务器名称只出现在状态文件名和记录内容中，不写死在通用监控逻辑里。

### 4.1 完成状态 CSV

一行对应一个“服务器 + 能源 + 模式 + 区域 + 情景”。solar 和 wind 是两个不同作业，因此分别占一行：

```csv
更新时间,服务器,能源,模式,区域,情景,年份,作业ID,Slurm状态,输出状态,输出路径,运行时长,分配核数,申请内存GB,峰值内存GB,日志路径,重试次数,备注
```

`输出状态` 只使用：

```text
NOT_READY READY SUBMITTED RUNNING VALID MISSING BROKEN FAILED CANCELLED OOM CHECK_FAILED
```

单个能源作业只有在对应输出通过校验时才算完成。Slurm 作业为 `COMPLETED` 不能直接转换为 `VALID`。一个 `region × scenario` 的完整风光任务只有在 solar 行和 wind 行均为 `VALID` 时才算完成。

### 4.2 资源使用 CSV

每次监控追加一行快照：

```csv
检查时间,服务器,模式,统计起始时间,账号活动作业数,本项目等待作业数,本项目运行作业数,本项目完成作业数,本项目失败作业数,累计核时,备注
```

其中“账号活动作业数”统计该账号在 `squeue` 中的全部作业；本项目各状态只统计作业名以 `cf_S01_` 开头且属于当前模式的作业。

## 5. 作业生成器

### 5.1 功能边界

作业生成器为：

```text
scnet/create_cf_S01_jobs.py
```

它负责：

- 使用 argparse 接收模式、区域、情景、年份和资源参数；
- 按 `energy × region × scenario` 生成独立脚本；
- solar 脚本只运行 S01E01，wind 脚本只运行 S01E02；
- 在每个作业中激活指定的 conda 环境；
- 创建作业脚本目录和日志目录；
- 对空值、重复区域、重复情景和非法年份做配置校验；
- 生成安全且唯一的 Slurm 作业名；
- 在脚本成功结束时写出 `[CF_DONE]` 日志标记。

它不负责：

- 调用 `sbatch`；
- 判断输入 NetCDF 是否完整；
- 判断已有输出是否有效；
- 自动删除、覆盖或重提作业；
- 自动维护本地完成状态 CSV。

因此，生成成功不等于可以提交；提交前仍须完成第 7 节检查。

### 5.2 主要参数

| 参数 | 含义 | 默认值 |
|---|---|---|
| `--cmip6-model` | CMIP6 模式 | 必填 |
| `--regions` | 一个或多个区域 | 必填 |
| `--scenarios` | 一个或多个情景 | `ssp126 ssp245 ssp585` |
| `--years` | `YYYY` 或 `YYYY-YYYY` | `2015-2060` |
| `--energies` | `solar`、`wind` 或两者 | `solar wind` |
| `--data-dir` | BCSD 输入根目录 | `~/data/bcsd_outputs` |
| `--project-dir` | 远程项目目录 | 根据生成器位置推断 |
| `--output-root` | 输出根目录 | `<project-dir>/output` |
| `--jobs-dir` | 生成脚本目录 | `~/jobs/cf_S01` |
| `--logs-dir` | Slurm 日志目录 | `~/logs/cf_S01` |
| `--python-executable` | 作业使用的 Python | `python` |
| `--conda-activate` | conda 激活脚本 | `/work/home/acbpgywfpz/miniconda3/bin/activate` |
| `--conda-environment` | conda 环境名 | `climate` |
| `--partition` | Slurm 分区 | `wzhctest` |
| `--nodes` | 节点数 | `1` |
| `--cpus-per-task` | 单任务 CPU 数 | `2` |
| `--walltime` | 可选时间上限 | 不写入 |
| `--chunk-time` | 时间分块长度 | `512` |
| `--compress-level` | NetCDF 压缩等级 | `4` |
| `--overwrite` | 是否覆盖已有输出 | 关闭 |

生产运行必须使用 `wzhctest`。虽然生成器保留 `--partition` 便于测试和兼容，但正式生成后必须检查所有脚本的分区行。

表中的能源、情景和年份默认值只是生成器当前的命令行便利值，不代表生产范围。生产运行必须显式传入 `--cmip6-model`、`--regions`、`--scenarios`、`--years` 和 `--energies`，并把实际范围写入对应模式的本地状态记录。

### 5.3 作业命名

```text
作业名：cf_S01_<ENERGY>_<MODEL>_<REGION>_<SCENARIO>
脚本名：cf_S01_<ENERGY>_<MODEL>_<REGION>_<SCENARIO>.sh
日志名：cf_S01_<ENERGY>_<MODEL>_<REGION>_<SCENARIO>_<JOB_ID>.out
```

模型、区域和情景会转换成安全的 ASCII 作业标识；真正传给 S01E* 的区域参数仍保留原始名称。

## 6. 部署和环境检查

### 6.1 本地修改流程

在本地仓库完成修改和测试：

```bash
git status --short
.venv/bin/python -m py_compile scnet/create_cf_S01_jobs.py
git diff --check
```

确认差异后再提交并推送。不得直接在远程仓库修改 S01E*、生成器或监控器后形成未回传的代码分叉。

### 6.2 远程更新

在目标服务器检查：

```bash
cd '<REMOTE_PROJECT_DIR>'
git branch --show-current
git status --short
```

只有分支正确且工作区为空时才执行：

```bash
git pull --ff-only origin '<BRANCH>'
git rev-parse --short HEAD
```

如果工作区不干净，停止更新并先确认改动归属；不得使用 `git reset --hard` 或覆盖远程文件。

### 6.3 Python 环境

Slurm 作业直接复用已经存在的 `climate` 环境。每个生成脚本必须包含：

```bash
set -eo pipefail
source /work/home/acbpgywfpz/miniconda3/bin/activate climate
set -u
```

必须在激活环境后再开启 `set -u`；该共享环境的 activate hook 会读取可能尚未定义的变量，如果先启用 nounset，作业会在执行 S01E* 前退出。

提交前在目标服务器验证：

```bash
source /work/home/acbpgywfpz/miniconda3/bin/activate climate
python -c \
  "import numpy, xarray, netCDF4, windpowerlib; from global_land_mask import globe; print('环境正常')"
```

环境激活后必须确认 `which python` 和 `python --version`。完整计算不得使用服务器裸 Python 2.7。生成器保留环境路径和环境名参数以便检查或迁移，但当前生产默认值就是上述共享 `climate` 环境。

## 7. 生成与提交前检查

### 7.1 输入完整性

按能源检查每个 `energy × region × scenario` 的输入：

```text
solar：检查 rsds、tas、uas、vas
wind：检查 uas、vas
```

对该能源所需的每个变量检查：

1. 文件存在且大小大于 0；
2. NetCDF 可以打开；
3. 目标变量存在，兼容 `<var>` 和 `<var>_bcsd`；
4. 包含 `time/lat/lon` 维度；
5. 时间范围覆盖请求年份；
6. 时间轴严格递增且没有重复；
7. 首尾时间片均存在有效值；
8. 同一能源作业所需变量的空间网格一致；
9. solar 作业中 `tas/uas/vas` 的时间范围足以插值到 `rsds` 主时间轴。

不满足条件的组合在本地完成状态 CSV 中记为 `NOT_READY` 或 `BROKEN`，不得提交。

### 7.2 生成和静态检查

生成示例：

```bash
python scnet/create_cf_S01_jobs.py \
  --cmip6-model '<MODEL>' \
  --regions '<REGION_1>' '<REGION_2>' \
  --scenarios '<SCENARIO_1>' '<SCENARIO_2>' \
  --years '<YEARS>' \
  --data-dir "$HOME/data/bcsd_outputs" \
  --energies solar wind
```

预期脚本数为：

```text
能源类型数 × 区域数 × 情景数
```

逐个检查本次生成的明确脚本列表：

```bash
bash -n '<JOB_SCRIPT>'
grep -nE '^#SBATCH|source .*/activate climate|S01E01|S01E02|--data_dir|--model|--region|--scenario|--years' \
  '<JOB_SCRIPT>'
```

确认分区、节点、CPU、日志路径、Python、输入目录、输出目录和运行参数正确。不要用无法区分新旧文件的宽泛 glob 直接提交。

### 7.3 精确去重

提交某个组合前依次检查：

1. 当前能源对应的输出是否已经存在且通过第 10 节校验；
2. 当前 `squeue` 是否已有相同精确作业名；
3. `sacct` 是否已有相同作业的成功或失败历史；
4. 历史日志是否已有 `[CF_DONE]`；
5. 本地完成状态 CSV 是否已经记录同一能源组合的有效输出或正在运行的 Job ID。

只有未完成、未在队列且符合重试条件的组合才可提交。S01E* 的“文件存在则跳过”不能替代输出完整性检查。

## 8. 提交流程与 20 作业上限

### 8.1 取得共享锁

所有 CF、BCSD 自动或人工批量提交流程必须使用同一个锁：

```bash
exec 9>"$HOME/.bcsd_submit.lock"
if ! flock -n 9; then
    echo "已有其他提交流程持有锁，本轮禁止提交" >&2
    exit 1
fi
```

持锁期间不得从其他终端绕过锁执行 `sbatch`。

### 8.2 计算空槽

```bash
count_active_jobs() {
    squeue -h -u "$USER" -o '%i' | wc -l
}

active=$(count_active_jobs)
if (( active > 20 )); then
    echo "当前已有 $active 个作业，已超过上限，禁止继续提交" >&2
    exit 1
fi

slots=$((20 - active))
echo "账号全部活动作业: $active；可用空槽: $slots"
```

必须遵守：

- `slots == 0` 时不执行任何 `sbatch`；
- 本轮提交数不得超过 `slots`；
- 其他项目、BCSD、CF 和过渡状态作业都占用槽位；
- 每次 `sbatch` 前后都重新调用 `count_active_jobs`；
- 达到 20 后立即停止，剩余候选留到下一轮；
- 不允许先提交到超过 20，再取消作业降回 20。

### 8.3 显式提交

脚本列表必须来自本轮已经检查和去重的明确清单：

```bash
scripts=(
    '<JOB_SCRIPT_1>'
    '<JOB_SCRIPT_2>'
)

for script in "${scripts[@]}"; do
    [[ -f "$script" ]] || {
        echo "缺少作业脚本: $script" >&2
        exit 1
    }

    active=$(count_active_jobs)
    if (( active >= 20 )); then
        echo "已达到 20 个活动作业，停止提交"
        break
    fi

    job_id=$(sbatch --parsable "$script") || exit 1
    job_id=${job_id%%;*}
    printf '%s,%s\n' "$job_id" "$script"

    active=$(count_active_jobs)
    if (( active > 20 )); then
        echo "错误：提交后活动作业数为 $active" >&2
        exit 1
    fi
done
```

提交后立即把 Job ID、服务器、能源、模型、区域、情景、年份和日志路径写入本地对应模式的 `completion_<SERVER>.csv`。

## 9. 本地监控流程

监控程序在本地运行，通过已配置的 SSH 主机别名只读查询远程服务器。监控器职责是“采集 + 校验 + 报告”，默认不调用 `sbatch`、不删除文件、不自动覆盖输出。

建议每个自然整点执行一次。每轮必须：

1. 读取目标账号的 `squeue`；
2. 使用 `sacct -X` 读取顶层作业状态，避免 `.batch/.extern` 重复计数；
3. 对运行作业用 `sstat` 采集实时资源，对完成作业记录 `MaxRSS`；
4. 扫描 `~/logs/cf_S01/` 中匹配 Job ID 的日志；
5. 解析包含 `energy=<ENERGY>` 的 `[CF_DONE]`；
6. 校验该作业对应能源的 NetCDF；
7. 更新 `completion_<SERVER>.csv` 和 `usage_<SERVER>.csv`；
8. 原子更新 `latest_snapshot.json`；
9. 更新与本轮快照时间一致的 `progress_summary.md`；
10. 把采集成功或失败写入 `monitor.log`。

同一能源组合可能有多个重试 Job ID。监控器应保留实际提供有效输出证据的 Job ID，不能因为更新的重试仍在运行，就把旧 Job ID 已经验证的输出降级。

远程 SSH、Slurm 或 NetCDF 检查失败时记录 `CHECK_FAILED` 或具体错误，不得把旧状态无依据覆盖成 `FAILED`。

## 10. 输出完整性检查

只有同时满足以下条件，单个 solar 或 wind 输出才标记为 `VALID`：

1. 对应 S01 程序成功结束，日志含相应完成标记；
2. 文件存在且大小大于 0；
3. NetCDF 可以打开；
4. 目标变量为 `solar_cf` 或 `wind_cf`；
5. 变量包含 `time/lat/lon` 维度；
6. 时间范围覆盖请求年份；
7. 时间轴严格递增且没有重复；
8. 输出空间坐标与输入一致；
9. 首尾时间片均包含有限值；
10. 抽样的有限容量因子均位于 `[0, 1]`；
11. 海洋格点为 `NaN`，陆地至少存在有效值；
12. 日志中不存在 traceback、OOM、permission denied 等错误。

时间步数应与实际日历和对应输入一致，不为所有模式写死同一个 Gregorian 或 `noleap` 步数。

单个作业完成条件为：

```text
对应能源的输出状态 == VALID
```

完整风光组合完成条件为：

```text
solar 记录 == VALID 且 wind 记录 == VALID
```

## 11. 资源校准与失败恢复

### 11.1 资源校准

首次运行新模式或明显更大网格时，先提交少量组合：

1. 一个小网格组合作功能验证；
2. 一个较大网格组合作资源验证；
3. 无异常后逐步扩大并发度。

初始默认 `--cpus-per-task 2`，约取得 7 GB 内存。由于程序是单进程，多核主要用于内存。若 I/O 竞争导致总吞吐下降，不必机械填满 20 个槽位，应采用实测总吞吐更高的并发度。

### 11.2 OOM

查看真实峰值内存时应检查 `${job_id}.batch`：

```bash
sacct -j '<JOB_ID>' --units=G \
  -o 'JobID,JobName%50,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,MaxVMSize'
```

优先减小 `--chunk-time`，再按约 25% 余量估算 CPU 数：

```text
建议 CPU 数 = ceil(MaxRSS_GB × 1.25 ÷ 3.5)
```

修改资源参数后重新生成并检查目标脚本，只重提失败组合。

### 11.3 其他失败

| 情况 | 处理 |
|---|---|
| 输入缺失或损坏 | 标记 `NOT_READY/BROKEN`，等待上游修复，不提交 |
| Python 依赖错误 | 在当前账号的 Python user-site 补齐依赖，或经授权修复共享 `climate` 环境；完成 import/小测试后重提 |
| 节点故障或抢占 | 日志无代码错误时可按原参数重提 |
| solar 成功、wind 失败 | 保留 solar 作业和输出，只重提独立的 wind 作业 |
| 输出存在但校验失败 | 记录坏文件和原因，获得明确授权后只处理该组合，再用 `--overwrite` |
| 日志长时间不更新 | 先检查进程、I/O 和文件大小是否推进，不能只凭日志安静就取消 |

每次重试前重新精确去重，并在本地 CSV 记录旧 Job ID、失败原因、参数变化、新 Job ID 和重试次数。

## 12. 完成标准

一次生产批次只有同时满足以下条件才完成：

1. 参数指定的所有 `energy × region × scenario` 组合均有记录；
2. 参数指定的每个能源组合均为 `VALID`；如要求风光两类，则每个 `region × scenario` 都同时具备 solar 和 wind 两条 `VALID` 记录；
3. 本地完成状态 CSV、最新快照和进度汇总一致；
4. 所有失败和重试均有可追溯记录；
5. 已汇总作业耗时、CPU、MaxRSS 和输出大小；
6. 没有遗留同组合的重复 `PENDING/RUNNING` 作业。

完成后保留本地模式级状态记录和远程日志，不清理或修改上游 BCSD 数据目录。
