#!/usr/bin/env python3
"""为 S01 BCSD 风光容量因子计算生成 Slurm 作业脚本。

每个 ``energy × region × scenario`` 组合生成一个独立作业，光伏和风电
不会放进同一个 Slurm 作业。本脚本只生成作业文件，不会调用 ``sbatch``。
"""

from __future__ import annotations

import argparse
import os
import re
import shlex
import unicodedata
from pathlib import Path
from typing import Sequence


DEFAULT_DATA_DIR = "~/data/bcsd_outputs"
DEFAULT_JOBS_DIR = "~/jobs/cf_S01"
DEFAULT_LOGS_DIR = "~/logs/cf_S01"
DEFAULT_SCENARIOS = ["ssp126", "ssp245", "ssp585"]
DEFAULT_ENERGIES = ["solar", "wind"]
DEFAULT_CONDA_ACTIVATE = "/work/home/acbpgywfpz/miniconda3/bin/activate"
DEFAULT_CONDA_ENVIRONMENT = "climate"
SOLAR_SCRIPT = "S01E01_Simulate_Solar_CF_BCSD.py"
WIND_SCRIPT = "S01E02_Simulate_Wind_CF_BCSD.py"


def positive_int(value: str) -> int:
    """解析正整数参数。"""
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("必须是正整数")
    return parsed


def compress_level(value: str) -> int:
    """解析 NetCDF 压缩等级。"""
    parsed = int(value)
    if not 0 <= parsed <= 9:
        raise argparse.ArgumentTypeError("压缩等级必须在 0-9 之间")
    return parsed


def years_range(value: str) -> str:
    """校验 YYYY 或 YYYY-YYYY 年份参数。"""
    match = re.fullmatch(r"(\d{4})(?:-(\d{4}))?", value)
    if match is None:
        raise argparse.ArgumentTypeError("年份必须是 YYYY 或 YYYY-YYYY")
    start = int(match.group(1))
    end = int(match.group(2) or match.group(1))
    if end < start:
        raise argparse.ArgumentTypeError("结束年份不能早于开始年份")
    return value


def safe_token(value: str, label: str) -> str:
    """检查用于模型和情景的安全标识。"""
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
        raise ValueError(
            f"{label} 只能包含英文字母、数字、点、下划线和连字符：{value!r}"
        )
    return value


def unique_values(values: Sequence[str], label: str) -> list[str]:
    """保留顺序并拒绝空值或重复值。"""
    result: list[str] = []
    seen: set[str] = set()
    for raw_value in values:
        value = raw_value.strip()
        if not value:
            raise ValueError(f"{label} 不能包含空值")
        if value in seen:
            raise ValueError(f"{label} 存在重复值：{value}")
        seen.add(value)
        result.append(value)
    if not result:
        raise ValueError(f"{label} 不能为空")
    return result


def job_token(value: str) -> str:
    """把模型、区域或情景转换为适合 Slurm 作业名的 ASCII 标识。"""
    normalized = unicodedata.normalize("NFKD", value)
    ascii_value = normalized.encode("ascii", "ignore").decode("ascii")
    token = re.sub(r"[^A-Za-z0-9_.-]+", "-", ascii_value).strip("-._")
    if not token:
        raise ValueError(f"无法从名称生成安全作业标识：{value!r}")
    return token


def expanded_path(value: str) -> Path:
    """展开用户目录并返回绝对路径，但不解引用符号链接。"""
    return Path(os.path.abspath(os.path.expanduser(value)))


def shell_command(value: str) -> str:
    """展开带路径的 Python 命令；普通命令名原样保留。"""
    if "/" in value or value.startswith("~"):
        return str(expanded_path(value))
    return value


def build_parser() -> argparse.ArgumentParser:
    """创建命令行解析器。"""
    repository_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description=(
            "按 energy × region × scenario 分别生成 S01E01/S01E02 Slurm 作业；"
            "只生成脚本，不提交作业。"
        )
    )
    parser.add_argument(
        "--cmip6-model",
        required=True,
        help="CMIP6 模式名称，例如 CANESM5 或 MIROC-ES2H",
    )
    parser.add_argument(
        "--regions",
        nargs="+",
        required=True,
        metavar="REGION",
        help="一个或多个区域；含空格的名称需要加引号",
    )
    parser.add_argument(
        "--scenarios",
        nargs="+",
        default=DEFAULT_SCENARIOS,
        metavar="SCENARIO",
        help="一个或多个情景（默认：ssp126 ssp245 ssp585）",
    )
    parser.add_argument(
        "--years",
        type=years_range,
        default="2015-2060",
        help="计算年份 YYYY 或 YYYY-YYYY（默认：2015-2060）",
    )
    parser.add_argument(
        "--energies",
        nargs="+",
        choices=DEFAULT_ENERGIES,
        default=DEFAULT_ENERGIES,
        metavar="ENERGY",
        help="需要生成的能源类型（默认：solar wind）",
    )
    parser.add_argument(
        "--data-dir",
        default=DEFAULT_DATA_DIR,
        help=f"BCSD 输入根目录（默认：{DEFAULT_DATA_DIR}）",
    )
    parser.add_argument(
        "--project-dir",
        default=str(repository_root),
        help="calculate_bcsd_-cfs 仓库目录（默认：根据本脚本位置推断）",
    )
    parser.add_argument(
        "--output-root",
        default=None,
        help="输出根目录（默认：<project-dir>/output）",
    )
    parser.add_argument(
        "--jobs-dir",
        default=DEFAULT_JOBS_DIR,
        help=f"Slurm 脚本目录（默认：{DEFAULT_JOBS_DIR}）",
    )
    parser.add_argument(
        "--logs-dir",
        default=DEFAULT_LOGS_DIR,
        help=f"Slurm 日志目录（默认：{DEFAULT_LOGS_DIR}）",
    )
    parser.add_argument(
        "--python-executable",
        default="python",
        help="激活 conda 环境后使用的 Python 命令或绝对路径（默认：python）",
    )
    parser.add_argument(
        "--conda-activate",
        default=DEFAULT_CONDA_ACTIVATE,
        help=f"conda activate 脚本（默认：{DEFAULT_CONDA_ACTIVATE}）",
    )
    parser.add_argument(
        "--conda-environment",
        default=DEFAULT_CONDA_ENVIRONMENT,
        help=f"conda 环境名（默认：{DEFAULT_CONDA_ENVIRONMENT}）",
    )
    parser.add_argument(
        "--partition",
        default="wzhctest",
        help="Slurm 分区（默认：wzhctest）",
    )
    parser.add_argument(
        "--nodes",
        type=positive_int,
        default=1,
        help="每个作业的节点数（默认：1）",
    )
    parser.add_argument(
        "--cpus-per-task",
        type=positive_int,
        default=2,
        help="每个作业的 CPU 数；可用于获得按核分配的内存（默认：2）",
    )
    parser.add_argument(
        "--walltime",
        default=None,
        help="可选 Slurm 时间上限，例如 2-00:00:00；默认不写 --time",
    )
    parser.add_argument(
        "--chunk-time",
        type=positive_int,
        default=512,
        help="S01E01/S01E02 的时间分块长度（默认：512）",
    )
    parser.add_argument(
        "--compress-level",
        type=compress_level,
        default=4,
        help="NetCDF 压缩等级 0-9（默认：4）",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="在生成的 S01 命令中加入 --overwrite；默认保留并跳过已有输出",
    )
    return parser


def render_job_script(
    *,
    project_dir: Path,
    data_dir: Path,
    output_root: Path,
    logs_dir: Path,
    python_executable: str,
    conda_activate: Path,
    conda_environment: str,
    model: str,
    region: str,
    scenario: str,
    energy: str,
    years: str,
    partition: str,
    nodes: int,
    cpus_per_task: int,
    walltime: str | None,
    chunk_time: int,
    compression: int,
    overwrite: bool,
    job_name: str,
) -> str:
    """渲染一个 energy × region × scenario 对应的 Slurm 作业。"""
    q = shlex.quote
    slurm_lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --nodes={nodes}",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={cpus_per_task}",
        f"#SBATCH --output={logs_dir}/{job_name}_%j.out",
        f"#SBATCH --error={logs_dir}/{job_name}_%j.out",
    ]
    if walltime:
        slurm_lines.append(f"#SBATCH --time={walltime}")

    common_arguments = [
        "--data_dir",
        str(data_dir),
        "--model",
        model,
        "--region",
        region,
        "--scenario",
        scenario,
        "--years",
        years,
        "--chunk_time",
        str(chunk_time),
        "--compress_level",
        str(compression),
    ]
    if overwrite:
        common_arguments.append("--overwrite")

    if energy == "solar":
        s01_script = project_dir / SOLAR_SCRIPT
        output_dir = output_root / "CFs_of_solar"
    elif energy == "wind":
        s01_script = project_dir / WIND_SCRIPT
        output_dir = output_root / "CFs_of_wind"
    else:  # pragma: no cover - create_jobs 已限制可选值
        raise ValueError(f"不支持的能源类型：{energy}")

    s01_command = [
        python_executable,
        str(s01_script),
        *common_arguments,
        "--output_dir",
        str(output_dir),
    ]

    body = [
        "",
        "set -eo pipefail",
        "",
        f"PROJECT_DIR={q(str(project_dir))}",
        f"REGION={q(region)}",
        f"SCENARIO={q(scenario)}",
        f"ENERGY={q(energy)}",
        "",
        f"source {q(str(conda_activate))} {q(conda_environment)}",
        "set -u",
        f"{q(python_executable)} -c \"import numpy, xarray, netCDF4, windpowerlib; from global_land_mask import globe; print('Python 环境正常')\"",
        "",
        f"mkdir -p {q(str(output_dir))}",
        'cd "${PROJECT_DIR}"',
        "",
        'echo "[INFO] 开始时间: $(date -Iseconds)"',
        'echo "[INFO] 节点: ${SLURMD_NODENAME:-$(hostname)}"',
        f'echo "[INFO] energy=${{ENERGY}} model={model} region=${{REGION}} scenario=${{SCENARIO}} years={years}"',
        f'echo "[INFO] Python: {q(python_executable)}"',
        "",
        'echo "[INFO] 开始计算 ${ENERGY} 容量因子"',
        shlex.join(s01_command),
        f'echo "[CF_DONE] energy=${{ENERGY}} model={model} region=${{REGION}} scenario=${{SCENARIO}}"',
        'echo "[INFO] 完成时间: $(date -Iseconds)"',
        "",
    ]
    return "\n".join([*slurm_lines, *body])


def write_script(path: Path, content: str) -> None:
    """原子写入作业脚本并添加可执行权限。"""
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(content, encoding="utf-8")
        temporary.chmod(0o755)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def create_jobs(args: argparse.Namespace) -> list[Path]:
    """根据命令行配置生成全部 Slurm 作业脚本。"""
    model = safe_token(args.cmip6_model.strip(), "CMIP6 model")
    regions = unique_values(args.regions, "regions")
    scenarios = unique_values(args.scenarios, "scenarios")
    energies = unique_values(args.energies, "energies")
    for scenario in scenarios:
        safe_token(scenario, "scenario")

    project_dir = expanded_path(args.project_dir)
    data_dir = expanded_path(args.data_dir)
    output_root = (
        expanded_path(args.output_root)
        if args.output_root
        else project_dir / "output"
    )
    jobs_dir = expanded_path(args.jobs_dir)
    logs_dir = expanded_path(args.logs_dir)
    python_executable = shell_command(args.python_executable)
    conda_activate = expanded_path(args.conda_activate)
    conda_environment = args.conda_environment.strip()
    if not conda_environment:
        raise ValueError("conda environment 不能为空")

    required_scripts = {
        "solar": project_dir / SOLAR_SCRIPT,
        "wind": project_dir / WIND_SCRIPT,
    }
    missing_scripts = [
        required_scripts[energy]
        for energy in energies
        if not required_scripts[energy].is_file()
    ]
    if missing_scripts:
        missing = ", ".join(str(path) for path in missing_scripts)
        raise FileNotFoundError(f"项目目录缺少 S01 脚本：{missing}")

    jobs_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    created: list[Path] = []
    used_job_names: dict[str, tuple[str, str, str]] = {}
    for region in regions:
        for scenario in scenarios:
            for energy in energies:
                job_name = "cf_S01_{}_{}_{}_{}".format(
                    job_token(energy),
                    job_token(model),
                    job_token(region),
                    job_token(scenario),
                )
                previous = used_job_names.get(job_name)
                if previous is not None:
                    raise ValueError(
                        "不同组合生成了相同作业名："
                        f"{previous} 和 {(energy, region, scenario)} -> {job_name}"
                    )
                used_job_names[job_name] = (energy, region, scenario)

                script_path = jobs_dir / f"{job_name}.sh"
                content = render_job_script(
                    project_dir=project_dir,
                    data_dir=data_dir,
                    output_root=output_root,
                    logs_dir=logs_dir,
                    python_executable=python_executable,
                    conda_activate=conda_activate,
                    conda_environment=conda_environment,
                    model=model,
                    region=region,
                    scenario=scenario,
                    energy=energy,
                    years=args.years,
                    partition=args.partition,
                    nodes=args.nodes,
                    cpus_per_task=args.cpus_per_task,
                    walltime=args.walltime,
                    chunk_time=args.chunk_time,
                    compression=args.compress_level,
                    overwrite=args.overwrite,
                    job_name=job_name,
                )
                write_script(script_path, content)
                created.append(script_path)
                print(f"已创建：{script_path}")

    print(
        "生成完成：{} 个能源类型 × {} 个区域 × {} 个情景 = {} 个 Slurm 作业".format(
            len(energies),
            len(regions),
            len(scenarios),
            len(created),
        )
    )
    print("本生成器不会自动提交作业，请检查脚本后再使用 sbatch。")
    return created


def main(argv: Sequence[str] | None = None) -> int:
    """命令行入口。"""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        create_jobs(args)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
