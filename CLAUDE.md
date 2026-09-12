# CLAUDE.md

本仓库只维护 global_bcsd patchify 场站 CF 入口。每个作业是
`model × scenario × patch × tech`，入口为 `patchify_station_cf.py`，物理核为
`cf_physics.py`；禁止恢复已删除的旧批处理入口或网格输出。

使用仓库 `.venv/bin/python`。BCSD 最终文件必须有 sidecar，且 patch_id/variable 与命令参数一致。
输出为 station-only NetCDF，临时文件原子替换；作业生成器位于 `scnet/create_cf_patch_jobs.py`，只生成脚本、不提交。

太阳能采用 Erbs 双轴跟踪和温度修正；风电采用 GE120/2500、100 m、1/7 幂律和 25 m/s 切出。
所有输出和 SCNet 运行态目录由 `.gitignore` 忽略。

每次回复用户称呼“小凯”，结尾使用“希望对你有帮助，小凯！”。

