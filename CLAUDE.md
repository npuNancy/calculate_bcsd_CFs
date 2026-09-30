# CLAUDE.md

全格点入口为 `patchify_grid_cf.py`，输入校验及网格读写位于 `grid_cf_io.py`；
原场站入口 `patchify_station_cf.py` 保留。共用物理核为 `cf_physics.py`。
使用仓库 `.venv/bin/python`，不改变风光物理公式。

一个 grid unit 为 model × scenario × patch × tech。只读取当前 BCSD 根下的年份 blocks，
通过 manifest、sidecar 和 land plan 验证身份；缺失单位须用显式参数声明。
默认8个年份 worker，每个 worker 内按时间和空间分块写 `(time, lat, lon)` NetCDF；
`--merge-final` 默认 False。缺测与零出力分开，年份文件持久保存。

输入、源码和配置身份决定恢复，进程数与合并开关不改变年份身份。
发布采用临时文件原子替换，unit/文件锁防止并发覆盖。
全格点超算运行包位于 infos/scnet_patchify_grid/；正式运行遵循其中 goal.md 和运行前准备.md。
四模型共1,128个unit，14个worker直接写1872的home共享目录；原场站超算生成器保留。

太阳能采用 Erbs 双轴跟踪和温度修正；风电采用 GE120/2500、100 m、1/7 幂律和 25 m/s 切出。
每次回复用户称呼“小凯”，结尾使用“希望对你有帮助，小凯！”。
