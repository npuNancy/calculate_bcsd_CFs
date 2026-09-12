# global_bcsd patchify 场站容量因子

本目录只保留基于 `global_bcsd` patch 文件的场站 CF 流程。每个作业处理一个
`model × scenario × patch × tech`，输出只包含 `(time, station)`，不生成格点或全球合并文件。

```bash
python patchify_station_cf.py \
  --bcsd-root ~/project_climate/bcsd/global_bcsd \
  --model CANESM5 --scenario ssp126 --patch R01C01 \
  --patch-manifest /path/to/patch_manifest.json \
  --stations-csv /path/to/stations_SSP1-2.6.csv \
  --tech wind --years 2015-2060 \
  --output-root /path/to/cf_patchify
```

BCSD 最终文件及其 `.json` sidecar 必须声明相同的 `patch_id` 和 `variable`。
场站经度统一到 `[-180, 180)`，支持最近邻或规则网格双线性抽取；风电使用
`GE120/2500`、100 m 幂律外推，光伏使用 Erbs 辐射分解和双轴跟踪物理核。

超算脚本由 `scnet/create_cf_patch_jobs.py`（guide 说明在
`infos/scnet_patchify/`）生成，生成器不提交作业。输出目录默认应设置为新的、
被 `.gitignore` 忽略的 patchify 路径。
