# grid_of_regions —— 各区域网格定义

本目录存放各国 / 各区域的**目标网格**（仅含坐标，无数据变量），用于在统计时把
全球数据（如 ERA5-Land 0.1° 场）截取到与各区域 BCSD 降尺度一致的空间范围。

> 主要消费者：`S06E02_Annual_Mean_CF_ERA5Land.py`（ERA5-Land 各国 & 全球逐年平均容量因子）。
> 读取这些网格是为了**保证 ERA5-Land 的各国区域与 S06E01（BCSD 降尺度）一致、可直接对比**。

## 文件命名

```
{region}_grid.nc
```

每个文件只含 `lat` / `lon` 坐标（`Data variables: *empty*`），属性里记录：
- `region`：区域名；
- `grid_source`：该网格的来源文件（BCSD 降尺度输出 / CORDEX 原始网格等）。

## 网格类型

共 **28** 个区域，分两类（脚本按 `lat`/`lon` 维度自动识别）：

### 1. 规则 1D 网格（27 个）

`lat (lat)` 与 `lon (lon)` 各为 1D 坐标，构成规则经纬网。
- 经度约定为 **`[0, 360)`**；部分跨 0° 经线的国家（如 France）的 `lon` 同时含 ~0–8° 与 ~355–360°。
- 分辨率随来源不同（多为 0.1°；china 为其 NESM3 native 网格，约 0.1–0.2°）。
- 来源多为 `data/bcsd_outputs/{model}/{region}/...`（china 为 `output/CFs_of_*_china/...`）。

### 2. 曲线 / 旋转极点 2D 网格（1 个：`NAM-12`）

`lat (rlat, rlon)` 与 `lon (rlat, rlon)` 为 2D 坐标（CORDEX NAM-12 旋转极点网格，0.11°）。
- 经度约定为 **`[-180, 180]`**；
- 统计时把 2D 网格**展平为散点**，逐点求在目标全球网格上的最近邻后截取。

## 区域清单

| 文件 | 类型 | 维度 (lat×lon) | 格点数 | 经度约定 |
|---|---|---|---:|---|
| Australia_grid.nc | 规则 1D | 345×407 | 140,415 | [0,360) |
| Austria_grid.nc | 规则 1D | 26×77 | 2,002 | [0,360) |
| Brazil_grid.nc | 规则 1D | 391×392 | 153,272 | [0,360) |
| Chile_grid.nc | 规则 1D | 384×431 | 165,504 | [0,360) |
| Denmark_grid.nc | 规则 1D | 32×71 | 2,272 | [0,360) |
| Egypt_grid.nc | 规则 1D | 97×122 | 11,834 | [0,360) |
| France_grid.nc | 规则 1D | 98×148 | 14,504 | [0,360) |
| Germany_grid.nc | 规则 1D | 81×101 | 8,181 | [0,360) |
| Greece_grid.nc | 规则 1D | 70×88 | 6,160 | [0,360) |
| India_grid.nc | 规则 1D | 311×301 | 93,611 | [0,360) |
| Ireland_grid.nc | 规则 1D | 40×45 | 1,800 | [0,360) |
| Italy_grid.nc | 规则 1D | 121×131 | 15,851 | [0,360) |
| Japan_grid.nc | 规则 1D | 213×229 | 48,777 | [0,360) |
| México_grid.nc | 规则 1D | 182×317 | 57,694 | [0,360) |
| NAM-12_grid.nc | **2D 曲线** | 628×655 | 411,340 | [-180,180] |
| Netherlands_grid.nc | 规则 1D | 28×39 | 1,092 | [0,360) |
| Poland_grid.nc | 规则 1D | 58×100 | 5,800 | [0,360) |
| Portugal_grid.nc | 规则 1D | 52×251 | 13,052 | [0,360) |
| Romania_grid.nc | 规则 1D | 46×94 | 4,324 | [0,360) |
| South-Africa_grid.nc | 规则 1D | 62×165 | 10,230 | [0,360) |
| South-Korea_grid.nc | 规则 1D | 54×63 | 3,402 | [0,360) |
| Spain_grid.nc | 规则 1D | 91×151 | 13,741 | [0,360) |
| Sweden_grid.nc | 规则 1D | 138×130 | 17,940 | [0,360) |
| Turkey_grid.nc | 规则 1D | 63×191 | 12,033 | [0,360) |
| Ukraine_grid.nc | 规则 1D | 80×180 | 14,400 | [0,360) |
| United-Kingdom_grid.nc | 规则 1D | 109×104 | 11,336 | [0,360) |
| Vietnam_grid.nc | 规则 1D | 148×74 | 10,952 | [0,360) |
| china_grid.nc | 规则 1D | 348×442 | 153,816 | [0,360) |

## ⚠️ 注意：NAM-12 与各国的空间重叠

`NAM-12` 是**北美域**，与 `México`（及美国陆地）在空间上**重叠**，且其格点数（~41 万，
陆地约 17 万）远大于其它任何区域。若 `NAM-12` 与 `México` **同时计入全球加权平均**，
会造成北美重复计数，并使 NAM-12 主导权重。

在 `S06E02_Annual_Mean_CF_ERA5Land.py` 中，默认全部计入、由使用者决定口径；
如需避免重复计数，可用 `--exclude-regions NAM-12`（或 `México`）排除其一。

## 读取示例

```python
import xarray as xr

# 规则 1D 网格
g = xr.open_dataset("data/grid_of_regions/France_grid.nc")
lat = g["lat"].values            # (98,)
lon = g["lon"].values            # (148,)

# 2D 曲线网格（NAM-12）
g = xr.open_dataset("data/grid_of_regions/NAM-12_grid.nc")
lat2d = g["lat"].values          # (628, 655)
lon2d = g["lon"].values          # (628, 655)
```
