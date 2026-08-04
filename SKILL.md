---
name: geoskill-sar-landslide-detection
description: '融合InSAR形变速率、后向散射变化与DEM坡度综合加权评分，双门限提取疑似滑坡连通域并矢量化分级（high/medium/low），输出滑坡GeoJSON、形变速率/风险评分GeoTIFF与风险汇总JSON。SAR landslide detection fusing InSAR deformation, backscatter change and slope.'
---

# SAR 滑坡检测 | SAR Landslide Detection

Identifies suspected landslide bodies by fusing multi-source SAR-derived factors:

- **InSAR deformation rate** (mm/yr): landslide bodies exhibit high deformation along the line of sight (absolute value taken).
- **Backscatter change**: sliding / churning alters surface roughness, producing large σ⁰ differences between before and after.
- **Slope** (derived from a DEM via Horn's gradient method): landslides mostly occur on steep slopes, making slope a key constraint.

Method: the three factors are normalized using robust percentile scaling (`--normalize robust|minmax`), then weighted to compute an integrated risk score; suspected areas are extracted with the dual thresholds `score ≥ --score-threshold` and `slope ≥ --slope-threshold`, cleaned morphologically, vectorized into polygons by connected components, and classified into risk levels.

## Dependencies / 依赖

```bash
pip install numpy rasterio geopandas shapely scipy
```

## Usage / 使用方法

### Basic usage (bbox only, synthetic data auto-generated)

```bash
python geoskill-sar-landslide-detection.py --bbox 116.0 39.0 117.0 40.0 --slope-threshold 15 --output-dir ./out
```

### Example 1: synthetic data (offline)

```bash
python geoskill-sar-landslide-detection.py --bbox 116 39 117 40 --synthetic --output-dir ./syn
```

### Example 2: real deformation rate + DEM

```bash
python geoskill-sar-landslide-detection.py --input deform_rate.tif --dem dem.tif --output-dir ./real
```

### Example 3: with before/after σ⁰ imagery

```bash
python geoskill-sar-landslide-detection.py --input deform.tif --dem dem.tif --sigma-before s0_before.tif --sigma-after s0_after.tif --output-dir ./full
```

### Example 4: lower thresholds for higher sensitivity

```bash
python geoskill-sar-landslide-detection.py --bbox 103 30 104 31 --slope-threshold 10 --score-threshold 0.4 --output-dir ./sensitive --quiet
```

## Output / 输出

| File | Format | Description |
|---|---|---|
| `landslides.geojson` | GeoJSON | Suspected landslide polygons (with area, score, risk level), EPSG:4326 |
| `deformation_rate.tif` | GeoTIFF (float32) | Deformation rate (mm/yr) |
| `risk_score.tif` | GeoTIFF (float32) | Integrated risk score [0,1] |
| `risk_summary.json` | JSON | Count / area / class summary |
| `output-manifest.json` | JSON | Run manifest |

## Data Source / 数据源 / Source

- **Real mode**: local InSAR deformation-rate GeoTIFF, with optional DEM and before/after σ⁰ imagery.
- **Synthetic mode**: locally generated DEM slopes + localized high-deformation patches + σ⁰ anomalies.

## Privacy / 隐私声明 / Privacy

- Fully offline by default; `--synthetic` mode makes no network calls.
- All processing is done locally; no user data is uploaded.

## License / License

MIT

---

<!-- ===== 中文原文 (Chinese Original) ===== -->

---
name: geoskill-sar-landslide-detection
description: '融合InSAR形变速率、后向散射变化与DEM坡度综合加权评分，双门限提取疑似滑坡连通域并矢量化分级（high/medium/low），输出滑坡GeoJSON、形变速率/风险评分GeoTIFF与风险汇总JSON。SAR landslide detection fusing InSAR deformation, backscatter change and slope.'
---

# SAR 滑坡检测 | SAR Landslide Detection

融合多源 SAR 派生因子识别疑似滑坡体：

- **InSAR 形变速率**（mm/yr）：滑坡体在视线向上表现为高形变（取绝对值）。
- **后向散射变化**：滑动 / 翻搅使地表粗糙度改变，σ⁰ 前后差异大。
- **坡度**（由 DEM 经 Horn 梯度法求）：滑坡多发生在陡坡，是关键约束。

方法：三因子稳健百分位归一化（`--normalize robust|minmax`）后加权求综合风险
评分，再用 `score ≥ --score-threshold` 且 `slope ≥ --slope-threshold` 双门限
提取疑似区，形态学清理后按连通域矢量化为多边形并分级。

## 依赖

```bash
pip install numpy rasterio geopandas shapely scipy
```

## 使用方法

### 基本用法（仅给 bbox，自动合成）

```bash
python geoskill-sar-landslide-detection.py --bbox 116.0 39.0 117.0 40.0 --slope-threshold 15 --output-dir ./out
```

### 示例 1：合成数据（离线）

```bash
python geoskill-sar-landslide-detection.py --bbox 116 39 117 40 --synthetic --output-dir ./syn
```

### 示例 2：真实形变速率 + DEM

```bash
python geoskill-sar-landslide-detection.py --input deform_rate.tif --dem dem.tif --output-dir ./real
```

### 示例 3：含 σ⁰ 前后影像

```bash
python geoskill-sar-landslide-detection.py --input deform.tif --dem dem.tif --sigma-before s0_before.tif --sigma-after s0_after.tif --output-dir ./full
```

### 示例 4：降低门限提高灵敏度

```bash
python geoskill-sar-landslide-detection.py --bbox 103 30 104 31 --slope-threshold 10 --score-threshold 0.4 --output-dir ./sensitive --quiet
```

## 输出

| 文件 | 格式 | 说明 |
|---|---|---|
| `landslides.geojson` | GeoJSON | 疑似滑坡多边形（含面积、评分、风险等级），EPSG:4326 |
| `deformation_rate.tif` | GeoTIFF (float32) | 形变速率（mm/yr） |
| `risk_score.tif` | GeoTIFF (float32) | 综合风险评分 [0,1] |
| `risk_summary.json` | JSON | 数量 / 面积 / 分级汇总 |
| `output-manifest.json` | JSON | 运行清单 |

## 数据源 / Source

- **真实模式**：本地 InSAR 形变速率 GeoTIFF，可选 DEM 与 σ⁰ 前后影像。
- **合成模式**：本地生成 DEM 斜坡 + 局部高形变斑块 + σ⁰ 异常。

## 隐私声明 / Privacy

- 默认完全离线运行，`--synthetic` 无任何网络。
- 所有处理本地完成，不上传用户数据。

## License

MIT
