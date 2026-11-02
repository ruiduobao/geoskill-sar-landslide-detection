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

> 中文完整版：需要时再读 [references/SKILL.zh-CN.md](references/SKILL.zh-CN.md)，不要默认全量加载。
