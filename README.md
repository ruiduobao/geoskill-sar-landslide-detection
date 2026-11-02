# SAR Landslide Detection (geoskill-sar-landslide-detection)

> Fuses InSAR deformation rate, backscatter change and DEM slope into a weighted composite score, extracts suspected landslide connected components with dual thresholds, vectorizes and grades them (high/medium/low), and outputs landslide GeoJSON, deformation-rate/risk-score GeoTIFFs and a risk summary JSON. SAR landslide detection fusing InSAR deformation, backscatter change and slope.

---

## 1. Overview

Identifies suspected landslides by fusing multi-source SAR-derived factors: - **InSAR deformation rate** (mm/yr): landslide bodies show high deformation along the line of sight (absolute value taken). - **Backscatter change**: sliding/turbulence alters surface roughness, so σ⁰ differs strongly between before and after. - **Slope** (derived from a DEM via the Horn gradient method): landslides mostly occur on steep slopes, making slope a key constraint. Method: after robust percentile normalization (`--normalize robust|minmax`) of the three factors, a weighted composite risk score is computed; suspected areas are then extracted with the dual thresholds `score ≥ --score-threshold` and `slope ≥ --slope-threshold`, cleaned morphologically, vectorized into polygons by connected components and graded.

## 3. Quick Start

```bash
pip install -r requirements.txt
python geoskill-sar-landslide-detection.py --bbox 116 39 117 40 --synthetic --output-dir ./out
```

## 4. CLI Parameters

Run `python <skill>.py --help` for the full list. Common parameters:

| Parameter | Type | Description |
|---|---|---|
| `--bbox` | `float[4]` | WGS84 bounding box `min_lon min_lat max_lon max_lat` |
| `--input` | `path` | Local input file (GeoJSON/GeoTIFF/etc.) |
| `--output-dir` | `path` | Output directory (default `./output`) |
| `--synthetic` | `flag` | Use synthetic data instead of real input |
| `--quiet` | `flag` | Suppress non-essential stdout |

## 5. Input / Output

| File | Format | Description |
|---|---|---|
| `landslides.geojson` | GeoJSON | Suspected landslide polygons (with area, score, risk level), EPSG:4326 |
| `deformation_rate.tif` | GeoTIFF (float32) | Deformation rate (mm/yr) |
| `risk_score.tif` | GeoTIFF (float32) | Composite risk score [0,1] |
| `risk_summary.json` | JSON | Count / area / grade summary |
| `output-manifest.json` | JSON | Run manifest |

## 6. Technical Principle

(see SKILL.md for details)

## 8. License

MIT License. See [`LICENSE`](./LICENSE) for full text.

---

> 中文完整版：[references/SKILL.zh-CN.md](references/SKILL.zh-CN.md)
