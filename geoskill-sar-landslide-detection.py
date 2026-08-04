#!/usr/bin/env python3
"""sar-landslide-detection — SAR 滑坡检测

融合多源 SAR 派生因子识别疑似滑坡体：

- **InSAR 形变速率**（mm/yr）：滑坡体在视线向上表现为高形变（取绝对值）。
- **后向散射变化**：滑动 / 翻搅使地表粗糙度改变，σ⁰ 前后差异大。
- **坡度**（由 DEM 求 Horn 梯度）：滑坡多发生在陡坡，坡度是重要约束。

方法：三因子分别做稳健百分位归一化后加权求综合风险评分，再用
``score ≥ score_threshold`` 且 ``slope ≥ slope_threshold`` 双门限提取疑似
区，形态学清理后按连通域矢量化为多边形。

数据源：本地形变速率 GeoTIFF（可选 DEM、σ⁰ 前后影像），或 ``--synthetic``
生成 DEM 斜坡 + 局部高形变斑块 + σ⁰ 异常的模拟场景。

隐私声明 / Privacy：
- 默认完全离线，``--synthetic`` 无网络。
- 所有处理本地完成，不上传用户数据。

Usage:
    python sar-landslide-detection.py --input deform.tif --dem dem.tif --output-dir ./out
    python sar-landslide-detection.py --bbox 116 39 117 40 --slope-threshold 15 --output-dir ./out

License: MIT
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

VERSION = "1.0.0"
SKILL_NAME = "sar-landslide-detection"

try:
    from _geoskill_core.errors import (
        GeoSkillError, UsageError, ValidationError, ProcessError, to_exit_code,
    )
    from _geoskill_core.manifest import OutputManifest, OutputFile
except ImportError:  # pragma: no cover - fallback minimal definitions
    class GeoSkillError(Exception):
        def __init__(self, message: str, code: int = 7, kind: str = "EGeo", **kw):
            super().__init__(message)
            self.message, self.code, self.kind = message, code, kind

    class UsageError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=2, kind="EUsage", **k)

    class ValidationError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=6, kind="EValidate", **k)

    class ProcessError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=7, kind="EProcess", **k)

    def to_exit_code(exc):
        return getattr(exc, "code", 7)

    OutputManifest = None
    OutputFile = None


def _utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat().replace("+00:00", "Z")


# ---------------------------------------------------------------------------
# 核心算法
# ---------------------------------------------------------------------------
def compute_slope(dem: np.ndarray, resolution_m: float = 30.0) -> np.ndarray:
    """由 DEM 计算坡度（度），Horn 梯度法。

    Horn 法用 3×3 加权差分估计 dz/dx、dz/dy，scipy 的 Sobel 卷积核恰好等于
    Horn 分子（未除 8·cellsize），故先除 8 再除像元尺寸得到无量纲梯度，
    ``arctan`` 后转角度。
    """
    from scipy.ndimage import sobel
    z = np.nan_to_num(np.asarray(dem, dtype=np.float64), nan=0.0)
    res = max(float(resolution_m), 1e-6)
    gx = sobel(z, axis=1) / 8.0
    gy = sobel(z, axis=0) / 8.0
    grad = np.sqrt(gx ** 2 + gy ** 2) / res
    return np.degrees(np.arctan(grad)).astype(np.float32)


def backscatter_change(before: np.ndarray, after: np.ndarray) -> np.ndarray:
    """相对后向散射变化 |σ⁰_after − σ⁰_before| / σ⁰_before（线性功率）。"""
    b = np.clip(np.nan_to_num(before, nan=0.0).astype(np.float64), 1e-4, None)
    a = np.clip(np.nan_to_num(after, nan=0.0).astype(np.float64), 1e-4, None)
    return np.abs(a - b) / b


def _normalize(
    x: np.ndarray, lo_pct: float = 2.0, hi_pct: float = 98.0, method: str = "robust"
) -> np.ndarray:
    """归一化到 [0,1]。

    ``robust``：用 [lo_pct, hi_pct] 分位数裁剪异常值后缩放（抗离群点）；
    ``minmax``：用全局最小 / 最大值缩放。
    """
    x = np.nan_to_num(np.asarray(x, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    if x.size == 0:
        return x
    if method == "minmax":
        lo, hi = float(x.min()), float(x.max())
    else:
        lo, hi = np.percentile(x, [lo_pct, hi_pct])
    if hi <= lo:
        return np.zeros_like(x)
    return np.clip((x - lo) / (hi - lo), 0.0, 1.0)


def detect_landslides(
    deform_rate: np.ndarray,
    slope: np.ndarray,
    bs_change: np.ndarray,
    slope_threshold: float = 15.0,
    score_threshold: float = 0.5,
    weights: Tuple[float, float, float] = (0.5, 0.25, 0.25),
    normalize_method: str = "robust",
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """综合风险评分 + 双门限提取疑似滑坡掩膜。

    返回 ``(mask_uint8, score_float32, params)``。评分 = 加权的（|形变|、
    后向散射变化、坡度）归一化值。掩膜 = ``score ≥ score_threshold`` 且
    ``slope ≥ slope_threshold``，经形态学闭 / 开运算清理。
    """
    wd, wb, ws = weights
    d = _normalize(np.abs(deform_rate), method=normalize_method)
    b = _normalize(bs_change, method=normalize_method)
    s = _normalize(slope, method=normalize_method)
    score = (wd * d + wb * b + ws * s).astype(np.float32)

    mask = (score >= score_threshold) & (slope >= slope_threshold)
    from scipy.ndimage import binary_opening, binary_closing
    mask = binary_closing(mask, structure=np.ones((3, 3), dtype=bool))
    mask = binary_opening(mask, structure=np.ones((3, 3), dtype=bool))

    params = {
        "slope_threshold_deg": float(slope_threshold),
        "score_threshold": float(score_threshold),
        "normalize_method": normalize_method,
        "weights": {"deformation": float(wd), "backscatter": float(wb), "slope": float(ws)},
        "mean_score": float(score.mean()),
        "max_score": float(score.max()),
    }
    return mask.astype(np.uint8), score, params


def mask_to_polygons(
    mask: np.ndarray, transform: Any
) -> List[Tuple[Any, np.ndarray]]:
    """把二值掩膜按连通域矢量化为 (shapely 几何, 区域布尔数组) 列表。"""
    from rasterio.features import shapes
    from scipy.ndimage import label
    from shapely.geometry import shape

    labels, n = label(mask.astype(bool))
    feats: List[Tuple[Any, np.ndarray]] = []
    for i in range(1, n + 1):
        region = labels == i
        polys = []
        for geom, val in shapes(region.astype("uint8"), mask=region, transform=transform):
            if val == 1:
                polys.append(shape(geom))
        if not polys:
            continue
        geom = max(polys, key=lambda p: p.area)
        feats.append((geom, region))
    return feats


def build_risk_features(
    feats: List[Tuple[Any, np.ndarray]],
    score: np.ndarray,
    deform_rate: np.ndarray,
    pixel_km2: float,
) -> List[Dict[str, Any]]:
    """为每个连通域计算属性并划分风险等级。"""
    records: List[Dict[str, Any]] = []
    for idx, (geom, region) in enumerate(feats, start=1):
        n_px = int(region.sum())
        mean_score = float(score[region].mean())
        max_deform = float(np.abs(deform_rate[region]).max())
        if mean_score >= 0.75:
            level = "high"
        elif mean_score >= 0.6:
            level = "medium"
        else:
            level = "low"
        records.append({
            "id": idx,
            "geometry": geom,
            "area_km2": n_px * pixel_km2,
            "pixels": n_px,
            "mean_score": mean_score,
            "max_deformation_mm_yr": max_deform,
            "risk_level": level,
        })
    records.sort(key=lambda r: r["mean_score"], reverse=True)
    for new_id, rec in enumerate(records, start=1):
        rec["id"] = new_id
    return records


def pixel_area_km2(bbox: List[float], height: int, width: int) -> float:
    lat_mid = 0.5 * (bbox[1] + bbox[3])
    km_per_deg_lon = 111.0 * float(np.cos(np.deg2rad(lat_mid)))
    px_w = (bbox[2] - bbox[0]) / max(width, 1) * km_per_deg_lon
    px_h = (bbox[3] - bbox[1]) / max(height, 1) * 111.0
    return float(abs(px_w * px_h))


def risk_summary(
    records: List[Dict[str, Any]], params: Dict[str, Any], bbox: List[float]
) -> Dict[str, Any]:
    levels = {"high": 0, "medium": 0, "low": 0}
    area_by_level = {"high": 0.0, "medium": 0.0, "low": 0.0}
    for r in records:
        levels[r["risk_level"]] += 1
        area_by_level[r["risk_level"]] += r["area_km2"]
    return {
        "n_landslides": len(records),
        "total_area_km2": sum(r["area_km2"] for r in records),
        "count_by_level": levels,
        "area_by_level_km2": area_by_level,
        "max_deformation_mm_yr": max((r["max_deformation_mm_yr"] for r in records), default=0.0),
        "max_score": max((r["mean_score"] for r in records), default=0.0),
        "bbox_wgs84": bbox,
        **params,
    }


# ---------------------------------------------------------------------------
# 合成数据
# ---------------------------------------------------------------------------
def generate_synthetic(
    bbox: List[float],
    width: int = 64,
    height: int = 64,
    seed: int = 42,
    resolution_m: float = 30.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
    """生成 DEM + 形变速率 + σ⁰ 前后影像，内含两处高形变陡坡滑坡斑块。

    返回 ``(dem, deform_rate, sigma_before, sigma_after, truth, info)``。
    """
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:height, 0:width]

    # 缓坡背景（约 5°）：2.6 m/pixel 高差 / 30 m → arctan ≈ 5°
    dem = (yy * 2.6 + xx * 0.5).astype(np.float32)
    dem = dem + rng.normal(0, 0.8, (height, width)).astype(np.float32)

    deform = rng.normal(0, 3, (height, width)).astype(np.float32)  # 背景 ~0 mm/yr
    sigma_before = np.full((height, width), 0.05, dtype=np.float32)
    sigma_after = sigma_before.copy()

    truth = np.zeros((height, width), dtype=np.uint8)
    patches = [(12, 12, 26, 26), (38, 36, 52, 50)]
    for (r0, c0, r1, c1) in patches:
        r1, c1 = min(r1, height), min(c1, width)
        rr, cc = r1 - r0, c1 - c0
        if rr <= 0 or cc <= 0:
            continue
        # 滑坡后壁陡坎：沿行快速抬升 → 局部高坡度（约 25°）
        scarp = (np.arange(rr) * 14.0).astype(np.float32)
        dem[r0:r1, c0:c1] = dem[r0:r1, c0:c1] + scarp[:, None]
        # 高形变（视线向位移）
        deform[r0:r1, c0:c1] = -80.0 + rng.normal(0, 5, (rr, cc)).astype(np.float32)
        # 地表翻搅 → 后向散射增强
        sigma_after[r0:r1, c0:c1] = 0.05 + rng.uniform(0.03, 0.06, (rr, cc)).astype(np.float32)
        truth[r0:r1, c0:c1] = 1

    # 乘性斑点
    speckle_b = np.exp(rng.normal(0, 0.08, (height, width))).astype(np.float32)
    speckle_a = np.exp(rng.normal(0, 0.08, (height, width))).astype(np.float32)
    sigma_before = sigma_before * speckle_b
    sigma_after = sigma_after * speckle_a

    info = {
        "bbox": bbox, "width": width, "height": height, "seed": seed,
        "resolution_m": resolution_m,
        "patches": [list(p) for p in patches],
        "truth_fraction": float(truth.mean()),
    }
    return dem, deform, sigma_before, sigma_after, truth, info


# ---------------------------------------------------------------------------
# GeoTIFF / 矢量 I/O
# ---------------------------------------------------------------------------
def write_geotiff(
    path: str, cube: np.ndarray, bbox: List[float],
    nodata: float = -9999.0, dtype: str = "float32",
):
    import rasterio
    from rasterio.transform import from_bounds
    if cube.ndim == 2:
        cube = cube[np.newaxis, ...]
    nb, h, w = cube.shape
    transform = from_bounds(bbox[0], bbox[1], bbox[2], bbox[3], w, h)
    profile = {
        "driver": "GTiff", "height": h, "width": w, "count": nb,
        "dtype": dtype, "crs": "EPSG:4326", "transform": transform,
        "nodata": nodata, "compress": "deflate",
    }
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with rasterio.open(path, "w", **profile) as dst:
        for b in range(nb):
            dst.write(cube[b].astype(dtype), b + 1)
    return transform


def read_geotiff(path: str) -> Tuple[np.ndarray, List[float]]:
    import rasterio
    if not os.path.exists(path):
        raise UsageError(f"input raster not found: {path}", path=path)
    with rasterio.open(path) as src:
        cube = src.read().astype(np.float32)
        b = src.bounds
        bbox = [b.left, b.bottom, b.right, b.top]
    return cube, bbox


def read_geotiff_full(path: str) -> Tuple[np.ndarray, List[float], Optional[float]]:
    """扩展版 read：同时返回 nodata 值（若无则为 None）。"""
    import rasterio
    if not os.path.exists(path):
        raise UsageError(f"input raster not found: {path}", path=path)
    with rasterio.open(path) as src:
        cube = src.read().astype(np.float32)
        b = src.bounds
        bbox = [b.left, b.bottom, b.right, b.top]
        nodata = src.nodata
        if nodata is not None:
            nodata = float(nodata)
    return cube, bbox, nodata


# ---------------------------------------------------------------------------
# 输入校验
# ---------------------------------------------------------------------------
def validate_bbox(bbox: Optional[List[float]]) -> List[float]:
    """校验地理 bbox 合法性，失败抛 ValidationError（exit 6）。"""
    if bbox is None or len(bbox) != 4:
        raise ValidationError("bbox must be 4 floats [W S E N]")
    try:
        w, s, e, n = [float(x) for x in bbox]
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"bbox entries must be numeric: {exc}")
    if not (-90.0 <= s <= 90.0 and -90.0 <= n <= 90.0):
        raise ValidationError(f"latitude out of [-90,90]: S={s}, N={n}")
    if not (-180.0 <= w <= 180.0 and -180.0 <= e <= 180.0):
        raise ValidationError(f"longitude out of [-180,180]: W={w}, E={e}")
    if s >= n:
        raise ValidationError(
            f"S >= N (S={s}, N={n}); bbox inverted (S must be < N)"
        )
    if w >= e:
        raise ValidationError(
            f"W >= E (W={w}, E={e}); cross-180° bbox not supported. "
            f"Split into two non-antipodal bboxes."
        )
    if (e - w) < 0.001 or (n - s) < 0.001:
        raise ValidationError(
            f"bbox too small ({(e-w):.6f}°×{(n-s):.6f}°); min span is 0.001°"
        )
    return [w, s, e, n]


def write_geojson(path: str, records: List[Dict[str, Any]]) -> None:
    import geopandas as gpd
    from shapely.geometry import mapping  # noqa: F401
    rows = []
    for r in records:
        rows.append({
            "id": r["id"],
            "area_km2": r["area_km2"],
            "pixels": r["pixels"],
            "mean_score": r["mean_score"],
            "max_deformation_mm_yr": r["max_deformation_mm_yr"],
            "risk_level": r["risk_level"],
            "geometry": r["geometry"],
        })
    if rows:
        gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    else:
        import pandas as pd
        gdf = gpd.GeoDataFrame(
            pd.DataFrame(columns=["id", "area_km2", "pixels", "mean_score",
                                  "max_deformation_mm_yr", "risk_level"]),
            geometry=[], crs="EPSG:4326",
        )
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    gdf.to_file(path, driver="GeoJSON")


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
def write_manifest(
    output_dir, args, outputs, qa, started_at, exit_code, bbox,
    input_nodata: Optional[float] = None,
) -> Optional[str]:
    if OutputManifest is None:
        return None
    cmd = " ".join([SKILL_NAME] + sys.argv[1:])
    man = OutputManifest(
        skill=SKILL_NAME, skill_version=VERSION, command=cmd,
        started_at=started_at, finished_at=_utc_now(), exit_code=exit_code,
        inputs={
            "input": getattr(args, "input", None),
            "dem": getattr(args, "dem", None),
            "slope_threshold": getattr(args, "slope_threshold", None),
            "score_threshold": getattr(args, "score_threshold", None),
            "synthetic": bool(getattr(args, "synthetic", False)),
            "input_nodata": input_nodata,
        },
        outputs=[OutputFile(**o) for o in outputs],
        qa=qa,
        software={"python": sys.version.split()[0], "skill": SKILL_NAME},
    )
    path = os.path.join(output_dir, "output-manifest.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(man.to_dict(), f, ensure_ascii=False, indent=2, default=str)
    return path


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def process(args: argparse.Namespace) -> int:
    started_at = _utc_now()
    output_dir = args.output_dir
    bbox = list(args.bbox) if args.bbox else None
    truth = None
    synth_info: Optional[Dict[str, Any]] = None
    input_nodata: Optional[float] = None
    n_valid_pixels: Optional[int] = None

    # 校验 CLI 参数（前置）
    if not (0.0 <= args.slope_threshold <= 90.0):
        raise ValidationError(
            f"--slope-threshold must be in [0, 90] degrees (got {args.slope_threshold})"
        )
    if not (0.0 <= args.score_threshold <= 1.0):
        raise ValidationError(
            f"--score-threshold must be in [0, 1] (got {args.score_threshold})"
        )
    if args.resolution <= 0:
        raise ValidationError(
            f"--resolution must be > 0 (got {args.resolution})"
        )

    if args.input and not args.synthetic:
        cube, file_bbox, src_nodata = read_geotiff_full(args.input)
        input_nodata = src_nodata
        if bbox is not None:
            bbox = validate_bbox(bbox)
        else:
            bbox = validate_bbox(file_bbox)
        deform = cube[0]
        # NoData 处理
        if src_nodata is not None:
            n_total = int(deform.size)
            n_nd = int(np.count_nonzero(deform == src_nodata))
            n_valid_pixels = n_total - n_nd
            if n_valid_pixels == 0:
                raise ValidationError(
                    f"input raster has no valid pixels "
                    f"(all {n_nd}/{n_total} are NoData={src_nodata})",
                    path=args.input, nodata=src_nodata,
                )
            deform = np.where(deform == src_nodata, np.nan, deform).astype(np.float32)
        else:
            n_valid_pixels = int(deform.size)
        h, w = deform.shape
        if args.dem:
            dem_cube, _, _ = read_geotiff_full(args.dem)
            dem = dem_cube[0]
        else:
            dem = None
        if args.sigma_before and args.sigma_after:
            sb, _, _ = read_geotiff_full(args.sigma_before)
            sa, _, _ = read_geotiff_full(args.sigma_after)
            bs_change = backscatter_change(sb[0], sa[0])
        else:
            bs_change = np.zeros_like(deform)
        source_note = args.input
    else:
        if bbox is None:
            raise UsageError("provide --bbox (synthetic mode) or --input <deformation-rate raster>")
        bbox = validate_bbox(bbox)
        dem, deform, sb, sa, truth, synth_info = generate_synthetic(
            bbox, resolution_m=args.resolution,
        )
        bs_change = backscatter_change(sb, sa)
        n_valid_pixels = int(deform.size)
        source_note = "synthetic"

    if deform.size == 0:
        raise ValidationError("input raster is empty")

    # 校验通过后再建输出目录
    os.makedirs(output_dir, exist_ok=True)

    h, w = deform.shape
    if dem is not None:
        slope = compute_slope(dem, resolution_m=args.resolution)
    else:
        # 无 DEM 时用常数坡度 0，则坡度门限需 ≤0 才会命中；退化为仅靠形变/后向散射
        slope = np.zeros_like(deform, dtype=np.float32)

    mask, score, params = detect_landslides(
        deform, slope, bs_change,
        slope_threshold=args.slope_threshold,
        score_threshold=args.score_threshold,
        normalize_method=args.normalize,
    )

    transform = write_geotiff(
        os.path.join(output_dir, "deformation_rate.tif"), deform, bbox,
    )
    write_geotiff(os.path.join(output_dir, "risk_score.tif"), score, bbox)
    if dem is not None:
        write_geotiff(os.path.join(output_dir, "slope.tif"), slope, bbox)

    px = pixel_area_km2(bbox, h, w)
    feats = mask_to_polygons(mask, transform)
    records = build_risk_features(feats, score, deform, px)
    summary = risk_summary(records, params, bbox)

    geojson_path = os.path.join(output_dir, "landslides.geojson")
    write_geojson(geojson_path, records)
    summary_path = os.path.join(output_dir, "risk_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    qa: Dict[str, Any] = {
        "source": source_note,
        "n_landslides": summary["n_landslides"],
        "total_area_km2": summary["total_area_km2"],
        "count_by_level": summary["count_by_level"],
        "max_deformation_mm_yr": summary["max_deformation_mm_yr"],
        "n_valid_pixels": int(n_valid_pixels) if n_valid_pixels is not None else None,
        "input_nodata": input_nodata,
    }
    if synth_info is not None:
        qa["synthetic_truth_fraction"] = synth_info["truth_fraction"]

    outputs = [
        {"path": geojson_path, "kind": "vector", "crs_epsg": 4326, "bbox_wgs84": bbox},
        {"path": os.path.join(output_dir, "deformation_rate.tif"), "kind": "raster",
         "crs_epsg": 4326, "bbox_wgs84": bbox, "band_count": 1},
        {"path": os.path.join(output_dir, "risk_score.tif"), "kind": "raster",
         "crs_epsg": 4326, "bbox_wgs84": bbox, "band_count": 1},
        {"path": summary_path, "kind": "json"},
    ]
    man_path = write_manifest(output_dir, args, outputs, qa, started_at, 0, bbox,
                              input_nodata=input_nodata)

    if not args.quiet:
        print(f"[{SKILL_NAME}] source: {source_note}")
        print(f"[{SKILL_NAME}] landslides detected: {summary['n_landslides']}  "
              f"area: {summary['total_area_km2']:.3f} km²")
        print(f"[{SKILL_NAME}] levels: {summary['count_by_level']}")
        print(f"[{SKILL_NAME}] output: {geojson_path}")
        if man_path:
            print(f"[{SKILL_NAME}] manifest: {man_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=SKILL_NAME,
        description="SAR landslide detection fusing InSAR deformation, backscatter change and slope.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bbox", nargs=4, type=float, metavar=("W", "S", "E", "N"),
                   help="geographic extent [minLon minLat maxLon maxLat]")
    p.add_argument("--input", help="InSAR deformation-rate GeoTIFF (mm/yr)")
    p.add_argument("--dem", help="DEM GeoTIFF for slope (optional)")
    p.add_argument("--sigma-before", help="pre-event σ⁰ GeoTIFF (optional, with --sigma-after)")
    p.add_argument("--sigma-after", help="post-event σ⁰ GeoTIFF (optional, with --sigma-before)")
    p.add_argument("--slope-threshold", type=float, default=15.0,
                   help="minimum slope in degrees for a candidate (default: 15)")
    p.add_argument("--score-threshold", type=float, default=0.5,
                   help="minimum composite risk score (default: 0.5)")
    p.add_argument("--normalize", default="robust", choices=["robust", "minmax"],
                   help="factor normalization scheme (default: robust)")
    p.add_argument("--resolution", type=float, default=30.0,
                   help="pixel size in meters for slope computation (default: 30)")
    p.add_argument("--synthetic", action="store_true",
                   help="generate a synthetic DEM + deformation scene (offline)")
    p.add_argument("--output-dir", default="./output", help="output directory")
    p.add_argument("--quiet", action="store_true", help="suppress console output")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return process(args)
    except GeoSkillError as exc:
        print(f"[{SKILL_NAME}] ERROR [{exc.kind}] {exc.message}", file=sys.stderr)
        return exc.code
    except KeyboardInterrupt:
        return 130
    except Exception as exc:  # noqa: BLE001
        print(f"[{SKILL_NAME}] ERROR {exc}", file=sys.stderr)
        return to_exit_code(exc)


if __name__ == "__main__":
    sys.exit(main())
