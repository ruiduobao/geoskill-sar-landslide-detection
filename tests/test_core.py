"""Core algorithm tests for sar-landslide-detection."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from conftest import mod as ld


def _transform(bbox, h, w):
    from rasterio.transform import from_bounds
    return from_bounds(bbox[0], bbox[1], bbox[2], bbox[3], w, h)


class TestSlope:
    def test_flat_zero(self):
        dem = np.full((32, 32), 100.0, dtype=np.float32)
        slope = ld.compute_slope(dem, resolution_m=30.0)
        assert slope.max() < 0.5

    def test_tilted_plane(self):
        # 每像元抬升 2.62 m，像元 30 m → arctan(2.62/30) ≈ 5°
        yy, _ = np.mgrid[0:32, 0:32]
        dem = (yy * 2.62).astype(np.float32)
        slope = ld.compute_slope(dem, resolution_m=30.0)
        interior = slope[5:-5, 5:-5]
        assert 3.5 < interior.mean() < 6.5

    def test_steeper_gives_larger(self):
        yy, _ = np.mgrid[0:32, 0:32]
        gentle = ld.compute_slope((yy * 2.0).astype(np.float32), 30.0)
        steep = ld.compute_slope((yy * 14.0).astype(np.float32), 30.0)
        assert steep[10:-10, 10:-10].mean() > gentle[10:-10, 10:-10].mean()


class TestBackscatterChange:
    def test_no_change_zero(self):
        s = np.full((16, 16), 0.05, dtype=np.float32)
        ch = ld.backscatter_change(s, s)
        np.testing.assert_allclose(ch, 0.0, atol=1e-6)

    def test_change_positive(self):
        before = np.full((16, 16), 0.05, dtype=np.float32)
        after = np.full((16, 16), 0.10, dtype=np.float32)
        ch = ld.backscatter_change(before, after)
        np.testing.assert_allclose(ch, 1.0, atol=1e-3)


class TestDetect:
    def test_synthetic_detection_overlaps_truth(self):
        dem, deform, sb, sa, truth, info = ld.generate_synthetic([116, 39, 117, 40], seed=5)
        bs_change = ld.backscatter_change(sb, sa)
        slope = ld.compute_slope(dem, resolution_m=30.0)
        mask, score, params = ld.detect_landslides(
            deform, slope, bs_change, slope_threshold=15.0, score_threshold=0.5,
        )
        recall = mask[truth == 1].mean()
        assert recall > 0.5, f"recall too low: {recall:.3f}"
        fp = mask[truth == 0].mean()
        assert fp < 0.15, f"false positive too high: {fp:.3f}"

    def test_no_signal_no_detection(self):
        rng = np.random.default_rng(0)
        dem = np.full((32, 32), 100.0, dtype=np.float32)  # 平地
        deform = rng.normal(0, 0.1, (32, 32)).astype(np.float32)
        bs = np.zeros((32, 32), dtype=np.float32)
        slope = ld.compute_slope(dem, 30.0)
        mask, score, _ = ld.detect_landslides(deform, slope, bs, slope_threshold=15.0)
        assert mask.sum() == 0

    def test_score_range(self):
        dem, deform, sb, sa, truth, info = ld.generate_synthetic([116, 39, 117, 40])
        slope = ld.compute_slope(dem, 30.0)
        bs = ld.backscatter_change(sb, sa)
        _, score, _ = ld.detect_landslides(deform, slope, bs)
        assert score.min() >= 0.0
        assert score.max() <= 1.0


class TestVectorize:
    def test_polygons_from_mask(self):
        mask = np.zeros((32, 32), dtype=np.uint8)
        mask[5:15, 5:15] = 1
        mask[20:28, 20:28] = 1
        transform = _transform([116, 39, 117, 40], 32, 32)
        feats = ld.mask_to_polygons(mask, transform)
        assert len(feats) == 2
        score = np.random.default_rng(0).uniform(0.5, 1.0, (32, 32)).astype(np.float32)
        deform = np.random.default_rng(1).uniform(-100, 0, (32, 32)).astype(np.float32)
        records = ld.build_risk_features(feats, score, deform, pixel_km2=1.0)
        assert len(records) == 2
        assert all(r["area_km2"] > 0 for r in records)
        assert all(r["risk_level"] in ("high", "medium", "low") for r in records)

    def test_empty_mask_no_features(self):
        mask = np.zeros((16, 16), dtype=np.uint8)
        transform = _transform([116, 39, 117, 40], 16, 16)
        feats = ld.mask_to_polygons(mask, transform)
        assert len(feats) == 0

    def test_write_geojson_empty(self, tmp_path):
        path = str(tmp_path / "empty.geojson")
        ld.write_geojson(path, [])
        assert os.path.exists(path)
        import json
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        assert data["type"] == "FeatureCollection"
        assert data["features"] == []


class TestRiskSummary:
    def test_counts_consistent(self):
        dem, deform, sb, sa, truth, info = ld.generate_synthetic([116, 39, 117, 40], seed=11)
        slope = ld.compute_slope(dem, 30.0)
        bs = ld.backscatter_change(sb, sa)
        mask, score, params = ld.detect_landslides(deform, slope, bs, score_threshold=0.4)
        transform = _transform([116, 39, 117, 40], *deform.shape)
        feats = ld.mask_to_polygons(mask, transform)
        records = ld.build_risk_features(feats, score, deform, 1.0)
        summary = ld.risk_summary(records, params, [116, 39, 117, 40])
        assert summary["n_landslides"] == len(records)
        assert summary["n_landslides"] >= 1
        assert summary["total_area_km2"] > 0
        assert sum(summary["count_by_level"].values()) == summary["n_landslides"]

    def test_empty_records(self):
        summary = ld.risk_summary([], {"slope_threshold_deg": 15.0}, [116, 39, 117, 40])
        assert summary["n_landslides"] == 0
        assert summary["max_score"] == 0.0


class TestSynthetic:
    def test_shapes(self):
        dem, deform, sb, sa, truth, info = ld.generate_synthetic([116, 39, 117, 40])
        assert dem.shape == (64, 64)
        assert deform.shape == (64, 64)
        assert truth.sum() > 0
        assert info["truth_fraction"] > 0

    def test_patches_high_deformation(self):
        dem, deform, sb, sa, truth, info = ld.generate_synthetic([116, 39, 117, 40])
        assert np.abs(deform[truth == 1]).mean() > 40
        assert np.abs(deform[truth == 0]).mean() < 10


class TestIO:
    def test_geotiff_roundtrip(self, tmp_path):
        arr = np.random.uniform(-100, 0, (16, 16)).astype(np.float32)
        path = str(tmp_path / "d.tif")
        ld.write_geotiff(path, arr, [116.0, 39.0, 117.0, 40.0])
        back, bbox = ld.read_geotiff(path)
        np.testing.assert_allclose(back[0], arr, atol=1e-4)

    def test_missing_raises(self):
        with pytest.raises(ld.UsageError):
            ld.read_geotiff("/nonexistent/x.tif")
