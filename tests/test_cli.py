"""CLI argument parsing tests for sar-landslide-detection."""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
from conftest import SCRIPT


def run_cli(args, timeout=120):
    return subprocess.run(
        [sys.executable, SCRIPT] + args,
        capture_output=True, text=True, timeout=timeout,
    )


class TestCLIBasics:
    def test_help(self):
        r = run_cli(["--help"])
        assert r.returncode == 0
        assert "--bbox" in r.stdout
        assert "--slope-threshold" in r.stdout
        assert "--synthetic" in r.stdout

    def test_version(self):
        r = run_cli(["--version"])
        assert r.returncode == 0
        assert "1.0.0" in r.stdout

    def test_no_args_exit_2(self):
        r = run_cli([])
        assert r.returncode == 2

    def test_bad_normalize_exit_2(self):
        r = run_cli(["--bbox", "116", "39", "117", "40", "--normalize", "zscore"])
        assert r.returncode == 2

    def test_input_not_found_exit_2(self):
        r = run_cli(["--input", "nonexistent.tif"])
        assert r.returncode == 2

    def test_synthetic_without_bbox_exit_2(self):
        r = run_cli(["--synthetic"])
        assert r.returncode == 2


class TestCLIEndToEnd:
    def test_bbox_only_runs(self, tmp_path):
        out = str(tmp_path / "o1")
        r = run_cli(["--bbox", "116", "39", "117", "40", "--output-dir", out, "--quiet"])
        assert r.returncode == 0, r.stderr
        assert os.path.exists(os.path.join(out, "landslides.geojson"))
        assert os.path.exists(os.path.join(out, "deformation_rate.tif"))
        assert os.path.exists(os.path.join(out, "risk_summary.json"))
        assert os.path.exists(os.path.join(out, "output-manifest.json"))

    def test_custom_thresholds(self, tmp_path):
        out = str(tmp_path / "o2")
        r = run_cli([
            "--bbox", "116", "39", "117", "40", "--synthetic",
            "--slope-threshold", "10", "--score-threshold", "0.4",
            "--output-dir", out, "--quiet",
        ])
        assert r.returncode == 0, r.stderr

    def test_geojson_detects_patches(self, tmp_path):
        """合成场景应检测到注入的滑坡斑块。"""
        import json
        out = str(tmp_path / "o3")
        r = run_cli(["--bbox", "116", "39", "117", "40", "--synthetic",
                     "--output-dir", out, "--quiet"])
        assert r.returncode == 0, r.stderr
        with open(os.path.join(out, "risk_summary.json"), encoding="utf-8") as f:
            summary = json.load(f)
        assert summary["n_landslides"] >= 1
        with open(os.path.join(out, "landslides.geojson"), encoding="utf-8") as f:
            gj = json.load(f)
        assert gj["type"] == "FeatureCollection"
        assert len(gj["features"]) >= 1
