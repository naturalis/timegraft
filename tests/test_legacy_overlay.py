"""The legacy arm of the benchmark is the run before the consistency layer: new defaults must not reach it."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "data/iberia/consistency/benchmark"


def _settings(*overlays):
    spec = importlib.util.spec_from_file_location("backbone_harness", BENCH / "backbone.py")
    harness = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(harness)
    return harness.settings_for(ROOT, [str(BENCH / o) for o in overlays])["backbone_build"]


def test_legacy_calibrates_with_timetree_raw_time():
    assert _settings()["ages"]["timetree_time"]["use"] == "replace"
    assert _settings("legacy.yaml")["ages"]["timetree_time"]["use"] == "raw"


def test_legacy_uses_every_timetree_answer():
    assert _settings()["ages"]["node_check"]["mode"] == "reject"
    assert _settings("legacy.yaml")["ages"]["node_check"]["mode"] == "off"
