"""The paired PD metrics of gate version 3 (#104)."""

import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1] / "data/iberia/consistency/benchmark"
sys.path.insert(0, str(BENCH))
from benchlib import DESIGNS, PDTree, draw_designs, paired_pd  # noqa: E402

GENERA = [f"G{g}" for g in range(110)]
UNITS = {f"G{g}_u{i}": f"G{g}" for g in range(110) for i in range(6)}


def _tree(tmp_path, name, stretch=1.0, drop=()):
    tips = [f"{u}:{1.0 + (hash(u) % 7) * stretch}" for u in sorted(UNITS) if UNITS[u] not in drop]
    path = tmp_path / f"{name}.tre"
    path.write_text("(" + ",".join(tips) + ");")
    return PDTree(path)


def test_draws_do_not_depend_on_the_order_of_the_units():
    a = draw_designs(list(UNITS), UNITS, 6, draws=5)
    b = draw_designs(list(reversed(list(UNITS))), UNITS, 6, draws=5)
    assert a == b and set(a) == set(DESIGNS) and all(len(v) == 5 for v in a.values())


def test_identical_settings_give_identical_values(tmp_path):
    clean = _tree(tmp_path, "clean")
    run = _tree(tmp_path, "run", stretch=1.3)
    values = paired_pd({"legacy": clean, "default": clean}, [{"legacy": run, "default": run}], UNITS, 6, 20)
    # the same samples for both settings: the same trees give exactly the same correlations
    assert values["legacy"] == values["default"]
    same = paired_pd({"legacy": clean, "default": clean}, [{"legacy": clean, "default": clean}], UNITS, 6, 20)
    assert all(abs(v - 1.0) < 1e-12 for v in same["default"].values())


def test_samples_come_from_the_units_in_all_four_trees(tmp_path):
    # the default's run lost five genera; its samples, and the legacy's, are drawn without them
    clean = _tree(tmp_path, "clean")
    dropped = _tree(tmp_path, "dropped", drop={f"G{g}" for g in range(5)})
    values = paired_pd({"legacy": clean, "default": clean}, [{"legacy": clean, "default": dropped}], UNITS, 6, 10)
    assert all(abs(v - 1.0) < 1e-12 for v in values["legacy"].values())
    assert all(abs(v - 1.0) < 1e-12 for v in values["default"].values())  # no missing unit was ever sampled
