"""The downstream harness of the gate (#81): subsetting units by genus, and the cache of build_genus() results."""

import csv
import importlib.util
import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parents[1] / "data/iberia/consistency/benchmark"
REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def ds():
    sys.path.insert(0, str(BENCH))
    try:
        spec = importlib.util.spec_from_file_location("downstream", BENCH / "downstream.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        yield module
    finally:
        sys.path.remove(str(BENCH))


def _table(path, columns, rows):
    with open(path, "w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        w.writerows(rows)


def _read(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def test_subset_keeps_the_genera_and_their_unassigned_units(ds, tmp_path):
    units = tmp_path / "units"
    units.mkdir()
    for name in ("screened_representatives.tsv", "representatives.tsv"):
        _table(units / name, ["unit_id", "genus"], [["u1", "A"], ["u2", "B"], ["u3", "C"]])
    # records: u1 holds one record of genus C, u3 one of genus A; x1 is an unassigned unit
    _table(units / "filtered.tsv", ["record_id", "unit_id", "genus"],
           [["r1", "u1", "A"], ["r2", "u1", "C"], ["r3", "u2", "B"], ["r4", "u3", "C"], ["r5", "u3", "A"],
            ["r6", "x1", "A"], ["r7", "x2", "C"]])
    _table(units / "screened_unassigned.tsv", ["unit_id", "reason", "genera"],
           [["x1", "genus_conflict", "A:3;B:2"], ["x2", "genus_conflict", "A:1;C:4"], ["x3", "other_genus", "B:"],
            ["x4", "unknown", ""]])
    paths = ds.subset_units(units, {"A", "B"}, tmp_path / "out")
    assert [r["unit_id"] for r in _read(paths["screened_representatives.tsv"])] == ["u1", "u2"]
    # records go with their unit, whatever their own genus: r2 (genus C) stays in u1, r5 (genus A) goes with u3
    assert [r["record_id"] for r in _read(paths["filtered.tsv"])] == ["r1", "r2", "r3", "r6"]
    # an unassigned unit stays only if all the genera of its records are kept; one with no genus goes
    assert [r["unit_id"] for r in _read(paths["screened_unassigned.tsv"])] == ["x1", "x3"]
    # without a subset the files are linked, not copied
    paths = ds.subset_units(units, None, tmp_path / "all")
    assert Path(paths["filtered.tsv"]).is_symlink() and len(_read(paths["screened_unassigned.tsv"])) == 4
    assert len(_read(paths["filtered.tsv"])) == 7


def test_memo_reuses_a_result_only_for_the_same_inputs(ds, tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "src"))
    import timegraft.build as build

    calls = []

    def fake(genus, ingroup, outgroup, stem_age, marker_cfg, known_crowns, timetree_crowns=None, rooting=None,
             rate=None, crown_basis=None):
        calls.append(genus)
        return f"tree of {genus} at {stem_age}", {"genus": genus, "stem_age": stem_age}

    monkeypatch.setattr(build, "build_genus", fake)
    stats = ds.install_memo(str(REPO), tmp_path / "memo")
    cfg, rooting = {"distance": {"method": "p"}}, {"rule": "tolerant"}
    args = ("A", [("u1", "ACGT"), ("u2", "ACGA")], [("o1", "TCGT")], 10.0, cfg, {}, {"A": 4.0}, rooting, None, {})
    first = build.build_genus(*args)
    assert build.build_genus(*args) == first and calls == ["A"] and stats == {"hits": 1, "misses": 1}
    # entries of the crown tables for other genera do not matter; the genus's own entry does
    assert build.build_genus(*args[:6], {"A": 4.0, "B": 9.0}, *args[7:]) == first and calls == ["A"]
    build.build_genus(*args[:6], {"A": 5.0}, *args[7:])
    # so do the stem age, the outgroup and the settings
    build.build_genus(*args[:3], 11.0, *args[4:])
    build.build_genus(args[0], args[1], [("o2", "TCGT")], *args[3:])
    build.build_genus(*args[:7], {"rule": "strict"}, *args[8:])
    assert len(calls) == 5 and stats["hits"] == 2

