"""What made each backbone node age (#84): the new columns of node_ages.tsv, and the comparison of two runs."""

import csv
import importlib.util
from pathlib import Path

import dendropy
import pytest

from timegraft.backbone_build import REPORT_COLUMNS, date_topology, interpolate
from timegraft.timetree import TimeTree

COMPARE = Path(__file__).resolve().parents[1] / "data/iberia/consistency/benchmark/compare_backbones.py"


def fake_timetree(ages, missing=(), adjusted=None):
    def fetch(url, payload=None):
        ids = [int(i) for i in url.split("/mrca/id/")[1].split("/")[0].split("+")]
        found = [i for i in ids if i not in missing]
        return {"found_ids": found, "precomputed_age": ages.get(frozenset(found)),
                "adjusted_age": (adjusted or {}).get(frozenset(found))}

    return TimeTree(fetch=fetch)


def by_genera(report):
    return {r["example_genera"]: r for r in report}


def test_interpolate_records_its_anchors():
    tree = dendropy.Tree.get(data="(((A,B),C),D);", schema="newick", rooting="force-rooted")
    abc = tree.find_node_with_taxon_label("C").parent_node
    ab = tree.find_node_with_taxon_label("A").parent_node
    anchors = {}
    interpolate(tree, {tree.seed_node: 90.0, ab: 20.0}, method="even", anchors=anchors)
    assert anchors == {abc: (tree.seed_node, ab)}
    anchors = {}
    interpolate(tree, {tree.seed_node: 90.0}, method="even", anchors=anchors)
    assert set(anchors) == {abc, ab}
    assert anchors[ab][0] is tree.seed_node and anchors[ab][1].is_leaf()  # no dated descendant: the tips
    assert anchors[abc][0] is tree.seed_node


def test_age_source_interpolated_from_and_query_genera():
    tree = dendropy.Tree.get(data="(((A,B),C),(D,E));", schema="newick", rooting="force-rooted")
    ncbi = {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5}
    ages = {frozenset({1, 3, 4, 5}): 100.0, frozenset({1, 2, 4, 5}): 100.0, frozenset({1, 2, 3}): 60.0}
    _, report = date_topology(tree, ncbi, fake_timetree(ages))
    assert set(REPORT_COLUMNS) >= {"age_source", "interpolated_from", "query_genera"}
    by = by_genera(report)
    names = {r["example_genera"]: r["node"] for r in report}
    root = report[0]
    assert root["age_source"] == "root" and root["interpolated_from"] == ""
    assert by["A;B;C"]["age_source"] == "timetree" and by["A;B;C"]["interpolated_from"] == ""
    assert set(by["A;B;C"]["query_genera"].split(";")) <= {"A", "B", "C"} and "C" in by["A;B;C"]["query_genera"]
    # (A,B) hangs below the dated (A,B,C) and above the tips; (D,E) below the root
    assert by["A;B"]["age_source"] == "interpolated"
    assert by["A;B"]["interpolated_from"] == f"{names['A;B;C']};tips"
    assert by["D;E"]["interpolated_from"] == f"{root['node']};tips"
    assert by["A;B"]["query_genera"] == ""  # TimeTree did not date it


def test_fit_that_moves_a_calibration():
    # (A,B) is given 55 Ma under a root of 40: the fit moves it, the precedence rule drops it
    tree = dendropy.Tree.get(data="((A,B),C);", schema="newick", rooting="force-rooted")
    ncbi = {"A": 1, "B": 2, "C": 3}
    ages = {frozenset({1, 2, 3}): 40.0, frozenset({1, 3}): 40.0, frozenset({1, 2}): 55.0}
    _, fitted = date_topology(tree, ncbi, fake_timetree(ages))
    assert by_genera(fitted)["A;B"]["age_source"] == "fit"
    tree = dendropy.Tree.get(data="((A,B),C);", schema="newick", rooting="force-rooted")
    _, kept = date_topology(tree, ncbi, fake_timetree(ages), method="precedence")
    row = by_genera(kept)["A;B"]
    assert row["age_source"] == "interpolated" and row["interpolated_from"] == "n0;tips"
    assert row["query_genera"] == "A;B"  # it was queried, but its age was dropped


def test_query_with_enclosing_taxa_lists_no_genera():
    tree = dendropy.Tree.get(data="(((A,B),C),D);", schema="newick", rooting="force-rooted")
    ncbi = {"A": 1, "B": 2, "C": 3, "D": 4}
    ages = {frozenset({1, 2, 4}): 100.0, frozenset({1, 3, 4}): 100.0, frozenset({1, 2}): 30.0,
            frozenset({11, 13}): 20.0}
    higher = {"A": [(101, 11)], "B": [(101, 11)], "C": [(103, 13)], "D": [(104, 14)]}
    _, report = date_topology(tree, ncbi, fake_timetree(ages, missing={3}), higher=higher, method="precedence")
    row = by_genera(report)["A;B;C"]
    assert row["timetree_attempt"] == "higher_taxa" and row["query_ids"] and row["query_genera"] == ""


def _write_backbone(folder, ncbi, ages, adjusted=None, timetree_time=None, columns=REPORT_COLUMNS):
    folder.mkdir()
    tree = dendropy.Tree.get(data="(((A,B),C),(D,E));", schema="newick", rooting="force-rooted")
    tree, report = date_topology(tree, ncbi, fake_timetree(ages, adjusted=adjusted), interpolation="even",
                                 timetree_time=timetree_time)
    tree.write(path=str(folder / "backbone.tre"), schema="newick", suppress_rooting=True)
    with open(folder / "node_ages.tsv", "w", newline="") as handle:
        w = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(report)
    with open(folder / "genera.tsv", "w", newline="") as handle:
        w = csv.DictWriter(handle, ["genus", "join_ncbi_id", "anchor_node", "placement_action"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows({"genus": g, "join_ncbi_id": str(i or ""), "anchor_node": f"ott{n}", "placement_action": ""}
                    for n, (g, i) in enumerate(ncbi.items()))
    return report


def test_compare_backbones_traces_a_changed_identifier(tmp_path):
    spec = importlib.util.spec_from_file_location("compare_backbones", COMPARE)
    cb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cb)
    ages = {frozenset({1, 2, 4, 5}): 100.0, frozenset({1, 3, 4, 5}): 100.0, frozenset({1, 2, 3}): 60.0}
    # in A, C has no identifier, so ((A,B),C) cannot be queried; in B it has one, and gets 60 Ma
    _write_backbone(tmp_path / "a", {"A": 1, "B": 2, "C": 0, "D": 4, "E": 5}, ages)
    _write_backbone(tmp_path / "b", {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5}, ages)
    a, b = cb.Backbone(tmp_path / "a"), cb.Backbone(tmp_path / "b")
    causes, chain = cb.explain(a, b, b.stem("C"))
    assert set(causes) == {"calibration_gained", "identifier"}
    assert any("C - to 3" in step for step in chain)
    # (A,B) is interpolated in both; its parent moved, and the chain leads there
    causes, chain = cb.explain(a, b, b.stem("A"))
    assert set(causes) == {"calibration_gained", "identifier"} and len(chain) >= 2
    # the whole script: every changed stem traced
    import sys
    sys.argv = ["compare_backbones.py", str(tmp_path / "a"), str(tmp_path / "b"), str(tmp_path / "out")]
    cb.main()
    rows = list(csv.DictReader(open(tmp_path / "out/stem_causes.tsv"), delimiter="\t"))
    assert {r["genus"] for r in rows} >= {"A", "B", "C"} and all(r["traced"] == "yes" for r in rows)
    assert all(r["causes"] == "identifier;calibration_gained" for r in rows)


def test_compare_backbones_finds_nothing_in_identical_runs(tmp_path):
    spec = importlib.util.spec_from_file_location("compare_backbones", COMPARE)
    cb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cb)
    ages = {frozenset({1, 2, 4, 5}): 100.0, frozenset({1, 3, 4, 5}): 100.0, frozenset({1, 2, 3}): 60.0}
    ncbi = {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5}
    _write_backbone(tmp_path / "a", ncbi, ages)
    _write_backbone(tmp_path / "b", ncbi, ages)
    a, b = cb.Backbone(tmp_path / "a"), cb.Backbone(tmp_path / "b")
    assert all(cb.explain(a, b, b.stem(g)) == ([], []) for g in "ABCDE")


def test_compare_backbones_needs_the_new_columns(tmp_path):
    spec = importlib.util.spec_from_file_location("compare_backbones", COMPARE)
    cb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cb)
    _write_backbone(tmp_path / "a", {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5},
                    {frozenset({1, 2, 4, 5}): 100.0, frozenset({1, 3, 4, 5}): 100.0})
    rows = list(csv.DictReader(open(tmp_path / "a/node_ages.tsv"), delimiter="\t"))
    with open(tmp_path / "a/node_ages.tsv", "w", newline="") as handle:
        cols = [c for c in rows[0] if c != "age_source"]
        w = csv.DictWriter(handle, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    with pytest.raises(SystemExit):
        cb.Backbone(tmp_path / "a")


def _row(**kw):
    row = {"query_ids": "", "query_genera": "", "timetree_attempt": "", "timetree_age": "", "source": "interpolated",
           "precedence_source": ""}
    row.update(kw)
    return row


def test_a_gained_calibration_states_why(tmp_path):
    spec = importlib.util.spec_from_file_location("compare_backbones", COMPARE)
    cb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cb)
    ages = {frozenset({1, 2, 4, 5}): 100.0, frozenset({1, 3, 4, 5}): 100.0, frozenset({1, 2, 3}): 60.0}
    _write_backbone(tmp_path / "a", {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5}, ages)
    a = b = cb.Backbone(tmp_path / "a")
    node = b.stem("C")
    cal = _row(query_ids="1+2+3", query_genera="A;B;C", timetree_attempt="0", timetree_age="60", source="timetree")
    why = lambda without: cb._why_calibrated(a, b, cal, without, node, "A", a)[0]
    # TimeTree had no age for the other run's picks, its identifiers unchanged
    assert why(_row()) == ["no_answer"]
    # the other run asked about other species
    assert why(_row(query_ids="1+3", query_genera="A;C", timetree_attempt="0", timetree_age="60",
                    source="interpolated_unresolved_in_timetree")) == ["picks", "unresolved"]
    # or with enclosing taxa
    assert why(_row(query_ids="11+13", timetree_attempt="higher_taxa", timetree_age="70")) == ["higher_taxa"]
    # the same query, another answer from TimeTree
    assert why(_row(query_ids="1+2+3", query_genera="A;B;C", timetree_attempt="0", timetree_age="")) == \
        ["timetree_answer"]
    # the same calibration, dropped by the precedence rule
    assert why(_row(query_ids="1+2+3", query_genera="A;B;C", timetree_attempt="0", timetree_age="60",
                    source="interpolated_after_drop")) == ["fit"]


def test_a_bare_calibration_change_is_not_traced():
    spec = importlib.util.spec_from_file_location("compare_backbones", COMPARE)
    cb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cb)
    assert not cb.traced(["calibration_gained"]) and not cb.traced(["calibration_lost", "anchors"])
    assert cb.traced(["calibration_gained", "no_answer"]) and cb.traced(["richness"])


def test_a_change_of_time_is_traced_against_an_older_report(tmp_path):
    # A is a report before format 7 (no calibration_age: it calibrated with timetree_age); B calibrates ((A,B),C)
    # with TimeTree's adjusted time (#108). The stem change is the choice of time, not the fit
    spec = importlib.util.spec_from_file_location("compare_backbones", COMPARE)
    cb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cb)
    ages = {frozenset({1, 2, 4, 5}): 100.0, frozenset({1, 3, 4, 5}): 100.0, frozenset({1, 2, 3}): 30.0}
    adjusted = {frozenset({1, 2, 3}): 60.0}
    ncbi = {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5}
    old = [c for c in REPORT_COLUMNS if c not in ("calibration_age", "timetree_adjusted_age")]
    _write_backbone(tmp_path / "a", ncbi, ages, adjusted, {"use": "raw"}, columns=old)
    _write_backbone(tmp_path / "b", ncbi, ages, adjusted, {"use": "adjusted"})
    a, b = cb.Backbone(tmp_path / "a"), cb.Backbone(tmp_path / "b")
    causes, chain = cb.explain(a, b, b.stem("C"))
    assert "timetree_time" in causes and "fit" not in causes
    assert any("calibrated at 30" in step and "60" in step for step in chain)
