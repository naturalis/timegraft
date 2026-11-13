"""Re-querying conflicting nodes with disjoint picks (#53)."""

import dendropy
import pytest

from timegraft import requery as rq
from timegraft.backbone_build import age_fit_summary, date_topology
from timegraft.timetree import TimeTree
from timegraft.webapi import JsonCache


def test_disjoint_side():
    ncbi = {"A": [1, 11, 12, 13, 14, 15, 16, 17, 18], "B": [2, 21], "C": [3], "D": None}
    # unused genera first
    assert rq.disjoint_side(["A", "B", "C", "D"], ncbi, {1, 11, 12, 13}, per_side=2) == [2, 21, 3]
    # every genus used: other species of the same genera
    assert rq.disjoint_side(["A"], ncbi, {1, 11, 12, 13}) == [14, 15, 16, 17]
    assert rq.disjoint_side(["A"], ncbi, {1, 11, 12, 13}, attempt=1) == [18]
    assert rq.disjoint_side(["C"], ncbi, {3}) == []
    for attempt in range(3):
        assert not set(rq.disjoint_side(["A", "B", "C"], ncbi, {1, 11, 12, 13}, attempt=attempt)) & {1, 11, 12, 13}


def test_triggers():
    out = rq.triggers({"a": 2.5, "b": -0.4, "c": -2.1}, {"a"}, {"d": 1.6, "e": 1.2}, {"d": "b", "e": "f"})
    assert out == {"a": "z;outlier", "c": "z", "d": "conflict", "b": "conflict"}


@pytest.mark.parametrize("first, again, reference, expected", [
    (54.9, 156.8, 168.6, ("under_spanning", "moved_older")),   # the crab ancestor of the analysis
    (721.5, 720.0, 360.0, ("topology", "stable_old")),        # Onchidella's Bilateria age, asked again
    (19.2, 19.2, 29.6, ("topology", "stable_young")),
    (300.0, 120.0, 110.0, ("pick_specific", "moved_younger")),
    (156.8, 156.8, 166.9, ("consistent", "close_to_reference")),
    (50.0, 20.0, 100.0, ("unresolved", "contrary")),
    (50.0, None, 100.0, ("unresolved", "no_answer")),
])
def test_classify(first, again, reference, expected):
    assert rq.classify(first, again, reference, 4, 4)[:2] == expected


def test_confidence():
    assert rq.classify(54.9, 156.8, 168.6, 4, 6)[2] == "high"
    assert rq.classify(54.9, 156.8, 168.6, 1, 1)[2] == "medium"  # decisive move, thin evidence
    assert rq.classify(54.9, 70.0, 168.6, 1, 1)[2] == "low"  # small move, thin evidence
    assert rq.classify(19.2, 22.5, 29.6, 1, 1)[2] == "low"


def _fake(ages, calls=None):
    def fetch(url, payload=None):
        ids = [int(i) for i in url.split("/mrca/id/")[1].split("/")[0].split("+")]
        if calls is not None:
            calls.append(tuple(ids))
        for key, age in ages.items():
            if set(ids) == set(key):
                return {"found_ids": ids, "precomputed_age": age, "time_estimates": "{" + f"{age},{age}" + "}"}
        return {"found_ids": ids, "precomputed_age": None}

    return fetch


# ((((A1,A2),(B1,B2)),D),C): the first query of node AB picks A1 and B1 and dates a shallower node at 20 Ma
# (under-spanning), younger than both its children (30 and 25 Ma); asked with A2 and B2, TimeTree dates it
# at 90 Ma.
NEWICK = "((((A1,A2),(B1,B2)),D),C);"
NCBI = {"A1": [1], "A2": [2], "B1": [3], "B2": [4], "C": [5], "D": [6]}
AGES = {(1, 3): 20.0, (2, 4): 90.0, (1, 2): 30.0, (3, 4): 25.0, (1, 5): 100.0, (1, 6): 95.0}


def _run(mode, cache=None, calls=None):
    tree = dendropy.Tree.get(data=NEWICK, schema="newick", rooting="force-rooted")
    tt = TimeTree(cache or JsonCache(None), fetch=_fake(AGES, calls))
    _, report = date_topology(tree, NCBI, tt, per_side=1, per_genus=0, requery={"mode": mode, "z_trigger": 1.0})
    return tree, {(r["example_genera"], r["n_genera"]): r for r in report}, report


def test_requery_classifies_and_can_apply():
    _, by, report = _run("classify")
    ab = by[("A1;A2;B1", 4)]
    assert ab["timetree_age"] == "20.0000" and ab["requery_ids"] == "2+4"
    assert ab["requery_age"] == "90.0000" and ab["conflict_class"] == "under_spanning"
    assert ab["conflict_detail"] == "moved_older" and ab["requery_applied"] == ""
    assert "z" in ab["requery_trigger"]
    untouched = [r for r in report if not r["requery_trigger"]]
    assert all(r["conflict_class"] == "" for r in untouched)
    classified_age = float(ab["age"])
    _, by, report = _run("apply")
    ab = by[("A1;A2;B1", 4)]
    assert ab["requery_applied"] == "yes" and float(ab["age"]) > classified_age
    assert float(ab["age"]) == pytest.approx(90.0, rel=0.1)
    summary = age_fit_summary(report)
    assert summary["final_against_requery"]["weighted_rms"] < summary["final"]["weighted_rms"]
    _, by, _ = _run("off")
    assert by[("A1;A2;B1", 4)]["conflict_class"] == "" and by[("A1;A2;B1", 4)]["requery_ids"] == ""


def test_requery_is_reproducible_from_the_cache():
    cache, calls = JsonCache(None), []
    _, first, _ = _run("classify", cache, calls)
    n = len(calls)
    _, second, _ = _run("classify", cache, calls)
    assert len(calls) == n  # nothing asked again
    assert {k: (r["conflict_class"], r["requery_age"]) for k, r in first.items()} == {
        k: (r["conflict_class"], r["requery_age"]) for k, r in second.items()}


def test_requery_bad_mode():
    with pytest.raises(ValueError):
        _run("sometimes")
    assert _run(False)[1][("A1;A2;B1", 4)]["conflict_class"] == ""  # YAML off
