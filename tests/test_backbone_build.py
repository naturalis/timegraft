"""Backbone assembly with mocked OToL and TimeTree responses."""

import json

import dendropy
import pytest

from timegraft.backbone_build import date_topology
from timegraft.dating import is_ultrametric, root_age
from timegraft.opentree import OpenTree, genus_topology
from timegraft.timetree import TimeTree
from timegraft.webapi import ApiError


def fake_timetree(ages, missing=()):
    """TimeTree stand-in: ``ages`` maps frozensets of ids to the age of their MRCA."""

    def fetch(url, payload=None):
        ids = [int(i) for i in url.split("/mrca/id/")[1].split("/")[0].split("+")]
        found = [i for i in ids if i not in missing]
        age = ages.get(frozenset(found))
        return {"found_ids": found, "precomputed_age": age}

    return TimeTree(fetch=fetch)


def test_split_age_needs_taxa_from_two_sides():
    tt = fake_timetree({frozenset({1, 2}): 50.0, frozenset({2}): None}, missing={1})
    assert tt.split_age([[1], [2]]) is None
    tt = fake_timetree({frozenset({1, 2}): 50.0})
    assert tt.split_age([[1], [2]]) == 50.0


def test_date_topology_calibrates_and_interpolates():
    tree = dendropy.Tree.get(data="(((A,B),C),(D,E));", schema="newick", rooting="force-rooted")
    ncbi = {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5}
    # root and the (A,B,C) node are known; (A,B) and (D,E) are not
    ages = {
        frozenset({1, 3, 4, 5}): 100.0,
        frozenset({1, 2, 4, 5}): 100.0,
        frozenset({1, 2, 3}): 60.0,
    }
    tree, report = date_topology(tree, ncbi, fake_timetree(ages))
    assert is_ultrametric(tree)
    assert root_age(tree) == pytest.approx(100.0)
    by_genera = {r["example_genera"]: r for r in report}
    assert by_genera["A;B"]["source"] == "interpolated"
    assert float(by_genera["A;B"]["age"]) == pytest.approx(30.0)
    assert float(by_genera["D;E"]["age"]) == pytest.approx(50.0)


def test_date_topology_drops_ages_older_than_ancestor():
    tree = dendropy.Tree.get(data="((A,B),C);", schema="newick", rooting="force-rooted")
    ncbi = {"A": 1, "B": 2, "C": 3}
    ages = {frozenset({1, 2, 3}): 40.0, frozenset({1, 3}): 40.0, frozenset({1, 2}): 55.0}
    tree, report = date_topology(tree, ncbi, fake_timetree(ages))
    sources = {r["example_genera"]: r["source"] for r in report}
    assert sources["A;B"] == "interpolated_after_drop"
    assert is_ultrametric(tree)


def test_root_age_required_when_unknown():
    tree = dendropy.Tree.get(data="(A,B);", schema="newick", rooting="force-rooted")
    with pytest.raises(ValueError):
        date_topology(tree, {"A": 1, "B": 2}, fake_timetree({}))
    tree, _ = date_topology(tree, {"A": 1, "B": 2}, fake_timetree({}), root_age=12.0)
    assert root_age(tree) == pytest.approx(12.0)


def test_match_genera_picks_genus_rank():
    def fetch(url, payload):
        return {
            "results": [
                {
                    "name": "Mytilus",
                    "matches": [
                        {"is_synonym": False, "is_approximate_match": False,
                         "taxon": {"ott_id": 11, "rank": "genus", "tax_sources": ["worms:1", "ncbi:6548"]}},
                    ],
                },
                {
                    "name": "Homonym",
                    "matches": [
                        {"is_synonym": False, "is_approximate_match": False, "taxon": {"ott_id": 21, "rank": "genus"}},
                        {"is_synonym": False, "is_approximate_match": False, "taxon": {"ott_id": 22, "rank": "genus"}},
                    ],
                },
            ]
        }

    result = OpenTree(fetch=fetch).match_genera(["Mytilus", "Homonym", "Nothing"], context="Animals")
    assert result["Mytilus"] == {"ott_id": 11, "ncbi_id": 6548, "status": "matched"}
    assert result["Homonym"]["status"] == "ambiguous"
    assert result["Nothing"]["status"] == "unmatched"


def test_induced_subtree_retries_without_unknown_ids():
    calls = []

    def fetch(url, payload):
        calls.append(list(payload["ott_ids"]))
        if 3 in payload["ott_ids"]:
            raise ApiError(400, json.dumps({"message": "unknown", "unknown": {"ott3": "pruned_ott_id"}}))
        return {"newick": "((ott1,ott2)mrcaott1ott2,ott4)ott99;", "broken": {"ott5": "mrcaott1ott2"}}

    newick, missing, broken = OpenTree(fetch=fetch).induced_subtree([1, 2, 3, 4, 5])
    assert missing == {3: "pruned_ott_id"}
    assert broken == {5: "mrcaott1ott2"}
    tree = genus_topology(newick, {1: "A", 2: "B", 4: "D", 5: "E"}, broken)
    assert {l.taxon.label for l in tree.leaf_node_iter()} == {"A", "B", "D", "E"}
    # E was attached at the MRCA of A and B
    e = tree.find_node_with_taxon_label("E")
    assert {l.taxon.label for l in e.parent_node.leaf_iter()} == {"A", "B", "E"}


def test_equal_ages_mean_unresolved_not_inconsistent():
    tree = dendropy.Tree.get(data="((A,B),C);", schema="newick", rooting="force-rooted")
    ncbi = {"A": 1, "B": 2, "C": 3}
    ages = {frozenset({1, 2, 3}): 40.0, frozenset({1, 3}): 40.0, frozenset({1, 2}): 40.0}
    tree, report = date_topology(tree, ncbi, fake_timetree(ages))
    row = next(r for r in report if r["example_genera"] == "A;B")
    assert row["source"] == "interpolated_unresolved_in_timetree"
    assert row["conflict_ratio"] == ""


def test_direct_age_outranks_higher_taxa_ancestor():
    # (A,B) is dated directly at 30; its parent ((A,B),C) only via enclosing taxa, at 20
    tree = dendropy.Tree.get(data="(((A,B),C),D);", schema="newick", rooting="force-rooted")
    ncbi = {"A": 1, "B": 2, "C": 3, "D": 4}
    ages = {frozenset({1, 2, 4}): 100.0, frozenset({1, 3, 4}): 100.0, frozenset({1, 2}): 30.0,
            frozenset({11, 13}): 20.0}
    higher = {"A": [(101, 11)], "B": [(101, 11)], "C": [(103, 13)], "D": [(104, 14)]}
    tree, report = date_topology(tree, ncbi, fake_timetree(ages, missing={3}), higher=higher)
    by = {r["example_genera"]: r for r in report}
    assert by["A;B"]["source"] == "timetree"
    assert by["A;B;C"]["source"] == "interpolated_after_drop"
    assert float(by["A;B;C"]["conflict_ratio"]) == pytest.approx(1.5)
    assert is_ultrametric(tree)
