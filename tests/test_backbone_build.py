"""Backbone assembly with mocked OToL and TimeTree responses."""

import json

import dendropy
import pytest

import math

from timegraft.backbone_build import date_topology, interpolate, yule_height
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
    # without clade richness, pure birth falls back to even spacing
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
    tree, report = date_topology(tree, ncbi, fake_timetree(ages), method="precedence")
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
    tree, report = date_topology(tree, ncbi, fake_timetree(ages, missing={3}), higher=higher, method="precedence")
    by = {r["example_genera"]: r for r in report}
    assert by["A;B"]["source"] == "timetree"
    assert by["A;B;C"]["source"] == "interpolated_after_drop"
    assert float(by["A;B;C"]["conflict_ratio"]) == pytest.approx(1.5)
    assert is_ultrametric(tree)


def test_yule_height():
    assert yule_height(1) == 0.0
    assert yule_height(2) == pytest.approx(0.5)
    assert yule_height(3) == pytest.approx(0.5 + 1 / 3)
    # the asymptotic form joins the exact sum
    exact = sum(1 / k for k in range(2, 1002))
    assert yule_height(1001) == pytest.approx(exact, abs=1e-12)
    assert yule_height(10**6) == pytest.approx(math.log(10**6) + 0.5772156649 - 1, abs=1e-6)


def _chain():
    # root dated at 100; ((A,B),C) and (A,B) undated; D outside
    tree = dendropy.Tree.get(data="(((A,B),C),D);", schema="newick", rooting="force-rooted")
    abc = tree.find_node_with_taxon_label("C").parent_node
    ab = tree.find_node_with_taxon_label("A").parent_node
    return tree, abc, ab


def test_even_spacing_matches_bladj():
    tree, abc, ab = _chain()
    ages = interpolate(tree, {tree.seed_node: 90.0}, method="even")
    assert ages[abc] == pytest.approx(60.0)
    assert ages[ab] == pytest.approx(30.0)


def test_pure_birth_places_nodes_by_clade_richness():
    tree, abc, ab = _chain()
    richness = {tree.seed_node: 1000, abc: 100, ab: 10}
    ages = interpolate(tree, {tree.seed_node: 90.0}, method="pure_birth", richness=richness)
    assert ages[abc] == pytest.approx(90.0 * yule_height(100) / yule_height(1000))
    assert ages[ab] == pytest.approx(90.0 * yule_height(10) / yule_height(1000))
    # the same clade sizes give the same ages whatever else is sampled: drop C
    tree2 = dendropy.Tree.get(data="((A,B),D);", schema="newick", rooting="force-rooted")
    ab2 = tree2.find_node_with_taxon_label("A").parent_node
    again = interpolate(tree2, {tree2.seed_node: 90.0}, method="pure_birth", richness={tree2.seed_node: 1000, ab2: 10})
    assert again[ab2] == pytest.approx(ages[ab])


def test_pure_birth_between_dated_nodes():
    tree, abc, ab = _chain()
    richness = {tree.seed_node: 1000, abc: 100, ab: 10}
    ages = interpolate(tree, {tree.seed_node: 90.0, ab: 20.0}, method="pure_birth", richness=richness)
    h = yule_height
    assert ages[abc] == pytest.approx(20.0 + 70.0 * (h(100) - h(10)) / (h(1000) - h(10)))


def test_pure_birth_falls_back_to_even_where_richness_is_not_ordered():
    tree, abc, ab = _chain()
    ages = interpolate(tree, {tree.seed_node: 90.0}, method="pure_birth", richness={tree.seed_node: 50, abc: 80, ab: 10})
    assert ages[abc] == pytest.approx(60.0)
    assert 0 < ages[ab] < ages[abc]


def test_date_topology_uses_richness_and_reports_it():
    tree = dendropy.Tree.get(data="(((A,B),C),(D,E));", schema="newick", rooting="force-rooted")
    ncbi = {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5}
    ages = {frozenset({1, 3, 4, 5}): 100.0, frozenset({1, 2, 4, 5}): 100.0, frozenset({1, 2, 3}): 60.0}
    ab = tree.find_node_with_taxon_label("A").parent_node
    de = tree.find_node_with_taxon_label("D").parent_node
    abc = ab.parent_node
    richness = {tree.seed_node: 500, abc: 50, ab: 5, de: 20}
    tree, report = date_topology(tree, ncbi, fake_timetree(ages), richness=richness)
    by_genera = {r["example_genera"]: r for r in report}
    assert float(by_genera["A;B"]["age"]) == pytest.approx(60.0 * yule_height(5) / yule_height(50), abs=1e-3)
    assert float(by_genera["D;E"]["age"]) == pytest.approx(100.0 * yule_height(20) / yule_height(500), abs=1e-3)
    assert by_genera["A;B"]["clade_genera"] == 5
    assert is_ultrametric(tree)


def test_unknown_interpolation_method():
    tree, _, _ = _chain()
    with pytest.raises(ValueError):
        interpolate(tree, {tree.seed_node: 1.0}, method="bladj")
