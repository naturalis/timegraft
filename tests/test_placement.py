"""Matching and placement of ASVs."""

import shutil

import dendropy
import pytest

from timegraft.dating import is_ultrametric
from timegraft.placement import attach, classify, place_all, vsearch_hits

TREE = "((A:2,B:2):8,C:10);"


def load():
    return dendropy.Tree.get(data=TREE, schema="newick", rooting="force-rooted")


def depth(tree, label):
    node, path = tree.find_node_with_taxon_label(label), 0.0
    while node.parent_node is not None:
        path += node.edge.length
        node = node.parent_node
    return path


def test_attach_splits_the_edge_at_the_right_age():
    tree = load()
    attach(tree, "q", tree.find_node_with_taxon_label("A"), 1.0)
    assert is_ultrametric(tree)
    q = tree.find_node_with_taxon_label("q")
    assert {l.taxon.label for l in q.parent_node.leaf_iter()} == {"A", "q"}
    assert q.edge.length == pytest.approx(1.0)


def test_attach_walks_up_past_younger_nodes():
    tree = load()
    attach(tree, "q", tree.find_node_with_taxon_label("A"), 5.0)
    q = tree.find_node_with_taxon_label("q")
    assert {l.taxon.label for l in q.parent_node.leaf_iter()} == {"A", "B", "q"}
    assert depth(tree, "q") == pytest.approx(10.0)
    assert is_ultrametric(tree)


def test_attach_older_than_root_joins_root():
    tree = load()
    assert attach(tree, "q", tree.find_node_with_taxon_label("A"), 50.0) == pytest.approx(10.0)
    assert is_ultrametric(tree)


def test_classify_and_place():
    hits = {
        "m": [("rA1", 0.99), ("rA2", 0.98)],
        "tie": [("rA1", 0.99), ("rB1", 0.99)],
        "p": [("rC1", 0.90)],
    }
    record_unit = {"rA1": "A", "rA2": "A", "rB1": "B", "rC1": "C"}
    rows = classify(hits, ["m", "tie", "p", "none"], record_unit, {"A": "G1", "B": "G1", "C": "G2"},
                    {"A", "B", "C"}, match_identity=0.97)
    status = {r["asv_id"]: r["status"] for r in rows}
    assert status == {"m": "matched", "tie": "matched_ambiguous", "p": "placed", "none": "unplaced"}
    tree = load()
    out = {r["asv_id"]: r for r in place_all(tree, rows, height="rate", rate=0.01)}
    assert out["m"]["tip"] == "A"
    assert out["p"]["attach_age"] == "5.0000"  # (1 - 0.90) / (2 * 0.01)
    assert out["tie"]["attach_age"] == "2.0000"  # at the MRCA of A and B
    assert is_ultrametric(tree)


@pytest.mark.skipif(shutil.which("vsearch") is None, reason="vsearch not installed")
def test_vsearch_hits_fragment_against_full_length():
    ref = "ATGGCTAGCTAGGATCGATCGTTAGCTAGCTAGGCTAGCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCTAGCTAGCATGCATCG" * 2
    hits = vsearch_hits([("frag", ref[30:150])], [("BOLD:X", ref)], min_identity=0.9)
    assert hits["frag"][0] == ("BOLD:X", 1.0)
