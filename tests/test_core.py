import math

import dendropy
import numpy as np
import pytest

from timegraft import distances, filters
from timegraft.backbone import nearest_genera, stem_ages
from timegraft.bionj import bionj
from timegraft.dating import crown_age, is_ultrametric, mpl_ultrametric, root_age
from timegraft.graft import graft
from timegraft.rooting import root_on_outgroup
from timegraft.units import representatives


def patristic(tree, a, b):
    pdm = tree.phylogenetic_distance_matrix()
    ta = tree.taxon_namespace.get_taxon(label=a)
    tb = tree.taxon_namespace.get_taxon(label=b)
    return pdm(ta, tb)


# distances ---------------------------------------------------------------

def test_k2p_identical_is_zero():
    assert distances.k2p("ACGTACGT", "ACGTACGT") == 0.0


def test_k2p_known_value():
    # 1 transition in 10 sites: P = 0.1, Q = 0
    d = distances.k2p("AAAAAAAAAA", "GAAAAAAAAA")
    assert d == pytest.approx(-0.5 * math.log(0.8))


def test_k2p_skips_gaps_and_ambiguity():
    assert distances.k2p("AC-TN", "ACGTA") == 0.0


def test_k2p_saturated_returns_cap():
    assert distances.k2p("AAAA", "CCCC", max_distance=2.0) == 2.0


def test_kmer_identical_and_gap_insensitive():
    seq = "ACGTTGCAAGGCTTACGATCGATCGTAGCTAGCTAGG"
    assert distances.kmer(seq, seq, k=5) == 0.0
    assert distances.kmer(seq, "--" + seq, k=5) == 0.0


# BioNJ and rooting --------------------------------------------------------

def additive_matrix():
    # tree ((A:1,B:2):1,(C:3,D:1)) unrooted, internal edge 1
    labels = ["A", "B", "C", "D"]
    d = np.array(
        [
            [0, 3, 5, 3],
            [3, 0, 6, 4],
            [5, 6, 0, 4],
            [3, 4, 4, 0],
        ],
        dtype=float,
    )
    return labels, d


def test_bionj_recovers_additive_tree():
    labels, d = additive_matrix()
    tree = bionj(labels, d)
    for i, a in enumerate(labels):
        for j, b in enumerate(labels):
            if i < j:
                assert patristic(tree, a, b) == pytest.approx(d[i, j])


def test_root_on_outgroup_returns_ingroup_clade():
    labels, d = additive_matrix()
    tree = bionj(labels, d)
    rooted, how = root_on_outgroup(tree, ["A", "B"], ["C", "D"])
    assert how == "outgroup"
    assert {l.taxon.label for l in rooted.leaf_node_iter()} == {"A", "B"}


def test_root_without_outgroup_is_midpoint():
    labels, d = additive_matrix()
    tree = bionj(labels, d)
    rooted, how = root_on_outgroup(tree, labels, [])
    assert how == "midpoint"
    assert len(rooted.leaf_nodes()) == 4
    # the longest path is B-C (6), so both lie 3 from the root
    depth = {}
    for leaf in rooted.leaf_node_iter():
        node, path = leaf, 0.0
        while node.parent_node is not None:
            path += node.edge.length
            node = node.parent_node
        depth[leaf.taxon.label] = path
    assert depth["B"] == pytest.approx(3.0)
    assert depth["C"] == pytest.approx(3.0)


# dating ----------------------------------------------------------------

def test_mpl_ultrametric_root_age():
    tree = dendropy.Tree.get(data="((A:0.1,B:0.3):0.2,C:0.25);", schema="newick", rooting="force-rooted")
    mpl_ultrametric(tree, 12.0)
    assert is_ultrametric(tree)
    assert root_age(tree) == pytest.approx(12.0)
    assert all(e.length >= 0 for e in tree.postorder_edge_iter() if e.length is not None)


def test_crown_age_order_of_preference():
    tree = dendropy.Tree.get(data="(A:0.05,B:0.05);", schema="newick", rooting="force-rooted")
    assert crown_age("G", 20, tree, {"G": 8.0}).source == "table"
    by_rate = crown_age("G", 20, tree, {}, method="rate", rate=0.01)
    assert by_rate.source == "rate" and by_rate.age == pytest.approx(5.0)
    capped = crown_age("G", 4, tree, {}, method="rate", rate=0.01, cap_fraction=0.9)
    assert capped.source == "rate_capped" and capped.age == pytest.approx(3.6)
    frac = crown_age("G", 20, tree, {}, method="fraction", fraction=0.5)
    assert frac.source == "fraction" and frac.age == pytest.approx(10.0)


# filters ----------------------------------------------------------------

def test_coi_filter_rejects_stops_in_all_frames():
    clean = "ATG" * 80
    assert filters.check_coi(clean) is None
    dirty = ("TAA" + "TAG" + "ATAAG") * 30
    assert filters.check_coi(dirty) == "stop_codon"
    assert filters.check_coi("ATG" * 10) == "too_short"


def test_tga_is_a_stop_only_in_the_standard_code():
    # every frame of TGACTGAC... contains TGA, and no TAA or TAG
    seq = "TGAC" * 75
    assert filters.has_open_frame(seq, genetic_code=5)
    assert not filters.has_open_frame(seq, genetic_code=1)


# units -------------------------------------------------------------------

def test_representatives_assign_and_flag():
    recs = [
        {"record_id": "r1", "unit_id": "U1", "genus": "Alpha", "sequence": "ACGTACGTTAGC" * 5},
        {"record_id": "r2", "unit_id": "U1", "genus": "Alpha", "sequence": "ACGTACGTTAGC" * 5},
        {"record_id": "r3", "unit_id": "U2", "genus": "", "sequence": "ACGTACGTTAGC" * 5},
        {"record_id": "r4", "unit_id": "U3", "genus": "Alpha", "sequence": "ACGTACGTTAGC" * 5},
        {"record_id": "r5", "unit_id": "U3", "genus": "Beta", "sequence": "ACGTACGTTAGC" * 5},
    ]
    reps, unassigned = representatives(recs, k=4)
    assert [r["unit_id"] for r in reps] == ["U1"]
    reasons = {u["unit_id"]: u["reason"] for u in unassigned}
    assert reasons == {"U2": "no_genus", "U3": "genus_conflict"}


# backbone and grafting --------------------------------------------------

BACKBONE = "((Alpha:20,Beta:20):30,(Gamma:35,(Delta:15,Epsilon:15):20):15);"


def test_nearest_genera():
    tree = dendropy.Tree.get(data=BACKBONE, schema="newick", rooting="force-rooted")
    assert nearest_genera(tree, "Alpha", {"Beta", "Gamma"}) == ["Beta"]
    assert nearest_genera(tree, "Alpha", {"Gamma", "Delta"}) == ["Delta", "Gamma"]
    assert nearest_genera(tree, "Delta", {"Alpha"}) == ["Alpha"]


def test_graft_keeps_tree_ultrametric():
    backbone = dendropy.Tree.get(data=BACKBONE, schema="newick", rooting="force-rooted")
    assert stem_ages(backbone)["Alpha"] == 20
    sub = "(('BOLD:AAA0001':3,'BOLD:AAA0002':3):5,'BOLD:AAA0003':8);"
    tree, status = graft(backbone, {"Alpha": (sub, 8.0)}, {"Delta": "BOLD:AAD0001"})
    labels = {l.taxon.label for l in tree.leaf_node_iter()}
    assert labels == {"BOLD:AAA0001", "BOLD:AAA0002", "BOLD:AAA0003", "BOLD:AAD0001"}
    assert status == {
        "Alpha": "subtree",
        "Beta": "pruned",
        "Gamma": "pruned",
        "Delta": "single_unit",
        "Epsilon": "pruned",
    }
    assert is_ultrametric(tree)
    assert root_age(tree) == pytest.approx(50.0)
