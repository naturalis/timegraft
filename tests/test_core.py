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


@pytest.mark.skipif(__import__("shutil").which("lsd2") is None, reason="lsd2 not installed")
def test_lsd2_dating_fixes_root_and_tips():
    from timegraft.dating import lsd2_ultrametric

    tree = dendropy.Tree.get(
        data="(('BOLD:A1':0.05,'BOLD:A2':0.07):0.02,('BOLD:A3':0.03,'BOLD:A4':0.04):0.03);",
        schema="newick", rooting="force-rooted",
    )
    dated = lsd2_ultrametric(tree, 8.0, seq_length=300)
    assert {l.taxon.label for l in dated.leaf_node_iter()} == {"BOLD:A1", "BOLD:A2", "BOLD:A3", "BOLD:A4"}
    assert is_ultrametric(dated, tolerance=1e-4)
    assert root_age(dated) == pytest.approx(8.0, rel=1e-4)


def test_vectorised_k2p_matches_pairwise():
    rng = np.random.default_rng(4)
    seqs = ["".join(rng.choice(list("ACGT-N"), size=120, p=[0.24, 0.24, 0.24, 0.24, 0.02, 0.02])) for _ in range(12)]
    fast = distances.matrix(seqs, "k2p")
    for i in range(12):
        for j in range(12):
            if i != j:
                assert fast[i, j] == pytest.approx(distances.k2p(seqs[i], seqs[j]))


def test_bionj_on_a_larger_additive_tree():
    tree = dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0.0, num_extant_tips=40,
                                                      rng=__import__("random").Random(2))
    labels = [l.taxon.label for l in tree.leaf_node_iter()]
    pdm = tree.phylogenetic_distance_matrix()
    ns = tree.taxon_namespace
    d = np.array([[pdm(ns.get_taxon(label=a), ns.get_taxon(label=b)) if a != b else 0.0 for b in labels] for a in labels])
    rebuilt = bionj(labels, d)
    pdm2 = rebuilt.phylogenetic_distance_matrix()
    ns2 = rebuilt.taxon_namespace
    for a in labels[:10]:
        for b in labels[10:20]:
            assert pdm2(ns2.get_taxon(label=a), ns2.get_taxon(label=b)) == pytest.approx(
                pdm(ns.get_taxon(label=a), ns.get_taxon(label=b)), rel=1e-6)


def test_timetree_crown_age_is_an_upper_bound():
    tree = dendropy.Tree.get(data="(A:0.05,B:0.05);", schema="newick", rooting="force-rooted")
    # age from the barcodes: 0.05 / 0.01 = 5 Ma
    below = crown_age("G", 50, tree, {}, method="rate", rate=0.01, suggested={"G": 12.0})
    assert (below.source, below.age, below.timetree, below.barcode) == ("rate", 5.0, 12.0, 5.0)
    above = crown_age("G", 50, tree, {}, method="rate", rate=0.01, suggested={"G": 3.0})
    assert (above.source, above.age, above.barcode) == ("rate_timetree_max", 3.0, 5.0)
    # an inflated TimeTree age (a genus that is not a natural group) changes nothing, at any size
    for inflated in (14.9, 15.1, 40.0, 400.0):
        assert crown_age("G", 50, tree, {}, method="rate", rate=0.01, suggested={"G": inflated}).age == 5.0
    ignored = crown_age("G", 50, tree, {}, method="rate", rate=0.01, suggested={"G": 3.0}, timetree_bound="none")
    assert (ignored.source, ignored.age) == ("rate", 5.0)
    fraction = crown_age("G", 30, tree, {}, method="fraction", fraction=0.5, suggested={"G": 10.0})
    assert (fraction.source, fraction.age) == ("fraction_timetree_max", 10.0)
    user_wins = crown_age("G", 50, tree, {"G": 20.0}, method="rate", rate=0.01, suggested={"G": 12.0})
    assert user_wins.source == "table"
    with pytest.raises(ValueError):
        crown_age("G", 50, tree, {}, method="rate", rate=0.01, timetree_bound="min")


def test_crown_age_has_no_threshold_jump():
    """A small change in the barcode age moves the crown age by as little (the Idas case of #36)."""
    ages = []
    for length in (0.0900, 0.1029, 0.1031, 0.1200):
        tree = dendropy.Tree.get(data=f"(A:{length},B:{length});", schema="newick", rooting="force-rooted")
        ages.append(crown_age("Idas", 60, tree, {}, method="rate", rate=0.01, suggested={"Idas": 30.9}).age)
    assert ages == pytest.approx([9.0, 10.29, 10.31, 12.0])


def test_table_bounds_on_the_crown_age():
    tree = dendropy.Tree.get(data="(A:0.05,B:0.05);", schema="newick", rooting="force-rooted")
    kw = dict(method="rate", rate=0.01, suggested={"G": 3.0})  # barcode age 5 Ma
    raised = crown_age("G", 50, tree, {"G": {"min": 8.0}}, **kw)
    assert (raised.source, raised.age) == ("rate_table_min", 8.0)  # the table replaces TimeTree
    lowered = crown_age("G", 50, tree, {"G": {"max": 4.0}}, **kw)
    assert (lowered.source, lowered.age) == ("rate_table_max", 4.0)
    inside = crown_age("G", 50, tree, {"G": {"min": 2.0, "max": 9.0}}, **kw)
    assert (inside.source, inside.age) == ("rate", 5.0)
    # the stem age wins over a minimum above it
    stem_wins = crown_age("G", 6, tree, {"G": {"min": 8.0}}, **kw)
    assert (stem_wins.source, stem_wins.age) == ("rate_capped", pytest.approx(5.4))
    point = crown_age("G", 50, tree, {"G": {"age": 7.0, "min": 8.0}}, **kw)
    assert (point.source, point.age) == ("table", 7.0)
    beyond = crown_age("G", 6, tree, {"G": 61.0}, **kw)
    assert (beyond.source, beyond.age) == ("table_capped", pytest.approx(5.4))
