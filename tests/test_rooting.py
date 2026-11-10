import dendropy
import pytest

from timegraft.backbone import pick_outgroup
from timegraft.rooting import mad, root_genus, root_on_outgroup


def unrooted(newick):
    return dendropy.Tree.get(data=newick, schema="newick", rooting="force-unrooted")


def leaves(tree):
    return {l.taxon.label for l in tree.leaf_node_iter()}


def root_parts(tree):
    return {frozenset(l.taxon.label for l in c.leaf_iter()) for c in tree.seed_node.child_node_iter()}


# ingroup A-E, outgroup O1, O2
CLEAN = "((A:1,B:1):1,((C:1,D:1):1,E:2):1,(O1:1,O2:1):3);"
# ingroup unit E sits inside the outgroup (merely sister to it would still be a clean split)
INGROUP_OUT = "((A:1,B:1):1,(C:1,D:1):1,((O1:1,(O2:1,E:1):1):1,O3:2):3);"
# outgroup unit O2 sits inside the ingroup
OUTGROUP_IN = "((A:1,B:1):1,((C:1,O2:1):1,(D:1,E:1):1):1,O1:4);"
IN = ["A", "B", "C", "D", "E"]
OUT = ["O1", "O2"]


@pytest.mark.parametrize("rule", ["strict", "tolerant"])
def test_clean_split_is_rooted_on_the_outgroup(rule):
    tree, how, info = root_genus(unrooted(CLEAN), IN, OUT, rule=rule)
    assert how == "outgroup"
    assert info == {"violations": 0, "misplaced": []}
    assert leaves(tree) == set(IN)
    assert root_parts(tree) == {frozenset("AB"), frozenset("CDE")}


def test_tolerant_keeps_a_misplaced_ingroup_unit_as_earliest_branch():
    tree, how, info = root_genus(unrooted(INGROUP_OUT), IN, OUT + ["O3"], rule="tolerant")
    assert how == "outgroup_tolerant"
    assert info == {"violations": 1, "misplaced": ["E"]}
    assert leaves(tree) == set(IN)
    assert root_parts(tree) == {frozenset("E"), frozenset("ABCD")}
    assert tree.is_rooted


def test_tied_roots_are_decided_by_mad():
    # O2 nests with C, so either O1 or O2 is the misfit: rooting on O1 puts the
    # root between AB and CDE, rooting on O2 puts it next to C. On the genus tree
    # (AB:1+1, C, DE), MAD prefers the AB | CDE root.
    tree, how, info = root_genus(unrooted(OUTGROUP_IN), IN, OUT, rule="tolerant")
    assert how == "outgroup_tolerant"
    assert info == {"violations": 1, "misplaced": []}
    assert leaves(tree) == set(IN)
    assert root_parts(tree) == {frozenset("AB"), frozenset("CDE")}


@pytest.mark.parametrize("fallback", ["mad", "midpoint"])
def test_strict_falls_back_when_the_genus_is_not_a_clade(fallback):
    tree, how, info = root_genus(unrooted(INGROUP_OUT), IN, OUT + ["O3"], rule="strict", fallback=fallback)
    assert how == f"{fallback}_outgroup_conflict"
    assert info["violations"] == 1
    assert leaves(tree) == set(IN)


def test_without_outgroup_the_fallback_is_used():
    tree, how, info = root_genus(unrooted(CLEAN), IN, [], rule="tolerant", fallback="mad")
    assert how == "mad"
    assert leaves(tree) == set(IN)
    assert info["violations"] == ""


def test_root_on_outgroup_keeps_the_old_behaviour():
    tree, how = root_on_outgroup(unrooted(INGROUP_OUT), IN, OUT + ["O3"])
    assert how == "midpoint_outgroup_conflict"
    tree, how = root_on_outgroup(unrooted(CLEAN), IN, OUT)
    assert how == "outgroup"


def test_mad_on_two_and_three_leaves():
    assert leaves(mad(unrooted("(A:1,B:2);"))) == {"A", "B"}
    # a clock-like triple: MAD roots between C and the cherry
    tree = mad(unrooted("((A:1,B:1):1,C:2);"))
    assert root_parts(tree) == {frozenset("AB"), frozenset("C")}


def test_unknown_rule_or_fallback():
    with pytest.raises(ValueError):
        root_genus(unrooted(CLEAN), IN, OUT, rule="lenient")
    with pytest.raises(ValueError):
        root_genus(unrooted(CLEAN), IN, OUT, fallback="random")


# outgroup choice --------------------------------------------------------

BACKBONE = "((((G:1,S1:1):1,S2:2):1,(T1:1,T2:1):2):1,U:4);"


def units(genus, lengths):
    return [{"unit_id": f"{genus}_{i}", "genus": genus, "sequence": "A" * n} for i, n in enumerate(lengths)]


def test_pick_outgroup_spreads_over_sister_clades_longest_first():
    tree = dendropy.Tree.get(data=BACKBONE, schema="newick", rooting="force-rooted")
    by_genus = {"G": units("G", [5, 5]), "S1": units("S1", [3, 9, 6]), "S2": units("S2", [4]),
                "T1": units("T1", [5]), "T2": units("T2", [7, 2]), "U": units("U", [9])}
    available = set(by_genus)
    picked = [u["unit_id"] for u in pick_outgroup(tree, "G", by_genus, available, 5, "spread", levels=3)]
    # clades: [S1], [S2], [T1, T2]; round robin, longest sequence first
    assert picked == ["S1_1", "S2_0", "T1_0", "T2_0", "S1_2"]
    two_levels = [u["unit_id"] for u in pick_outgroup(tree, "G", by_genus, available, 10, "spread", levels=2)]
    assert set(two_levels) == {"S1_0", "S1_1", "S1_2", "S2_0"}


def test_pick_outgroup_sister_is_the_old_choice():
    tree = dendropy.Tree.get(data=BACKBONE, schema="newick", rooting="force-rooted")
    by_genus = {"G": units("G", [5, 5]), "S1": units("S1", [3, 9, 6]), "S2": units("S2", [4])}
    picked = [u["unit_id"] for u in pick_outgroup(tree, "G", by_genus, set(by_genus), 2, "sister")]
    assert picked == ["S1_0", "S1_1"]
