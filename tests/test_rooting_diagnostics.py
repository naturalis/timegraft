import itertools
import math
import random

import dendropy
import numpy as np
import pytest

from timegraft.rooting_diagnostics import (
    copy_tree,
    k2p_from_pq,
    leaf_sets,
    mad_root,
    pq_matrices,
    encode,
    reroot,
    rogue_counts,
    root_split,
)
from timegraft import distances


def random_unrooted(n, seed):
    rng = random.Random(seed)
    tree = dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0.0, num_extant_tips=n,
                                                      rng=rng)
    for i, leaf in enumerate(tree.leaf_node_iter()):
        leaf.taxon.label = f"t{i}"
    for edge in tree.postorder_edge_iter():
        if edge.tail_node is not None:
            edge.length = rng.uniform(0.01, 1.0)  # strongly non-clocklike
    tree.is_rooted = False
    tree.collapse_basal_bifurcation()
    return copy_tree(tree)


def patristic(tree):
    pdm = tree.phylogenetic_distance_matrix()
    leaves = list(tree.taxon_namespace)
    return {(a.label, b.label): pdm(a, b) for a in leaves for b in leaves}


def mad_brute(tree, steps=400):
    """MAD by explicit rerooting on a fine grid along every edge."""
    labels = [l.taxon.label for l in tree.leaf_node_iter()]
    best = (math.inf, None)
    for node, below in leaf_sets(tree).items():
        if node.parent_node is None:
            continue
        length = node.edge.length
        for k in range(steps + 1):
            to_child = length * k / steps
            t = reroot(copy_tree(tree), below, to_child)
            dist = patristic(t)
            depth = {}
            for nd in t.preorder_node_iter():
                depth[nd] = 0.0 if nd.parent_node is None else depth[nd.parent_node] + nd.edge.length
            leaf = {l.taxon.label: l for l in t.leaf_node_iter()}
            ss = n = 0
            for a, b in itertools.combinations(labels, 2):
                mrca = t.mrca(taxon_labels=[a, b])
                dab = depth[leaf[a]] - depth[mrca]
                ss += (2 * dab / dist[(a, b)] - 1) ** 2
                n += 1
            score = math.sqrt(ss / n)
            if score < best[0] - 1e-12:
                best = (score, root_split(t))
    return best


@pytest.mark.parametrize("seed", range(5))
def test_mad_matches_brute_force(seed):
    tree = random_unrooted(7, seed)
    side, to_child, score, rai = mad_root(copy_tree(tree))
    brute_score, brute_split = mad_brute(tree)
    assert score == pytest.approx(brute_score, rel=1e-3)
    assert root_split(reroot(copy_tree(tree), side, to_child)) == brute_split
    assert 0 < rai <= 1


def rogue_brute(tree, ingroup, outgroup):
    labels = ingroup | outgroup

    def is_split(tree, keep_in, keep_out):
        t = copy_tree(tree)
        t.prune_taxa_with_labels(list(labels - keep_in - keep_out))
        sets = leaf_sets(t)
        allk = keep_in | keep_out
        return any(s == keep_in or allk - s == keep_in for s in sets.values())

    rin = next(k for k in range(len(ingroup)) for drop in itertools.combinations(sorted(ingroup), k)
               if len(ingroup) - k >= 1 and is_split(tree, ingroup - set(drop), outgroup))
    return rin


@pytest.mark.parametrize("seed", range(8))
def test_rogue_in_matches_brute_force(seed):
    rng = random.Random(seed)
    tree = random_unrooted(9, seed)
    labels = sorted(l.taxon.label for l in tree.leaf_node_iter())
    outgroup = frozenset(rng.sample(labels, 2))
    ingroup = frozenset(labels) - outgroup
    assert rogue_counts(tree, ingroup, outgroup)["rogue_in"] == rogue_brute(tree, ingroup, outgroup)


def test_k2p_matrix_matches_package():
    rng = random.Random(1)
    base = "".join(rng.choice("ACGT") for _ in range(300))
    seqs = []
    for i in range(6):
        s = list(base)
        for j in rng.sample(range(300), 20 * i):
            s[j] = rng.choice("ACGT-N")
        seqs.append("".join(s))
    P, Q, S = pq_matrices(encode(seqs))
    ours = k2p_from_pq(P, Q)
    theirs = distances.matrix(seqs, method="k2p")
    assert np.allclose(ours, theirs)
    gamma = k2p_from_pq(P, Q, alpha=0.5)
    assert np.all(gamma >= ours - 1e-12)
