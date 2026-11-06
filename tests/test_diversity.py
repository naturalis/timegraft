"""Diversity metrics checked against brute-force calculations."""

import itertools
import random

import dendropy
import numpy as np
import pytest

from timegraft.diversity import alpha, beta, edge_table

TREE = "(((a:1,b:1):2,c:3):1,(d:2,e:2):2);"


def load():
    return dendropy.Tree.get(data=TREE, schema="newick", rooting="force-rooted")


def brute_patristic(tree):
    pdm = tree.phylogenetic_distance_matrix()
    ns = tree.taxon_namespace
    return lambda x, y: pdm(ns.get_taxon(label=x), ns.get_taxon(label=y))


def brute_pd(tree, units):
    edges = set()
    for leaf in tree.leaf_node_iter():
        if leaf.taxon.label in units:
            node = leaf
            while node.parent_node is not None:
                edges.add(node)
                node = node.parent_node
    return sum(n.edge.length for n in edges)


def test_faith_pd_and_mpd_by_hand():
    table, missing = edge_table(load(), {"s1": {"a": 1, "b": 1}, "s2": {"a": 1, "d": 1, "x": 5}})
    rows = {r["sample"]: r for r in alpha(table)}
    # s1: edges a, b, (a,b), ((a,b),c) = 1 + 1 + 2 + 1
    assert rows["s1"]["faith_pd"] == pytest.approx(5.0)
    assert rows["s1"]["mpd"] == pytest.approx(2.0)
    # s2: a path 1+2+1, d path 2+2 = 8; distance a-d = 8
    assert rows["s2"]["faith_pd"] == pytest.approx(8.0)
    assert rows["s2"]["mpd"] == pytest.approx(8.0)
    assert missing == {"s1": [], "s2": ["x"]}


def test_unifrac_and_phylosor_extremes():
    table, _ = edge_table(load(), {"s1": {"a": 2, "b": 1}, "s2": {"a": 2, "b": 1}, "s3": {"d": 1, "e": 1}})
    rows = {(r["sample_a"], r["sample_b"]): r for r in beta(table)}
    same = rows[("s1", "s2")]
    assert same["unifrac_unweighted"] == pytest.approx(0.0)
    assert same["unifrac_weighted"] == pytest.approx(0.0)
    assert same["phylosor"] == pytest.approx(1.0)
    apart = rows[("s1", "s3")]
    assert apart["unifrac_unweighted"] == pytest.approx(1.0)
    assert apart["unifrac_weighted"] == pytest.approx(1.0)
    assert apart["phylosor"] == pytest.approx(0.0)


def test_against_brute_force_on_random_samples():
    rng = random.Random(3)
    tree = dendropy.simulate.treesim.birth_death_tree(
        birth_rate=1.0, death_rate=0.2, num_extant_tips=30, rng=rng
    )
    labels = [leaf.taxon.label for leaf in tree.leaf_node_iter()]
    samples = {}
    for s in range(5):
        chosen = rng.sample(labels, rng.randint(2, 12))
        samples[f"s{s}"] = {u: rng.randint(1, 50) for u in chosen}
    table, _ = edge_table(tree, samples)
    dist = brute_patristic(tree)

    for row in alpha(table):
        units = samples[row["sample"]]
        total = sum(units.values())
        p = {u: a / total for u, a in units.items()}
        assert row["faith_pd"] == pytest.approx(brute_pd(tree, units))
        pairs = list(itertools.combinations(units, 2))
        assert row["mpd"] == pytest.approx(np.mean([dist(x, y) for x, y in pairs]))
        num = sum(p[x] * p[y] * dist(x, y) for x, y in itertools.permutations(units, 2))
        den = sum(p[x] * p[y] for x, y in itertools.permutations(units, 2))
        assert row["mpd_abundance"] == pytest.approx(num / den)

    for row in beta(table):
        a, b = samples[row["sample_a"]], samples[row["sample_b"]]
        pd_a, pd_b = brute_pd(tree, a), brute_pd(tree, b)
        pd_union = brute_pd(tree, set(a) | set(b))
        shared = pd_a + pd_b - pd_union
        assert row["unifrac_unweighted"] == pytest.approx((pd_union - shared) / pd_union)
        assert row["phylosor"] == pytest.approx(shared / (0.5 * (pd_a + pd_b)))
