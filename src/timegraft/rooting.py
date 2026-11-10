"""Rooting a genus tree: on outgroup units from nearby genera, or without them.

Two rules root on an outgroup and then drop it:

``strict``
    The ingroup must be a clade with respect to the outgroup; if one unit is out
    of place, the outgroup is dropped and the fallback is used.
``tolerant``
    The root goes on the edge with the fewest units on the wrong side (minimum
    violation): ingroup units among the outgroup, or outgroup units inside the
    ingroup. Ingroup units on the wrong side stay in the tree, as the earliest
    branches of the genus, and are reported.

Without an outgroup, the tree is rooted by MAD (minimal ancestor deviation,
Tria et al. 2017) or at its midpoint.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Iterable

import dendropy
import numpy as np

def _leaf_sets(tree: dendropy.Tree) -> dict[dendropy.Node, frozenset[str]]:
    sets: dict[dendropy.Node, frozenset[str]] = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            sets[node] = frozenset([node.taxon.label])
        else:
            sets[node] = frozenset().union(*(sets[c] for c in node.child_node_iter()))
    return sets


def _detach(tree: dendropy.Tree, node: dendropy.Node) -> dendropy.Tree:
    if node.parent_node is not None:
        node.parent_node.remove_child(node)
    node.edge.length = None
    out = dendropy.Tree(seed_node=node, taxon_namespace=tree.taxon_namespace, is_rooted=True)
    out.suppress_unifurcations()
    return out


def _neighbours(node: dendropy.Node) -> list[tuple[dendropy.Node, float]]:
    out = [(c, c.edge.length or 0.0) for c in node.child_node_iter()]
    if node.parent_node is not None:
        out.append((node.parent_node, node.edge.length or 0.0))
    return out


def _farthest(start: dendropy.Node) -> tuple[dendropy.Node, float, dict]:
    dist = {start: 0.0}
    prev: dict[dendropy.Node, dendropy.Node | None] = {start: None}
    stack = [start]
    while stack:
        node = stack.pop()
        for other, length in _neighbours(node):
            if other not in dist:
                dist[other] = dist[node] + length
                prev[other] = node
                stack.append(other)
    leaves = [n for n in dist if n.is_leaf()]
    end = max(leaves, key=lambda n: dist[n])
    return end, dist[end], prev


def midpoint(tree: dendropy.Tree) -> dendropy.Tree:
    """Root halfway along the longest leaf-to-leaf path."""
    leaves = tree.leaf_nodes()
    if len(leaves) < 2:
        tree.is_rooted = True
        return tree
    a, _, _ = _farthest(leaves[0])
    b, diameter, prev = _farthest(a)
    half = diameter / 2.0
    # walk from b back towards a until the midpoint is passed
    node, travelled = b, 0.0
    while prev[node] is not None:
        step_to = prev[node]
        length = node.edge.length if step_to is node.parent_node else step_to.edge.length
        length = length or 0.0
        if travelled + length >= half:
            edge = node.edge if step_to is node.parent_node else step_to.edge
            child = edge.head_node
            # distance from the edge's child end to the midpoint
            from_node = half - travelled
            to_child = from_node if child is node else length - from_node
            # DendroPy gives length2 to the edge's head (child) side
            tree.reroot_at_edge(edge, length1=length - to_child, length2=to_child, update_bipartitions=False)
            break
        travelled += length
        node = step_to
    tree.is_rooted = True
    tree.suppress_unifurcations()
    return tree


def copy_tree(tree: dendropy.Tree) -> dendropy.Tree:
    return dendropy.Tree.get(data=tree.as_string(schema="newick", suppress_rooting=True),
                             schema="newick", preserve_underscores=True)


def _all_node_distances(tree: dendropy.Tree):
    nodes = list(tree.preorder_node_iter())
    index = {node: i for i, node in enumerate(nodes)}
    adj = defaultdict(list)
    for node in nodes:
        if node.parent_node is not None:
            length = max(0.0, node.edge.length or 0.0)
            adj[index[node]].append((index[node.parent_node], length))
            adj[index[node.parent_node]].append((index[node], length))
    leaves = [index[n] for n in nodes if n.is_leaf()]
    L = np.zeros((len(leaves), len(nodes)))
    for row, start in enumerate(leaves):
        dist = {start: 0.0}
        stack = [start]
        while stack:
            u = stack.pop()
            for v, length in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + length
                    stack.append(v)
        L[row] = [dist[i] for i in range(len(nodes))]
    return nodes, index, leaves, L


def mad_scores(tree: dendropy.Tree) -> list[tuple[float, frozenset[str], float]]:
    """MAD score of the best root on every edge, best first.

    Each entry holds the score, the leaf set on the child side of the edge and
    the distance of the root from that child.
    """
    nodes, index, leaves, L = _all_node_distances(tree)
    n = len(leaves)
    leaf_col = np.array(leaves)
    D = L[:, leaf_col]
    iu = np.triu_indices(n, 1)
    usable = D[iu] > 0
    n_pairs = int(usable.sum())
    sets = _leaf_sets(tree)
    label_row = {nodes[c].taxon.label: r for r, c in enumerate(leaves)}

    def within(rows: np.ndarray, col: int) -> float:
        if len(rows) < 2:
            return 0.0
        lu = L[rows, col]
        sub = D[np.ix_(rows, rows)]
        k = np.triu_indices(len(rows), 1)
        d = sub[k]
        diff = (lu[:, None] - lu[None, :])[k]
        keep = d > 0
        return float(np.sum((diff[keep] / d[keep]) ** 2))

    scores = []
    for node in nodes:
        parent = node.parent_node
        if parent is None:
            continue
        length = max(0.0, node.edge.length or 0.0)
        u, v = index[parent], index[node]
        B = np.array(sorted(label_row[l] for l in sets[node]))
        mask = np.ones(n, bool)
        mask[B] = False
        A = np.nonzero(mask)[0]
        if len(A) == 0:
            continue
        ss = within(A, u) + within(B, v)
        dbc = D[np.ix_(A, B)]
        ok = dbc > 0
        a = (2.0 * L[A, u][:, None] / np.where(ok, dbc, 1.0)) - 1.0
        w = np.where(ok, 1.0 / np.where(ok, dbc, 1.0), 0.0)
        a = np.where(ok, a, 0.0)
        denom = 2.0 * np.sum(w ** 2)
        x = -np.sum(a * w) / denom if denom > 0 else 0.0
        x = min(max(x, 0.0), length)
        ss += float(np.sum(np.where(ok, (a + 2.0 * x * w) ** 2, 0.0)))
        scores.append((math.sqrt(ss / max(n_pairs, 1)), sets[node], length - x))
    scores.sort(key=lambda t: t[0])
    return scores


def mad_root(tree: dendropy.Tree) -> tuple[frozenset[str], float, float, float]:
    """Minimal ancestor deviation root (Tria et al. 2017) of an unrooted tree.

    Returns the leaf set on the child side of the root edge, the distance of the
    root from that child, the MAD score and the root ambiguity index (best score
    over second best; close to 1 means another root is about as good).
    """
    scores = mad_scores(tree)
    best = scores[0]
    rai = best[0] / scores[1][0] if len(scores) > 1 and scores[1][0] > 0 else float("nan")
    return best[1], best[2], best[0], rai


def reroot(tree: dendropy.Tree, child_side: frozenset[str], to_child: float) -> dendropy.Tree:
    labels = frozenset(l.taxon.label for l in tree.leaf_node_iter())
    for node, below in _leaf_sets(tree).items():
        if node.parent_node is None:
            continue
        if below == child_side or below == labels - child_side:
            length = max(0.0, node.edge.length or 0.0)
            d = to_child if below == child_side else length - to_child
            tree.reroot_at_edge(node.edge, length1=length - d, length2=d, update_bipartitions=False)
            tree.is_rooted = True
            tree.suppress_unifurcations()
            return tree
    raise ValueError("root edge not found")


leaf_sets = _leaf_sets


def mad(tree: dendropy.Tree) -> dendropy.Tree:
    """Root ``tree`` by minimal ancestor deviation."""
    if len(tree.leaf_nodes()) < 2:
        tree.is_rooted = True
        return tree
    side, to_child, _, _ = mad_root(tree)
    return reroot(tree, side, to_child)


FALLBACKS = {"mad": mad, "midpoint": midpoint}


def _drop_outgroup(tree: dendropy.Tree, outgroup: frozenset[str]) -> dendropy.Tree:
    if outgroup:
        tree.prune_taxa_with_labels(list(outgroup))
    while tree.seed_node.num_child_nodes() == 1:
        child = tree.seed_node.child_nodes()[0]
        tree.seed_node.remove_child(child)
        child.edge.length = None
        tree.seed_node = child
    tree.is_rooted = True
    tree.suppress_unifurcations()
    return tree


def minimum_violation_sides(
    tree: dendropy.Tree, ingroup: frozenset[str], outgroup: frozenset[str]
) -> tuple[int, list[tuple[frozenset[str], frozenset[str], float]]]:
    """The edges that best separate ingroup from outgroup (minimum violation).

    Returns the number of units on the wrong side, and for every edge that
    attains it: the leaf set below the edge, the side of the edge holding the
    ingroup, and the edge length.
    """
    labels = ingroup | outgroup
    best, ties = None, []
    for node, below in _leaf_sets(tree).items():
        if node.parent_node is None:
            continue
        for side in (below, labels - below):
            if not ingroup & side:
                continue
            cost = len(ingroup - side) + len(outgroup & side)
            if best is None or cost < best:
                best, ties = cost, []
            if cost == best:
                ties.append((below, side, max(0.0, node.edge.length or 0.0)))
    return best, ties


def root_genus(
    tree: dendropy.Tree,
    ingroup: Iterable[str],
    outgroup: Iterable[str],
    rule: str = "tolerant",
    fallback: str = "mad",
) -> tuple[dendropy.Tree, str, dict]:
    """Root an unrooted genus tree and return the ingroup only.

    Returns the rooted tree, how it was rooted and details:
    ``outgroup``: the ingroup is a clade with respect to the outgroup;
    ``outgroup_tolerant``: rooted at the minimum-violation edge, with units on the wrong side;
    ``mad`` or ``midpoint``: no outgroup was given;
    ``mad_outgroup_conflict`` or ``midpoint_outgroup_conflict``: strict rule, ingroup not a clade.
    The details hold ``violations`` (units on the wrong side) and ``misplaced``
    (ingroup units among the outgroup, kept in the tree).
    """
    if rule not in ("strict", "tolerant"):
        raise ValueError(f"unknown rooting rule {rule!r}")
    if fallback not in FALLBACKS:
        raise ValueError(f"unknown rooting fallback {fallback!r}")
    ingroup = frozenset(ingroup)
    outgroup = frozenset(outgroup) - ingroup
    labels = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    if not ingroup <= labels:
        raise ValueError("ingroup labels missing from tree")
    outgroup &= labels
    tree.prune_taxa_with_labels(list(labels - ingroup - outgroup))

    if not outgroup:
        return FALLBACKS[fallback](tree), fallback, {"violations": "", "misplaced": []}

    violations, ties = minimum_violation_sides(tree, ingroup, outgroup)
    if violations == 0:
        anchor = next(leaf for leaf in tree.leaf_node_iter() if leaf.taxon.label in outgroup)
        alen = anchor.edge.length or 0.0
        tree.reroot_at_edge(anchor.edge, length1=alen / 2.0, length2=alen / 2.0, update_bipartitions=False)
        for node, leaves in _leaf_sets(tree).items():
            if leaves == ingroup:
                return _detach(tree, node), "outgroup", {"violations": 0, "misplaced": []}
        raise AssertionError("ingroup split not found")  # pragma: no cover

    if rule == "strict":
        _drop_outgroup(tree, outgroup)
        return FALLBACKS[fallback](tree), f"{fallback}_outgroup_conflict", {"violations": violations, "misplaced": []}

    # each tied edge roots the genus somewhere; where they differ, MAD on the genus tree decides
    candidates = {}
    for below, side, length in ties:
        rooted = _drop_outgroup(reroot(copy_tree(tree), below, length / 2.0), outgroup)
        split = _root_split(rooted)
        if split not in candidates:
            candidates[split] = (rooted, sorted(ingroup - side))
    if len(candidates) > 1:
        genus_tree = _drop_outgroup(copy_tree(tree), outgroup)
        score = {}
        for value, child, _ in mad_scores(genus_tree):
            labels_all = frozenset(l.taxon.label for l in genus_tree.leaf_node_iter())
            key = frozenset([child, labels_all - child])
            score.setdefault(key, value)
        best = min(candidates, key=lambda k: (score.get(k, math.inf), sorted(map(sorted, k))))
    else:
        best = next(iter(candidates))
    rooted, misplaced = candidates[best]
    return rooted, "outgroup_tolerant", {"violations": violations, "misplaced": misplaced}


def _root_split(tree: dendropy.Tree) -> frozenset[frozenset[str]]:
    """The root bipartition; a basal polytomy is split off its first child."""
    kids = tree.seed_node.child_nodes()
    first = frozenset(l.taxon.label for l in kids[0].leaf_iter())
    rest = frozenset(l.taxon.label for c in kids[1:] for l in c.leaf_iter())
    return frozenset([first, rest])


def root_on_outgroup(
    tree: dendropy.Tree, ingroup: Iterable[str], outgroup: Iterable[str]
) -> tuple[dendropy.Tree, str]:
    """The strict rule with midpoint fallback, as in timegraft before #19.

    Returns the rooted ingroup tree and ``outgroup``, ``midpoint`` (no outgroup)
    or ``midpoint_outgroup_conflict``.
    """
    rooted, how, _ = root_genus(tree, ingroup, outgroup, rule="strict", fallback="midpoint")
    return rooted, how
