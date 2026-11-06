"""Rooting a genus tree on outgroup units from the sister genus, then dropping them."""

from __future__ import annotations

from typing import Iterable

import dendropy


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


def root_on_outgroup(
    tree: dendropy.Tree, ingroup: Iterable[str], outgroup: Iterable[str]
) -> tuple[dendropy.Tree, str]:
    """Root ``tree`` so that ``ingroup`` is a clade, and return only that clade.

    Returns the rooted ingroup tree and the method used:
    ``outgroup`` when the outgroup separates the ingroup cleanly,
    ``midpoint`` when no outgroup was supplied, and
    ``midpoint_outgroup_conflict`` when the ingroup is not a clade with respect
    to the outgroup (the outgroup is then removed and the ingroup midpoint-rooted).
    """
    ingroup = frozenset(ingroup)
    outgroup = frozenset(outgroup) - ingroup
    labels = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    if not ingroup <= labels:
        raise ValueError("ingroup labels missing from tree")
    outgroup &= labels

    if not outgroup:
        tree.prune_taxa_with_labels(labels - ingroup)
        return midpoint(tree), "midpoint"

    anchor = next(leaf for leaf in tree.leaf_node_iter() if leaf.taxon.label in outgroup)
    length = anchor.edge.length or 0.0
    tree.reroot_at_edge(anchor.edge, length1=length / 2.0, length2=length / 2.0, update_bipartitions=False)

    for node, leaves in _leaf_sets(tree).items():
        if leaves == ingroup:
            return _detach(tree, node), "outgroup"

    tree.prune_taxa_with_labels(outgroup)
    return midpoint(tree), "midpoint_outgroup_conflict"
