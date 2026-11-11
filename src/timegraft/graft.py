"""Grafting dated genus subtrees onto the dated backbone."""

from __future__ import annotations

import dendropy

PROVENANCE_COLUMNS = (
    "genus",
    "status",
    "n_units",
    "stem_age",
    "crown_age",
    "crown_source",
    "crown_lower",
    "crown_upper",
    "crown_upper_source",
    "barcode_age",
    "timetree_crown",
    "root_to_tip",
    "rate",
    "rate_source",
    "rooting",
    "violations",
    "misplaced",
    "root_ambiguity",
    "root_alternatives",
    "crown_age_min",
    "crown_age_max",
    "pd_rel_range",
    "dating",
    "outgroup",
)


def graft(
    backbone: dendropy.Tree,
    subtrees: dict[str, tuple[str, float]],
    singles: dict[str, str],
    prune_empty: bool = True,
) -> tuple[dendropy.Tree, dict[str, str]]:
    """Replace genus tips with their subtrees, in place.

    ``subtrees`` maps genus to (newick, crown age); ``singles`` maps genus to
    the id of its only unit. For a genus with stem age S and crown age C the
    tip edge is replaced by a stem of length S - C followed by the subtree,
    which keeps the tree ultrametric. Genera with no units are pruned unless
    ``prune_empty`` is False. Returns the tree and a status per backbone genus.
    """
    namespace = backbone.taxon_namespace
    status: dict[str, str] = {}
    empty = []
    for tip in list(backbone.leaf_node_iter()):
        genus = tip.taxon.label
        if genus in subtrees:
            newick, crown = subtrees[genus]
            stem = float(tip.edge.length)
            if not crown < stem:
                raise ValueError(f"{genus}: crown age {crown} is not below stem age {stem}")
            sub = dendropy.Tree.get(
                data=newick, schema="newick", taxon_namespace=namespace, preserve_underscores=True, rooting="force-rooted"
            )
            root = sub.seed_node
            parent = tip.parent_node
            parent.remove_child(tip)
            parent.add_child(root)
            root.edge.length = stem - crown
            status[genus] = "subtree"
        elif genus in singles:
            tip.taxon = namespace.require_taxon(label=singles[genus])
            status[genus] = "single_unit"
        else:
            empty.append(tip.taxon)
            status[genus] = "pruned" if prune_empty else "no_units"
    if prune_empty and empty:
        backbone.prune_taxa(empty, suppress_unifurcations=True)
    for taxon in empty if prune_empty else []:
        namespace.remove_taxon(taxon)
    return backbone, status


def rescale_crowns(
    tree: dendropy.Tree,
    genus_units: dict[str, list[str]],
    crown_ages: dict[str, float],
) -> dict[str, float]:
    """Move genus crown nodes to other ages, in place, keeping the tree ultrametric.

    For each genus in ``crown_ages`` with at least two of its ``genus_units``
    in the tree, the MRCA of those units is the crown node. All branches below
    it are multiplied by the new age over the present one, and the branch
    above it absorbs the difference. Mean path length and LSD2 place nodes at
    fixed fractions of the crown age, so this gives the tree that a run with
    those crown ages would give, without rebuilding the subtrees.

    A new age is kept just below the age of the node above the crown (which
    can be a placed unit on the stem, not the stem node itself). Returns the
    age each genus was moved to. Used for diversity at the bounds of the crown
    ages (``crown_lower``, ``crown_upper`` in provenance.tsv).
    """
    tips = {leaf.taxon.label: leaf for leaf in tree.leaf_node_iter()}
    depth = {}
    for node in tree.preorder_node_iter():
        depth[node] = 0 if node.parent_node is None else depth[node.parent_node] + 1

    def mrca(a, b):
        while a is not b:
            if depth[a] >= depth[b]:
                a = a.parent_node
            else:
                b = b.parent_node
        return a

    moved = {}
    crowns = {}
    for genus, target in crown_ages.items():
        nodes = [tips[u] for u in genus_units.get(genus, ()) if u in tips]
        if len(nodes) < 2:
            continue
        crown = nodes[0]
        for node in nodes[1:]:
            crown = mrca(crown, node)
        if crown.parent_node is not None:
            crowns[genus] = (crown, target)
    for genus, (crown, target) in crowns.items():
        leaf, age = crown, 0.0
        while leaf.child_nodes():  # the tree is ultrametric: any path down gives the age
            leaf = leaf.child_nodes()[0]
            age += leaf.edge.length or 0.0
        if age <= 0:
            continue
        above = age + (crown.edge.length or 0.0)
        new = max(0.0, min(float(target), above * (1 - 1e-9)))
        factor = new / age
        for node in crown.preorder_iter():
            if node is not crown and node.edge.length is not None:
                node.edge.length *= factor
        crown.edge.length = above - new
        moved[genus] = new
    return moved

