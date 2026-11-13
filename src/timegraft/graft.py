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
    "crown_basis",
    "crown_plausibility",
    "timetree_ratio",
    "ratio_allowed",
    "crown_decision",
    "root_to_tip",
    "rate_outliers",
    "rate_outlier_z",
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
    placed: set[str] | frozenset[str] = frozenset(),
) -> dict[str, float]:
    """Move genus crown nodes to other ages, in place, keeping the tree ultrametric.

    For each genus in ``crown_ages`` with at least two of its ``genus_units``
    in the tree, the MRCA of those units is the crown node. All branches below
    it are multiplied by the new age over the present one. Mean path length
    and LSD2 place nodes at fixed fractions of the crown age, so this gives
    the tree that a run with those crown ages would give, without rebuilding
    the subtrees.

    Tips in ``placed`` (units and sequences added by placement) may hang on
    the stem between the crown and the backbone node above it. Their joints
    keep their relative position between the crown and that node, so the
    crown can move up to the backbone node. Anything else on the stem stops
    the crown just below it. Returns the age each genus was moved to. Used
    for diversity at the bounds of the crown ages (``crown_lower``,
    ``crown_upper`` in provenance.tsv).
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

    def age_of(node):
        age = 0.0
        while node.child_nodes():  # the tree is ultrametric: any path down gives the age
            node = node.child_nodes()[0]
            age += node.edge.length or 0.0
        return age

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

    moved = {}
    for genus, (crown, target) in crowns.items():
        age = age_of(crown)
        if age <= 0:
            continue
        # joints on the stem whose other children are all placed tips
        chain, node = [], crown
        while node.parent_node is not None:
            others = [c for c in node.parent_node.child_nodes() if c is not node]
            if others and all(c.is_leaf() and c.taxon is not None and c.taxon.label in placed for c in others):
                chain.append(node.parent_node)
                node = node.parent_node
            else:
                break
        top = chain[-1] if chain else crown
        stem = age_of(top) + (top.edge.length or 0.0)  # age of the node above the stem
        old_ages = [age_of(j) for j in chain]  # before anything moves
        new = max(0.0, min(float(target), stem * (1 - 1e-9)))
        factor = new / age
        for node in crown.preorder_iter():
            if node is not crown and node.edge.length is not None:
                node.edge.length *= factor
        # joints keep their relative position between crown and stem
        span_old, span_new = stem - age, stem - new
        new_ages = [new + (a - age) * (span_new / span_old if span_old > 0 else 0.0) for a in old_ages]
        ages = [new] + new_ages
        lineage = [crown] + chain
        for i, node in enumerate(lineage):
            above = ages[i + 1] if i + 1 < len(ages) else stem
            node.edge.length = above - ages[i]
        for joint, joint_age in zip(chain, new_ages):
            for child in joint.child_nodes():
                if child.is_leaf():
                    child.edge.length = joint_age
        moved[genus] = new
    return moved
