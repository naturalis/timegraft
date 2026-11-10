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
    "timetree_crown",
    "rooting",
    "violations",
    "misplaced",
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
