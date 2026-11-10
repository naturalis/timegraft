"""The dated backbone: an ultrametric tree in millions of years with genera at the tips."""

from __future__ import annotations

import zlib
from pathlib import Path

import dendropy

from .dating import is_ultrametric


def read_backbone(path: str | Path) -> dendropy.Tree:
    tree = dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True, rooting="force-rooted")
    labels = [leaf.taxon.label for leaf in tree.leaf_node_iter()]
    if len(labels) != len(set(labels)):
        raise ValueError(f"{path}: duplicate genus labels in backbone")
    if any(leaf.edge.length is None for leaf in tree.leaf_node_iter()):
        raise ValueError(f"{path}: backbone needs branch lengths (ages in Ma)")
    if not is_ultrametric(tree, tolerance=1e-4):
        raise ValueError(f"{path}: backbone is not ultrametric")
    return tree


def stem_ages(tree: dendropy.Tree) -> dict[str, float]:
    """Stem age of each genus: the length of its tip edge in an ultrametric tree."""
    return {leaf.taxon.label: float(leaf.edge.length) for leaf in tree.leaf_node_iter()}


def nearest_genera(tree: dendropy.Tree, genus: str, available: set[str]) -> list[str]:
    """Genera in the closest clade beside ``genus`` that have units, sorted.

    Walks up from the genus tip and returns the available genera below the
    first ancestor that has any, i.e. the sister group in the backbone.
    """
    leaf = tree.find_node_with_taxon_label(genus)
    if leaf is None:
        raise KeyError(genus)
    node = leaf
    while node.parent_node is not None:
        node = node.parent_node
        found = sorted(
            l.taxon.label
            for l in node.leaf_iter()
            if l.taxon.label != genus and l.taxon.label in available
        )
        if found:
            return found
    return []


def shard_of(genus: str, n_shards: int) -> int:
    return zlib.crc32(genus.encode()) % n_shards


def sister_clades(tree: dendropy.Tree, genus: str, available: set[str], levels: int = 3) -> list[list[str]]:
    """Genera with units in successively more distant clades beside ``genus``."""
    node = tree.find_node_with_taxon_label(genus)
    seen = {genus}
    clades = []
    while node.parent_node is not None and len(clades) < levels:
        node = node.parent_node
        found = sorted(l.taxon.label for l in node.leaf_iter() if l.taxon.label in available and l.taxon.label not in seen)
        if found:
            clades.append(found)
            seen.update(found)
    return clades


def _sites(seq: str) -> int:
    return sum(ch in "ACGTacgt" for ch in seq)


def pick_spread(clades, by_genus, n: int, n_genera: int | None = None) -> list[tuple[str, str]]:
    """Round-robin over sister genera, nearest clade first, longest sequences first."""
    genera = [g for clade in clades for g in clade]
    if n_genera is not None:
        genera = genera[:n_genera]
    queues = {g: sorted(by_genus[g], key=lambda r: (-_sites(r["sequence"]), r["unit_id"])) for g in genera}
    picked = []
    while len(picked) < n and any(queues.values()):
        for g in genera:
            if queues[g] and len(picked) < n:
                r = queues[g].pop(0)
                picked.append((r["unit_id"], r["sequence"]))
    return picked


def pick_outgroup(
    tree: dendropy.Tree,
    genus: str,
    by_genus: dict[str, list[dict]],
    available: set[str],
    n: int,
    choice: str = "spread",
    levels: int = 3,
) -> list[dict]:
    """Outgroup units for ``genus``.

    ``spread``: up to ``n`` units taken in turn from the genera of the nearest
    ``levels`` sister clades, nearest clade first, longest sequences first.
    ``sister``: the first ``n`` units, by unit id, of the genera in the nearest
    sister clade, in alphabetical order (timegraft before #19).
    """
    if choice == "sister":
        picked = []
        for sister in nearest_genera(tree, genus, available):
            for unit in sorted(by_genus[sister], key=lambda r: r["unit_id"]):
                if len(picked) < n:
                    picked.append(unit)
        return picked
    if choice != "spread":
        raise ValueError(f"unknown outgroup choice {choice!r}")
    genera = [g for clade in sister_clades(tree, genus, available, levels) for g in clade]
    queues = {g: sorted(by_genus[g], key=lambda r: (-_sites(r["sequence"]), r["unit_id"])) for g in genera}
    picked = []
    while len(picked) < n and any(queues.values()):
        for g in genera:
            if queues[g] and len(picked) < n:
                picked.append(queues[g].pop(0))
    return picked
