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
