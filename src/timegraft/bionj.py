"""BioNJ tree building (Gascuel 1997) on a distance matrix.

Returns an unrooted DendroPy tree. Negative branch lengths are set to zero.
"""

from __future__ import annotations

from typing import Sequence

import dendropy
import numpy as np


def _leaf(label: str, namespace: dendropy.TaxonNamespace) -> dendropy.Node:
    return dendropy.Node(taxon=namespace.require_taxon(label=label))


def _attach(parent: dendropy.Node, child: dendropy.Node, length: float) -> None:
    parent.add_child(child)
    child.edge.length = max(0.0, float(length))


def bionj(labels: Sequence[str], dist: np.ndarray) -> dendropy.Tree:
    labels = list(labels)
    n = len(labels)
    if dist.shape != (n, n):
        raise ValueError("distance matrix does not match labels")
    if len(set(labels)) != n:
        raise ValueError("labels must be unique")
    namespace = dendropy.TaxonNamespace()
    nodes = [_leaf(label, namespace) for label in labels]

    if n == 1:
        return dendropy.Tree(seed_node=nodes[0], taxon_namespace=namespace, is_rooted=False)
    if n == 2:
        root = dendropy.Node()
        half = dist[0, 1] / 2.0
        _attach(root, nodes[0], half)
        _attach(root, nodes[1], half)
        return dendropy.Tree(seed_node=root, taxon_namespace=namespace, is_rooted=False)

    # room for the n leaves and the n - 3 internal nodes created by joins
    size = 2 * n
    d = np.zeros((size, size))
    v = np.zeros((size, size))
    d[:n, :n] = dist
    v[:n, :n] = dist  # variance estimates, initialised to the distances
    active = list(range(n))

    while len(active) > 3:
        m = len(active)
        idx = np.array(active)
        sub = d[np.ix_(idx, idx)]
        r = sub.sum(axis=1)
        q = (m - 2) * sub - r[:, None] - r[None, :]
        np.fill_diagonal(q, np.inf)
        a, b = np.unravel_index(np.argmin(q), q.shape)
        i, j = active[a], active[b]

        bi = 0.5 * d[i, j] + (r[a] - r[b]) / (2.0 * (m - 2))
        bj = d[i, j] - bi

        others = np.array([k for k in active if k not in (i, j)])
        if v[i, j] > 0:
            lam = 0.5 + (v[j, others] - v[i, others]).sum() / (2.0 * (m - 2) * v[i, j])
            lam = min(1.0, max(0.0, lam))
        else:
            lam = 0.5

        new = dendropy.Node()
        _attach(new, nodes[i], bi)
        _attach(new, nodes[j], bj)
        u = len(nodes)
        nodes.append(new)
        d[u, others] = d[others, u] = lam * d[i, others] + (1 - lam) * d[j, others] - lam * bi - (1 - lam) * bj
        v[u, others] = v[others, u] = lam * v[i, others] + (1 - lam) * v[j, others] - lam * (1 - lam) * v[i, j]
        active = list(others) + [u]

    x, y, z = active
    centre = dendropy.Node()
    _attach(centre, nodes[x], (d[x, y] + d[x, z] - d[y, z]) / 2.0)
    _attach(centre, nodes[y], (d[x, y] + d[y, z] - d[x, z]) / 2.0)
    _attach(centre, nodes[z], (d[x, z] + d[y, z] - d[x, y]) / 2.0)
    return dendropy.Tree(seed_node=centre, taxon_namespace=namespace, is_rooted=False)
