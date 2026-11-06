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

    d = np.array(dist, dtype=float)
    v = d.copy()  # variance estimates, initialised to the distances
    active = list(range(n))

    while len(active) > 3:
        m = len(active)
        sub = d[np.ix_(active, active)]
        r = sub.sum(axis=1)
        q = (m - 2) * sub - r[:, None] - r[None, :]
        np.fill_diagonal(q, np.inf)
        a, b = np.unravel_index(np.argmin(q), q.shape)
        i, j = active[a], active[b]

        bi = 0.5 * d[i, j] + (r[a] - r[b]) / (2.0 * (m - 2))
        bj = d[i, j] - bi

        others = [k for k in active if k not in (i, j)]
        if v[i, j] > 0:
            lam = 0.5 + sum(v[j, k] - v[i, k] for k in others) / (2.0 * (m - 2) * v[i, j])
            lam = min(1.0, max(0.0, lam))
        else:
            lam = 0.5

        new = dendropy.Node()
        _attach(new, nodes[i], bi)
        _attach(new, nodes[j], bj)

        u = len(nodes)
        nodes.append(new)
        size = u + 1
        d2 = np.zeros((size, size))
        v2 = np.zeros((size, size))
        d2[:u, :u] = d
        v2[:u, :u] = v
        for k in others:
            d2[u, k] = d2[k, u] = lam * d[i, k] + (1 - lam) * d[j, k] - lam * bi - (1 - lam) * bj
            v2[u, k] = v2[k, u] = lam * v[i, k] + (1 - lam) * v[j, k] - lam * (1 - lam) * v[i, j]
        d, v = d2, v2
        active = others + [u]

    x, y, z = active
    centre = dendropy.Node()
    _attach(centre, nodes[x], (d[x, y] + d[x, z] - d[y, z]) / 2.0)
    _attach(centre, nodes[y], (d[x, y] + d[y, z] - d[x, z]) / 2.0)
    _attach(centre, nodes[z], (d[x, z] + d[y, z] - d[x, y]) / 2.0)
    return dendropy.Tree(seed_node=centre, taxon_namespace=namespace, is_rooted=False)
