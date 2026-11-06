"""Phylogenetic diversity per sample and between samples.

All metrics are computed from one edge-by-sample summary of the tree, which
avoids pruning a separate tree for every sample:

- Faith's PD: total length of the edges that lead to at least one unit in the
  sample, up to and including the path to the root of the tree.
- MPD: mean patristic distance between pairs of units in the sample. Every
  edge contributes its length times the number of pairs it separates, so
  MPD = sum(l * n * (k - n)) / (k * (k - 1) / 2) for presence data, and the
  analogous abundance-weighted form with fractions.
- Unweighted UniFrac: branch length unique to one sample over branch length
  in either.
- Weighted UniFrac, normalised (Lozupone et al. 2007).
- PhyloSor (Bryant et al. 2008): shared branch length over the mean branch
  length of the two samples; a similarity.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import dendropy
import numpy as np


@dataclass
class EdgeTable:
    lengths: np.ndarray          # edge lengths, one per edge used by any sample
    present: np.ndarray          # edges x samples: number of the sample's units below the edge
    weight: np.ndarray           # edges x samples: share of the sample's abundance below the edge
    depth: dict[str, float]      # root-to-tip path length per unit, for weighted UniFrac
    samples: list[str]
    richness: np.ndarray         # units per sample that are in the tree
    tip_depth_weighted: np.ndarray  # per sample: sum over units of abundance share times tip depth
    sum_sq: np.ndarray           # per sample: sum of squared abundance shares


def edge_table(tree: dendropy.Tree, samples: dict[str, dict[str, float]]) -> tuple[EdgeTable, dict[str, list[str]]]:
    """Summarise ``samples`` (sample -> unit -> abundance) over the edges of ``tree``.

    Returns the table and, per sample, the units that are not in the tree.
    """
    leaves = {leaf.taxon.label: leaf for leaf in tree.leaf_node_iter()}
    names = sorted(samples)
    index: dict[dendropy.Node, int] = {}
    lengths: list[float] = []
    present = defaultdict(lambda: np.zeros(len(names)))
    weight = defaultdict(lambda: np.zeros(len(names)))
    depth: dict[str, float] = {}
    missing: dict[str, list[str]] = {}
    richness = np.zeros(len(names))
    tip_depth_weighted = np.zeros(len(names))
    sum_sq = np.zeros(len(names))

    for s, name in enumerate(names):
        units = {u: a for u, a in samples[name].items() if a > 0}
        missing[name] = sorted(u for u in units if u not in leaves)
        units = {u: a for u, a in units.items() if u in leaves}
        total = sum(units.values())
        richness[s] = len(units)
        for unit, abundance in units.items():
            share = abundance / total
            sum_sq[s] += share * share
            node, path = leaves[unit], 0.0
            while node.parent_node is not None:
                if node not in index:
                    index[node] = len(lengths)
                    lengths.append(node.edge.length or 0.0)
                e = index[node]
                present[e][s] += 1
                weight[e][s] += share
                path += node.edge.length or 0.0
                node = node.parent_node
            depth[unit] = path
            tip_depth_weighted[s] += share * path

    n_edges = len(lengths)
    pres = np.zeros((n_edges, len(names)))
    wt = np.zeros((n_edges, len(names)))
    for e in range(n_edges):
        pres[e] = present[e]
        wt[e] = weight[e]
    table = EdgeTable(np.array(lengths), pres, wt, depth, names, richness, tip_depth_weighted, sum_sq)
    return table, missing


def alpha(table: EdgeTable) -> list[dict]:
    """Richness, Faith's PD and MPD (presence and abundance-weighted) per sample."""
    l = table.lengths[:, None]
    k = table.richness[None, :]
    covered = table.present > 0
    pd = (l * covered).sum(axis=0)
    pairs = k * (k - 1) / 2
    sep = (l * table.present * (k - table.present)).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        mpd = np.where(pairs[0] > 0, sep / pairs[0], np.nan)
        # weighted: sum over ordered pairs i != j of p_i p_j d_ij, over sum of p_i p_j
        num = (l * 2 * table.weight * (1 - table.weight)).sum(axis=0)
        denom = 1 - table.sum_sq
        mpd_w = np.where(denom > 1e-12, num / np.where(denom > 1e-12, denom, 1.0), np.nan)
    out = []
    for s, name in enumerate(table.samples):
        out.append(
            {
                "sample": name,
                "richness": int(table.richness[s]),
                "faith_pd": float(pd[s]),
                "mpd": float(mpd[s]),
                "mpd_abundance": float(mpd_w[s]),
            }
        )
    return out


def beta(table: EdgeTable) -> list[dict]:
    """Pairwise unweighted UniFrac, weighted UniFrac and PhyloSor."""
    l = table.lengths
    covered = (table.present > 0).astype(float)
    bl = (l[:, None] * covered).sum(axis=0)
    shared = covered.T @ (l[:, None] * covered)
    rows = []
    n = len(table.samples)
    for i in range(n):
        diff = np.abs(table.weight - table.weight[:, [i]])
        weighted_num = (l[:, None] * diff).sum(axis=0)
        for j in range(i + 1, n):
            union = bl[i] + bl[j] - shared[i, j]
            unweighted = (union - shared[i, j]) / union if union > 0 else np.nan
            denom = table.tip_depth_weighted[i] + table.tip_depth_weighted[j]
            weighted = weighted_num[j] / denom if denom > 0 else np.nan
            mean_bl = 0.5 * (bl[i] + bl[j])
            phylosor = shared[i, j] / mean_bl if mean_bl > 0 else np.nan
            rows.append(
                {
                    "sample_a": table.samples[i],
                    "sample_b": table.samples[j],
                    "unifrac_unweighted": float(unweighted),
                    "unifrac_weighted": float(weighted),
                    "phylosor": float(phylosor),
                }
            )
    return rows
