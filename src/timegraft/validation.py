"""How far does a barcode distance track a trusted tree? The informative window.

Pairs of units are sampled per stratum, by the deepest taxonomic rank they
share (same genus, same family, same order, ...). For each pair, each
distance method is compared with the patristic distance in a reference tree.
Per method, the result is:

- the rank correlation with the reference, per stratum and overall;
- binned medians of the distance against the reference distance, as in a
  saturation plot;
- the plateau onset: the reference distance at which the binned median first
  reaches 90% of its maximum, a simple marker for the end of the window.

Neither the barcode distance nor a barcode-based reference tree is ground
truth; use a reference tree built from independent data where possible.
"""

from __future__ import annotations

import random
from collections import defaultdict
from itertools import combinations

import dendropy
import numpy as np

from . import distances

PAIR_COLUMNS_BASE = ("unit_a", "unit_b", "stratum", "reference")


def stratum(a: dict, b: dict, ranks: list[str]) -> str:
    """Name of the lowest rank two units share, e.g. 'genus', or 'above_<top rank>'."""
    for rank in ranks:
        if a.get(rank) and a.get(rank) == b.get(rank):
            return rank
    return f"above_{ranks[-1]}" if ranks else "any"


def sample_pairs(units: dict[str, dict], ranks: list[str], per_stratum: int, seed: int = 1) -> list[tuple[str, str, str]]:
    """Up to ``per_stratum`` pairs of unit ids per stratum, reproducibly."""
    by_stratum = defaultdict(list)
    for a, b in combinations(sorted(units), 2):
        by_stratum[stratum(units[a], units[b], ranks)].append((a, b))
    rng = random.Random(seed)
    out = []
    for name in sorted(by_stratum):
        pairs = by_stratum[name]
        if len(pairs) > per_stratum:
            pairs = rng.sample(pairs, per_stratum)
        out.extend((a, b, name) for a, b in sorted(pairs))
    return out


def patristic_function(tree: dendropy.Tree):
    pdm = tree.phylogenetic_distance_matrix()
    ns = tree.taxon_namespace
    return lambda a, b: pdm(ns.get_taxon(label=a), ns.get_taxon(label=b))


def spearman(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    rx = _ranks(x)
    ry = _ranks(y)
    return float(np.corrcoef(rx, ry)[0, 1])


def _ranks(v: np.ndarray) -> np.ndarray:
    order = np.argsort(v, kind="mergesort")
    ranks = np.empty(len(v))
    ranks[order] = np.arange(len(v))
    # average ranks for ties
    for value in np.unique(v):
        tied = v == value
        if tied.sum() > 1:
            ranks[tied] = ranks[tied].mean()
    return ranks


def binned_medians(reference, distance, n_bins: int = 10) -> list[dict]:
    """Medians and interquartile range of ``distance`` in quantile bins of ``reference``."""
    reference = np.asarray(reference, float)
    distance = np.asarray(distance, float)
    edges = np.unique(np.quantile(reference, np.linspace(0, 1, n_bins + 1)))
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        inside = (reference >= lo) & (reference <= hi if hi == edges[-1] else reference < hi)
        if inside.sum() == 0:
            continue
        d = distance[inside]
        rows.append(
            {
                "reference_median": float(np.median(reference[inside])),
                "median": float(np.median(d)),
                "q25": float(np.quantile(d, 0.25)),
                "q75": float(np.quantile(d, 0.75)),
                "n": int(inside.sum()),
            }
        )
    return rows


def plateau_onset(bins: list[dict], fraction: float = 0.9) -> float | None:
    """Reference distance at which the binned median first reaches ``fraction`` of its maximum.

    Beyond this point the barcode distance hardly grows with the reference
    distance, so it marks the end of the informative window. Returns None
    when there are too few bins to tell.
    """
    if len(bins) < 3:
        return None
    top = max(b["median"] for b in bins)
    for b in bins:
        if b["median"] >= fraction * top:
            return b["reference_median"]
    return None


def compare(
    tree: dendropy.Tree,
    units: dict[str, dict],
    methods: dict[str, dict],
    ranks: list[str],
    per_stratum: int = 200,
    seed: int = 1,
) -> tuple[list[dict], list[dict], dict[str, list[dict]]]:
    """Run the comparison.

    ``units`` maps unit id to a dict with ``sequence`` and the rank columns;
    only units that are tips of ``tree`` are used. ``methods`` maps a label to
    the options of a distance method (``method`` plus its keyword arguments).

    Returns the pair table, a summary per method and stratum, and the binned
    medians per method.
    """
    tips = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    units = {u: v for u, v in units.items() if u in tips}
    pairs = sample_pairs(units, ranks, per_stratum, seed)
    patristic = patristic_function(tree)
    rows = []
    for a, b, name in pairs:
        row = {"unit_a": a, "unit_b": b, "stratum": name, "reference": patristic(a, b)}
        for label, options in methods.items():
            opts = dict(options)
            func = distances.get_method(opts.pop("method"))
            row[label] = func(units[a]["sequence"], units[b]["sequence"], **opts)
        rows.append(row)

    summary, bins = [], {}
    strata = sorted({r["stratum"] for r in rows}, key=lambda s: ranks.index(s) if s in ranks else len(ranks))
    for label in methods:
        ref = [r["reference"] for r in rows]
        dist = [r[label] for r in rows]
        bins[label] = binned_medians(ref, dist) if rows else []
        point = plateau_onset(bins[label])
        summary.append(
            {"method": label, "stratum": "all", "n_pairs": len(rows), "spearman": spearman(ref, dist),
             "plateau_onset": "" if point is None else point}
        )
        for name in strata:
            sub = [r for r in rows if r["stratum"] == name]
            summary.append(
                {"method": label, "stratum": name, "n_pairs": len(sub),
                 "spearman": spearman([r["reference"] for r in sub], [r[label] for r in sub]),
                 "plateau_onset": ""}
            )
    return rows, summary, bins


def plot(rows: list[dict], bins: dict[str, list[dict]], ranks: list[str], path: str) -> None:
    """One panel per method: pairs coloured by stratum, binned medians on top."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    methods = list(bins)
    fig, axes = plt.subplots(1, len(methods), figsize=(5.5 * len(methods), 4.5), squeeze=False)
    strata = sorted({r["stratum"] for r in rows}, key=lambda s: ranks.index(s) if s in ranks else len(ranks))
    colours = plt.get_cmap("viridis")(np.linspace(0.1, 0.9, max(1, len(strata))))
    for ax, label in zip(axes[0], methods):
        for colour, name in zip(colours, strata):
            sub = [r for r in rows if r["stratum"] == name]
            ax.scatter([r["reference"] for r in sub], [r[label] for r in sub], s=8, alpha=0.6,
                       color=colour, label=f"same {name}" if name in ranks else name)
        b = bins[label]
        if b:
            x = [v["reference_median"] for v in b]
            ax.errorbar(x, [v["median"] for v in b],
                        yerr=[[v["median"] - v["q25"] for v in b], [v["q75"] - v["median"] for v in b]],
                        color="#c0392b", marker="o", capsize=3, label="binned median (IQR)")
        ax.set_xlabel("reference tree distance")
        ax.set_ylabel(f"{label} distance")
        ax.set_title(label)
        ax.grid(alpha=0.3)
    axes[0][0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
