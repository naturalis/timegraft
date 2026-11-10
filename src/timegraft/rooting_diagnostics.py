"""Diagnostics for outgroup rooting of genus subtrees (issue #19).

Reads the outputs of a finished run and, for every genus subtree, rebuilds the
unrooted BioNJ tree exactly as the pipeline did, then reports:

* whether the ingroup is a split with respect to the outgroup (the pipeline's test),
  and whether that reproduces the rooting recorded in the provenance;
* the rogue count: the fewest ingroup units whose removal makes the ingroup a split,
  and the fewest outgroup units whose removal does (outgroup nested in the genus);
* distance statistics that could explain failures: ingroup-to-outgroup distances,
  the share of saturated (capped) pairs, the smallest pairwise overlap;
* the root, crown age (barcodes and rate) and relative within-genus PD under outgroup,
  midpoint and MAD (minimal ancestor deviation) rooting of the same ingroup tree;
* optionally, for the genera that failed, the same test with other outgroup choices
  and a gamma-corrected distance.

Usage::

    python -m timegraft.rooting_diagnostics --run results/iberia --threads 8

``--run`` is the run's output folder. Settings come from ``config/config.yaml``
unless other config files are given; when the sequences turn out to be
unaligned and the marker says ``align: none``, they are aligned with MAFFT.
"""

from __future__ import annotations

import argparse
import math
import statistics
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import dendropy
import numpy as np

from . import distances
from .backbone import read_backbone, stem_ages
from .bionj import bionj
from .build import prepare_sequences
from .dating import mean_root_to_tip, mpl_ultrametric
from .rooting import midpoint, root_on_outgroup
from .seqio import read_tsv, write_tsv

NUC = {"A": 0, "C": 1, "G": 2, "T": 3}
PURINE = np.array([True, False, True, False])


# ---------------------------------------------------------------- config

def _merge(base: dict, extra: dict) -> dict:
    out = dict(base)
    for key, value in extra.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _merge(out[key], value)
        else:
            out[key] = value
    return out


def load_config(paths: list[str]) -> dict:
    import yaml

    cfg: dict = {}
    for path in paths:
        with open(path) as fh:
            cfg = _merge(cfg, yaml.safe_load(fh) or {})
    return cfg


# ---------------------------------------------------------------- distances

def encode(seqs: list[str]) -> np.ndarray:
    """Aligned sequences as an int8 array, -1 for gaps and ambiguity codes."""
    width = max(len(s) for s in seqs)
    arr = np.full((len(seqs), width), -1, dtype=np.int8)
    for i, seq in enumerate(seqs):
        for j, ch in enumerate(seq.upper()):
            arr[i, j] = NUC.get(ch, -1)
    return arr


def pq_matrices(arr: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transition and transversion proportions and comparable sites, pairwise deletion."""
    n = arr.shape[0]
    P = np.zeros((n, n))
    Q = np.zeros((n, n))
    S = np.zeros((n, n), dtype=int)
    valid = arr >= 0
    purine = np.where(valid, PURINE[np.clip(arr, 0, 3)], False)
    for i in range(n):
        both = valid[i] & valid
        diff = both & (arr[i] != arr)
        ts = diff & (purine[i] == purine)
        sites = both.sum(axis=1)
        S[i] = sites
        with np.errstate(invalid="ignore", divide="ignore"):
            P[i] = np.where(sites > 0, ts.sum(axis=1) / sites, np.nan)
            Q[i] = np.where(sites > 0, (diff.sum(axis=1) - ts.sum(axis=1)) / sites, np.nan)
    return P, Q, S


def k2p_from_pq(P, Q, alpha: float | None = None, cap: float = distances.DEFAULT_MAX_DISTANCE) -> np.ndarray:
    w1 = 1.0 - 2.0 * P - Q
    w2 = 1.0 - 2.0 * Q
    ok = (w1 > 0) & (w2 > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        if alpha is None:
            d = -0.5 * np.log(w1) - 0.25 * np.log(w2)
        else:  # Jin and Nei (1990)
            d = (alpha / 2.0) * (w1 ** (-1.0 / alpha) + 0.5 * w2 ** (-1.0 / alpha) - 1.5)
    d = np.where(ok, np.minimum(d, cap), cap)
    np.fill_diagonal(d, 0.0)
    return d


# ---------------------------------------------------------------- tree helpers

def leaf_sets(tree: dendropy.Tree) -> dict[dendropy.Node, frozenset[str]]:
    sets = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            sets[node] = frozenset([node.taxon.label])
        else:
            sets[node] = frozenset().union(*(sets[c] for c in node.child_node_iter()))
    return sets


def copy_tree(tree: dendropy.Tree) -> dendropy.Tree:
    return dendropy.Tree.get(data=tree.as_string(schema="newick", suppress_rooting=True),
                             schema="newick", preserve_underscores=True)


def rogue_counts(tree: dendropy.Tree, ingroup: frozenset, outgroup: frozenset) -> dict:
    """Fewest units to remove for the ingroup to be one side of an edge.

    ``rogue_in``: removing ingroup units only; ``rogue_out``: removing outgroup units
    only; ``min_violation``: either. Sides are taken over every edge of the tree.
    """
    labels = ingroup | outgroup
    best_in = (len(ingroup), None)
    best_out = (len(outgroup) + 1, None)
    best_mv = (len(labels), None)
    for node, below in leaf_sets(tree).items():
        if node.parent_node is None:
            continue
        for side in (below, labels - below):
            ins = len(ingroup & side)
            if ins == 0:
                continue
            miss_in = len(ingroup - side)
            extra_out = len(outgroup & side)
            if extra_out == 0 and miss_in < best_in[0]:
                best_in = (miss_in, ingroup - side)
            if miss_in == 0 and extra_out < best_out[0]:
                best_out = (extra_out, outgroup & side)
            if miss_in + extra_out < best_mv[0]:
                best_mv = (miss_in + extra_out, side)
    return {
        "rogue_in": best_in[0],
        "rogue_in_units": ";".join(sorted(best_in[1] or [])),
        "rogue_out": best_out[0] if best_out[1] is not None else "",
        "min_violation": best_mv[0],
    }


def _all_node_distances(tree: dendropy.Tree):
    nodes = list(tree.preorder_node_iter())
    index = {node: i for i, node in enumerate(nodes)}
    adj = defaultdict(list)
    for node in nodes:
        if node.parent_node is not None:
            length = max(0.0, node.edge.length or 0.0)
            adj[index[node]].append((index[node.parent_node], length))
            adj[index[node.parent_node]].append((index[node], length))
    leaves = [index[n] for n in nodes if n.is_leaf()]
    L = np.zeros((len(leaves), len(nodes)))
    for row, start in enumerate(leaves):
        dist = {start: 0.0}
        stack = [start]
        while stack:
            u = stack.pop()
            for v, length in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + length
                    stack.append(v)
        L[row] = [dist[i] for i in range(len(nodes))]
    return nodes, index, leaves, L


def mad_root(tree: dendropy.Tree) -> tuple[frozenset[str], float, float, float]:
    """Minimal ancestor deviation root (Tria et al. 2017) of an unrooted tree.

    Returns the leaf set on the child side of the root edge, the distance of the
    root from that child, the MAD score and the root ambiguity index (best score
    over second best; close to 1 means another root is about as good).
    """
    nodes, index, leaves, L = _all_node_distances(tree)
    n = len(leaves)
    leaf_col = np.array(leaves)
    D = L[:, leaf_col]
    iu = np.triu_indices(n, 1)
    usable = D[iu] > 0
    n_pairs = int(usable.sum())
    sets = leaf_sets(tree)
    label_row = {nodes[c].taxon.label: r for r, c in enumerate(leaves)}

    def within(rows: np.ndarray, col: int) -> float:
        if len(rows) < 2:
            return 0.0
        lu = L[rows, col]
        sub = D[np.ix_(rows, rows)]
        k = np.triu_indices(len(rows), 1)
        d = sub[k]
        diff = (lu[:, None] - lu[None, :])[k]
        keep = d > 0
        return float(np.sum((diff[keep] / d[keep]) ** 2))

    scores = []
    for node in nodes:
        parent = node.parent_node
        if parent is None:
            continue
        length = max(0.0, node.edge.length or 0.0)
        u, v = index[parent], index[node]
        B = np.array(sorted(label_row[l] for l in sets[node]))
        mask = np.ones(n, bool)
        mask[B] = False
        A = np.nonzero(mask)[0]
        if len(A) == 0:
            continue
        ss = within(A, u) + within(B, v)
        dbc = D[np.ix_(A, B)]
        ok = dbc > 0
        a = (2.0 * L[A, u][:, None] / np.where(ok, dbc, 1.0)) - 1.0
        w = np.where(ok, 1.0 / np.where(ok, dbc, 1.0), 0.0)
        a = np.where(ok, a, 0.0)
        denom = 2.0 * np.sum(w ** 2)
        x = -np.sum(a * w) / denom if denom > 0 else 0.0
        x = min(max(x, 0.0), length)
        ss += float(np.sum(np.where(ok, (a + 2.0 * x * w) ** 2, 0.0)))
        scores.append((math.sqrt(ss / max(n_pairs, 1)), sets[node], length - x))
    scores.sort(key=lambda t: t[0])
    best = scores[0]
    rai = best[0] / scores[1][0] if len(scores) > 1 and scores[1][0] > 0 else float("nan")
    return best[1], best[2], best[0], rai


def reroot(tree: dendropy.Tree, child_side: frozenset[str], to_child: float) -> dendropy.Tree:
    labels = frozenset(l.taxon.label for l in tree.leaf_node_iter())
    for node, below in leaf_sets(tree).items():
        if node.parent_node is None:
            continue
        if below == child_side or below == labels - child_side:
            length = max(0.0, node.edge.length or 0.0)
            d = to_child if below == child_side else length - to_child
            tree.reroot_at_edge(node.edge, length1=length - d, length2=d, update_bipartitions=False)
            tree.is_rooted = True
            tree.suppress_unifurcations()
            return tree
    raise ValueError("root edge not found")


def root_split(tree: dendropy.Tree, keep: frozenset[str] | None = None) -> frozenset | None:
    """The root bipartition as a frozenset of two frozensets, optionally restricted."""
    parts = [frozenset(l.taxon.label for l in c.leaf_iter()) for c in tree.seed_node.child_node_iter()]
    if keep is not None:
        parts = [p & keep for p in parts]
    parts = [p for p in parts if p]
    if len(parts) < 2:
        return None
    if len(parts) > 2:  # polytomy at the root: compare on the smallest part
        parts = sorted(parts, key=len)
        parts = [parts[0], frozenset().union(*parts[1:])]
    return frozenset(parts)


def describe(tree: dendropy.Tree, rate: float) -> dict:
    depth = mean_root_to_tip(tree)
    crown = depth / rate if rate else float("nan")
    dated = copy_tree(tree)  # a newick round trip keeps the root
    dated.is_rooted = True
    mpl_ultrametric(dated, 1.0)
    pd_rel = sum(e.length or 0.0 for e in dated.postorder_edge_iter() if e.tail_node is not None)
    return {"crown": crown, "pd_rel": pd_rel}


# ---------------------------------------------------------------- per genus

def analyse(task: dict) -> dict:
    genus, ingroup_items, outgroup_items = task["genus"], task["in"], task["out"]
    cfg = task["marker_cfg"]
    rate = cfg.get("crown_age", {}).get("rate") or 0.0
    method = cfg["distance"]["method"]
    options = {k: v for k, v in cfg["distance"].items() if k != "method"}
    ingroup = frozenset(n for n, _ in ingroup_items)
    outgroup = frozenset(n for n, _ in outgroup_items)

    items = prepare_sequences(ingroup_items + outgroup_items, method, cfg.get("align", "none"))
    labels = [n for n, _ in items]
    dist = distances.matrix([s for _, s in items], method=method, **options)
    row = {"genus": genus, "n_in": len(ingroup), "n_out": len(outgroup)}

    # distance statistics
    in_idx = [i for i, l in enumerate(labels) if l in ingroup]
    out_idx = [i for i, l in enumerate(labels) if l in outgroup]
    cap = options.get("max_distance", distances.DEFAULT_MAX_DISTANCE)
    iu = np.triu_indices(len(in_idx), 1)
    within = dist[np.ix_(in_idx, in_idx)][iu]
    row["median_in_in"] = float(np.median(within)) if within.size else 0.0
    row["max_in_in"] = float(np.max(within)) if within.size else 0.0
    if out_idx:
        cross = dist[np.ix_(in_idx, out_idx)]
        row["min_in_out"] = float(cross.min())
        row["mean_in_out"] = float(cross.mean())
        row["frac_capped_in_out"] = float(np.mean(cross >= cap))
    row["frac_capped_all"] = float(np.mean(dist[np.triu_indices(len(labels), 1)] >= cap)) if len(labels) > 1 else 0.0
    arr = None
    if distances.NEEDS_ALIGNMENT.get(method, False):
        arr = encode([s for _, s in items])
        P, Q, S = pq_matrices(arr)
        off = ~np.eye(len(labels), dtype=bool)
        row["min_overlap"] = int(S[off].min()) if off.any() else ""
        row["min_seq_sites"] = int((arr >= 0).sum(axis=1).min())
        if out_idx:
            row["min_overlap_in_out"] = int(S[np.ix_(in_idx, out_idx)].min())

    tree = bionj(labels, dist)
    row.update(_rooting_comparison(tree, ingroup, outgroup, rate))

    # variants for genera where the outgroup failed
    if task.get("variants") and outgroup and not row["split_ok"]:
        for name, out_items in task["variants"].items():
            row.update(_variant(name, ingroup_items, out_items, cfg, rate))
        if arr is not None:
            for alpha in task.get("gammas", []):
                d = k2p_from_pq(P, Q, alpha=alpha, cap=cap)
                t = bionj(labels, d)
                rc = rogue_counts(t, ingroup, outgroup)
                row[f"gamma{alpha}_split_ok"] = rc["rogue_in"] == 0
                row[f"gamma{alpha}_rogue_in"] = rc["rogue_in"]
    return row


def _rooting_comparison(tree, ingroup, outgroup, rate) -> dict:
    row = {}
    if outgroup:
        rc = rogue_counts(tree, ingroup, outgroup)
        row.update(rc)
        row["split_ok"] = rc["rogue_in"] == 0
    else:
        row["split_ok"] = False

    # outgroup-rooted ingroup tree, if it exists; otherwise after pruning the rogues
    keep = ingroup
    og_tree = None
    if outgroup:
        if row["split_ok"]:
            og_tree, how = root_on_outgroup(copy_tree(tree), ingroup, outgroup)
            row["reproduced"] = how
        else:
            row["reproduced"] = "midpoint_outgroup_conflict"
            rogues = frozenset(row["rogue_in_units"].split(";")) if row["rogue_in_units"] else frozenset()
            keep = ingroup - rogues
            if len(keep) >= 2:
                pruned = copy_tree(tree)
                pruned.prune_taxa_with_labels(list(rogues))
                og_tree, _ = root_on_outgroup(pruned, keep, outgroup)
    else:
        row["reproduced"] = "midpoint"

    # ingroup-only tree, as the fallback sees it
    ing = copy_tree(tree)
    ing.prune_taxa_with_labels([l for l in outgroup])
    if len(ingroup) < 3:
        return row
    mp = midpoint(copy_tree(ing))
    side, to_child, score, rai = mad_root(ing)
    md = reroot(copy_tree(ing), side, to_child)
    row["mad_score"], row["mad_rai"] = score, rai

    descs = {"midpoint": describe(mp, rate), "mad": describe(md, rate)}
    if og_tree is not None and len(keep) >= 3:
        descs["outgroup"] = describe(og_tree, rate)
        ref = root_split(og_tree, keep)
        row["midpoint_same_root"] = root_split(mp, keep) == ref
        row["mad_same_root"] = root_split(md, keep) == ref
    row["midpoint_mad_same_root"] = root_split(mp) == root_split(md)
    for name, d in descs.items():
        row[f"crown_{name}"] = d["crown"]
        row[f"pd_rel_{name}"] = d["pd_rel"]
    return row


def _variant(name, ingroup_items, out_items, cfg, rate) -> dict:
    if not out_items:
        return {f"{name}_split_ok": "", f"{name}_rogue_in": ""}
    method = cfg["distance"]["method"]
    options = {k: v for k, v in cfg["distance"].items() if k != "method"}
    items = prepare_sequences(ingroup_items + out_items, method, cfg.get("align", "none"))
    labels = [n for n, _ in items]
    tree = bionj(labels, distances.matrix([s for _, s in items], method=method, **options))
    rc = rogue_counts(tree, frozenset(n for n, _ in ingroup_items), frozenset(n for n, _ in out_items))
    return {f"{name}_split_ok": rc["rogue_in"] == 0, f"{name}_rogue_in": rc["rogue_in"],
            f"{name}_min_violation": rc["min_violation"]}


# ---------------------------------------------------------------- outgroup variants

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


def pick_longest_first_sister(clades, by_genus, n: int) -> list[tuple[str, str]]:
    if not clades:
        return []
    g = clades[0][0]
    rows = sorted(by_genus[g], key=lambda r: (-_sites(r["sequence"]), r["unit_id"]))[:n]
    return [(r["unit_id"], r["sequence"]) for r in rows]


# ---------------------------------------------------------------- report

def _rate(rows, key="split_ok"):
    vals = [r for r in rows if r.get(key) != ""]
    return (sum(1 for r in vals if not r[key]), len(vals))


def _fmt_rate(fail, total):
    return f"{fail}/{total} ({100 * fail / total:.0f}%)" if total else "0/0"


def _quantile_bins(rows, key, k=4):
    vals = sorted(float(r[key]) for r in rows if r.get(key) not in ("", None))
    if len(vals) < k:
        return []
    edges = [vals[int(len(vals) * i / k)] for i in range(1, k)]
    bins = defaultdict(list)
    for r in rows:
        if r.get(key) in ("", None):
            continue
        v = float(r[key])
        b = sum(v > e for e in edges)
        bins[b].append(r)
    labels = [f"<= {edges[0]:.3g}"] + [f"{edges[i-1]:.3g} to {edges[i]:.3g}" for i in range(1, len(edges))] + [f"> {edges[-1]:.3g}"]
    return [(labels[b], bins[b]) for b in sorted(bins)]


def _ratio_summary(rows, a, b):
    ratios = [float(r[a]) / float(r[b]) for r in rows if r.get(a) not in ("", None) and r.get(b) not in ("", None) and float(r[b]) > 0]
    if not ratios:
        return "n/a"
    q = statistics.quantiles(ratios, n=4) if len(ratios) > 1 else [ratios[0]] * 3
    off = sum(1 for x in ratios if x < 0.8 or x > 1.25)
    return f"median {q[1]:.2f} (IQR {q[0]:.2f} to {q[2]:.2f}); off by more than 25%: {off}/{len(ratios)}"


def report(rows: list[dict], variant_names: list[str], gammas: list[float]) -> str:
    out = ["# Rooting diagnostics", ""]
    with_og = [r for r in rows if r["n_out"] > 0]
    fails = [r for r in with_og if not r["split_ok"]]
    oks = [r for r in with_og if r["split_ok"]]
    mism = [r for r in rows if r.get("recorded") and r["recorded"] != r["reproduced"]]
    out += [f"Genus subtrees: {len(rows)}; with an outgroup: {len(with_og)}; "
            f"outgroup not a split: {_fmt_rate(len(fails), len(with_og))}.",
            f"Rooting recorded in the run but not reproduced here: {len(mism)}"
            + (f" ({', '.join(r['genus'] for r in mism[:10])}{'...' if len(mism) > 10 else ''})" if mism else "") + ".", ""]

    out += ["## How badly do the failures fail?", "",
            "Fewest ingroup units to remove for the ingroup to be a split (rogue_in), and fewest "
            "outgroup units to remove instead (rogue_out: the outgroup nested in the genus).", "",
            "| rogue_in | genera | | rogue_out | genera |", "| --- | --- | --- | --- | --- |"]
    def bucket(v):
        if v == "":
            return "none possible"
        v = int(v)
        return str(v) if v <= 2 else ("3 to 5" if v <= 5 else "more than 5")
    ci = Counter(bucket(r["rogue_in"]) for r in fails)
    co = Counter(bucket(r["rogue_out"]) for r in fails)
    order = ["1", "2", "3 to 5", "more than 5", "none possible"]
    for key in order:
        if ci.get(key) or co.get(key):
            out.append(f"| {key} | {ci.get(key, 0)} | | {key} | {co.get(key, 0)} |")
    frac = [int(r["rogue_in"]) / r["n_in"] for r in fails if r["rogue_in"] != ""]
    if frac:
        out += ["", f"Rogue units as a share of the genus: median {statistics.median(frac):.0%}."]
    out.append("")

    out += ["## What predicts failure?", "", "Failure rate by quartile of each variable.", ""]
    for key in ["n_in", "stem_age", "mean_in_out", "min_in_out", "max_in_in", "frac_capped_in_out", "min_overlap", "min_seq_sites"]:
        bins = _quantile_bins(with_og, key)
        if not bins:
            continue
        out += [f"**{key}**", "", "| bin | failed |", "| --- | --- |"]
        out += [f"| {label} | {_fmt_rate(*_rate(rs))} |" for label, rs in bins]
        out.append("")
    if any(r.get("backbone_status") for r in with_og):
        out += ["**backbone_status**", "", "| status | failed |", "| --- | --- |"]
        groups = defaultdict(list)
        for r in with_og:
            s = r.get("backbone_status", "") or "unknown"
            groups["attached_at_taxon" if s.startswith("attached_at_taxon") else s].append(r)
        out += [f"| {s} | {_fmt_rate(*_rate(rs))} |" for s, rs in sorted(groups.items())]
        out.append("")

    out += ["## How good are the fallbacks?", "",
            "On the genera where the outgroup worked, the ingroup tree was rooted again by midpoint "
            "and by MAD, and compared with the outgroup root (same root split) and with the crown age "
            "from barcodes and rate and the relative within-genus PD (tree made ultrametric with crown 1).", ""]
    comp = [r for r in oks if r.get("midpoint_same_root") not in ("", None)]
    if comp:
        mp = sum(1 for r in comp if r["midpoint_same_root"])
        md = sum(1 for r in comp if r["mad_same_root"])
        out += ["| | midpoint | MAD |", "| --- | --- | --- |",
                f"| same root as outgroup | {mp}/{len(comp)} ({100*mp/len(comp):.0f}%) | {md}/{len(comp)} ({100*md/len(comp):.0f}%) |",
                f"| crown age / outgroup crown | {_ratio_summary(comp, 'crown_midpoint', 'crown_outgroup')} | {_ratio_summary(comp, 'crown_mad', 'crown_outgroup')} |",
                f"| relative PD / outgroup | {_ratio_summary(comp, 'pd_rel_midpoint', 'pd_rel_outgroup')} | {_ratio_summary(comp, 'pd_rel_mad', 'pd_rel_outgroup')} |", ""]
        for n_lo, n_hi in [(3, 5), (6, 20), (21, 10**9)]:
            sub = [r for r in comp if n_lo <= r["n_in"] <= n_hi]
            if sub:
                a = sum(1 for r in sub if r["midpoint_same_root"]); b = sum(1 for r in sub if r["mad_same_root"])
                out.append(f"- {n_lo} to {n_hi if n_hi < 10**9 else 'more'} units: midpoint {a}/{len(sub)}, MAD {b}/{len(sub)}")
        for label, sub in _quantile_bins(comp, "stem_age"):
            a = sum(1 for r in sub if r["midpoint_same_root"]); b = sum(1 for r in sub if r["mad_same_root"])
            out.append(f"- stem age {label} Ma: midpoint {a}/{len(sub)}, MAD {b}/{len(sub)}")
        out += ["", "Agreement is not correctness: the outgroup root can pass the test and still be wrong, "
                "more often in old genera. Falling agreement with stem age points to saturation of the "
                "outgroup distances rather than to poor fallbacks."]
        out.append("")
    fcomp = [r for r in fails if r.get("midpoint_same_root") not in ("", None)]
    if fcomp:
        mp = sum(1 for r in fcomp if r["midpoint_same_root"]); md = sum(1 for r in fcomp if r["mad_same_root"])
        out += [f"On the failed genera, after pruning the rogue units and rooting on the outgroup: "
                f"midpoint agrees in {mp}/{len(fcomp)}, MAD in {md}/{len(fcomp)}.", ""]

    if variant_names or gammas:
        out += ["## Would other outgroups or distances help?", "",
                "Re-run on the failed genera only; 'rescued' means the ingroup is now a split.", "",
                "| variant | rescued | median rogue_in |", "| --- | --- | --- |"]
        for name in variant_names + [f"gamma{a}" for a in gammas]:
            vals = [r for r in fails if r.get(f"{name}_split_ok") not in ("", None)]
            if not vals:
                continue
            ok = sum(1 for r in vals if r[f"{name}_split_ok"])
            med = statistics.median(int(r[f"{name}_rogue_in"]) for r in vals)
            out.append(f"| {name} | {ok}/{len(vals)} | {med:g} |")
        out.append("")
    return "\n".join(out)


# ---------------------------------------------------------------- main

def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", help="the run's output folder (default: outdir from the config)")
    ap.add_argument("--configfile", nargs="+", default=["config/config.yaml"],
                    help="config files, as given to snakemake (default: config/config.yaml)")
    ap.add_argument("--marker", help="marker block to use (default: marker from the config)")
    ap.add_argument("--align", choices=["none", "mafft"], help="override the marker's align setting")
    ap.add_argument("--outdir", default=None, help="default: <run outdir>/diagnostics/rooting")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--genera", nargs="*", help="only these genera")
    ap.add_argument("--no-variants", action="store_true", help="skip the outgroup and distance variants")
    ap.add_argument("--gamma", type=float, nargs="*", default=[0.5], help="gamma shape(s) for K2P+G variants")
    args = ap.parse_args(argv)

    cfg = load_config(args.configfile)
    run = Path(args.run or cfg["outdir"])
    if not (run / "genera").is_dir():
        sys.exit(f"{run} has no genera/ folder; point --run at the run's output folder")
    outdir = Path(args.outdir) if args.outdir else run / "diagnostics" / "rooting"
    outdir.mkdir(parents=True, exist_ok=True)
    marker_cfg = dict(cfg["markers"][args.marker or cfg["marker"]])
    if args.align:
        marker_cfg["align"] = args.align

    groups = defaultdict(lambda: {"in": [], "out": [], "stem": None})
    for path in sorted((run / "genera").glob("shard_*.tsv")):
        for r in read_tsv(path):
            g = groups[r["genus"]]
            g[r["role"]].append((r["unit_id"], r["sequence"]))
            g["stem"] = float(r["stem_age"])
    if (marker_cfg.get("align", "none") == "none" and not args.align
            and distances.NEEDS_ALIGNMENT.get(marker_cfg["distance"]["method"], False)
            and any(len({len(seq) for _, seq in g["in"] + g["out"]}) > 1 for g in groups.values())):
        marker_cfg["align"] = "mafft"
        print("sequences are unaligned: aligning with MAFFT, as the run must have done", file=sys.stderr)
    recorded = {}
    for path in sorted((run / "subtrees").glob("shard_*.tsv")):
        for r in read_tsv(path):
            recorded[r["genus"]] = r
    status = {}
    if (run / "backbone" / "genera.tsv").exists():
        for r in read_tsv(run / "backbone" / "genera.tsv"):
            status[r["genus"]] = r.get("status", "")

    backbone_path = run / "backbone" / "backbone.tre"
    if not backbone_path.exists():
        backbone_path = Path(cfg["backbone"])
    backbone = read_backbone(backbone_path)
    reps = read_tsv(run / "units" / "representatives.tsv")
    by_genus = defaultdict(list)
    for r in reps:
        by_genus[r["genus"]].append(r)
    unit_genus = {r["unit_id"]: r["genus"] for r in reps}
    available = {g for g in by_genus if g in stem_ages(backbone)}

    variant_names = [] if args.no_variants else ["out2_longest", "out5_two_genera", "out5_spread", "out10_spread"]
    tasks = []
    for genus in sorted(groups):
        if args.genera and genus not in args.genera:
            continue
        g = groups[genus]
        variants = None
        if variant_names:
            clades = sister_clades(backbone, genus, available)
            variants = {
                "out2_longest": pick_longest_first_sister(clades, by_genus, 2),
                "out5_two_genera": pick_spread(clades, by_genus, 5, n_genera=2),
                "out5_spread": pick_spread(clades, by_genus, 5),
                "out10_spread": pick_spread(clades, by_genus, 10),
            }
        tasks.append({"genus": genus, "in": g["in"], "out": g["out"], "marker_cfg": marker_cfg,
                      "variants": variants, "gammas": [] if args.no_variants else args.gamma})

    if args.threads > 1:
        with Pool(args.threads) as pool:
            rows = list(pool.imap_unordered(analyse, tasks, chunksize=1))
    else:
        rows = [analyse(t) for t in tasks]

    for row in rows:
        g = groups[row["genus"]]
        row["stem_age"] = g["stem"]
        row["backbone_status"] = status.get(row["genus"], "")
        row["recorded"] = recorded.get(row["genus"], {}).get("rooting", "")
        row["outgroup_genera"] = ";".join(sorted({unit_genus.get(u, "?") for u, _ in g["out"]}))
        for k, v in list(row.items()):
            if isinstance(v, float):
                row[k] = round(v, 6)
    rows.sort(key=lambda r: r["genus"])
    columns = []
    for row in rows:
        for k in row:
            if k not in columns:
                columns.append(k)
    write_tsv(outdir / "per_genus.tsv", rows, columns)
    gammas = [] if args.no_variants else args.gamma
    (outdir / "report.md").write_text(report(rows, variant_names, gammas))
    print(f"wrote {outdir / 'per_genus.tsv'} and {outdir / 'report.md'}", file=sys.stderr)


if __name__ == "__main__":
    main()
