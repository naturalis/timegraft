"""Exploratory: prune suspect units from the rooted genus trees and recompute within-genus distances.

Usage: prune_impact.py DIR PROVENANCE_BEFORE
  DIR holds trees.tsv (../dating-benchmark/trees.tsv.gz, unpacked) and unit_matches.tsv (unpacked)
"""
import csv, sys, math, statistics as st
from pathlib import Path
import dendropy
csv.field_size_limit(10**9)
out, prov_path = Path(sys.argv[1]), sys.argv[2]
def read(p):
    with open(p) as fh: return list(csv.DictReader(fh, delimiter="\t"))
def f(x): return float(x) if x else 0.0
m = {r["unit_id"]: r for r in read(out / "unit_matches.tsv")}
other = {u for u, r in m.items() if f(r["best_other"]) - f(r["best_same_genus"]) >= 0.02}
isolated = {u for u, r in m.items() if r["genus_units"] != "1" and f(r["best_same_genus"]) < 0.80}
suspect = other | isolated
print(f"suspect units: nearest to another genus by 2% or more {len(other)}; no congener at 80% identity {len(isolated)}; either {len(suspect)}")
prov = {r["genus"]: r for r in read(prov_path) if r["status"] == "subtree"}
def stats(tree):
    leaves = tree.leaf_nodes()
    d = [l.distance_from_root() for l in leaves]
    return st.mean(d), 2 * max(d) if len(d) > 1 else 0.0   # mean root-to-tip, and twice the deepest tip as a stand-in for the largest distance
rows = []
for t in read(out / "trees.tsv"):
    tree = dendropy.Tree.get(data=t["newick"], schema="newick", rooting="force-rooted", preserve_underscores=True)
    labels = [l.taxon.label for l in tree.leaf_node_iter()]
    bad = [x for x in labels if x in suspect]
    outs = [float(x.rsplit("=", 1)[1]) for x in t["d_out"].split(";") if x]
    row = {"genus": t["genus"], "n": len(labels), "bad": len(bad), "mrt": float(t["mrt"]), "dmax": float(t["d_in_max"]),
           "out": st.median(outs) if outs else None, "stem": float(t["stem_age"])}
    if bad and len(labels) - len(bad) >= 2:
        tree.prune_taxa_with_labels(bad)
        # after pruning the root may have one child: move to the MRCA of what is left
        node = tree.seed_node
        while len(node.child_nodes()) == 1: node = node.child_nodes()[0]
        sub = dendropy.Tree(seed_node=node.extract_subtree() if hasattr(node, "extract_subtree") else node)
        d = []
        for l in node.leaf_iter():
            x, n2 = 0.0, l
            while n2 is not node: x += n2.edge.length or 0.0; n2 = n2.parent_node
            d.append(x)
        row["mrt_clean"] = st.mean(d); row["left"] = len(d)
    elif bad:
        row["mrt_clean"] = None; row["left"] = len(labels) - len(bad)
    else:
        row["mrt_clean"] = row["mrt"]; row["left"] = len(labels)
    rows.append(row)
aff = [r for r in rows if r["bad"]]
print(f"genera with a subtree: {len(rows)}; with suspect units: {len(aff)}; left with fewer than two units: {sum(r['mrt_clean'] is None for r in rows)}")
ok = [r for r in aff if r["mrt_clean"] is not None and r["mrt"] > 0]
ratio = [r["mrt_clean"] / r["mrt"] for r in ok]
q = st.quantiles(ratio, n=4)
print(f"mean root-to-tip after removing them, relative to before ({len(ok)} genera): median {st.median(ratio):.2f}, quartiles {q[0]:.2f}-{q[2]:.2f}; "
      f"falls by more than a quarter in {sum(x < .75 for x in ratio)}, by more than half in {sum(x < .5 for x in ratio)}")
allr = [r for r in rows if r["mrt_clean"] is not None and r["mrt"] > 0]
b = [r["mrt"] for r in allr]; c = [r["mrt_clean"] for r in allr]
print(f"all genera: median root-to-tip {st.median(b):.3f} -> {st.median(c):.3f}; 10-90%: {sorted(b)[len(b)//10]:.3f}-{sorted(b)[9*len(b)//10]:.3f} -> {sorted(c)[len(c)//10]:.3f}-{sorted(c)[9*len(c)//10]:.3f}; "
      f"below 0.05: {sum(x < .05 for x in b)/len(b):.0%} -> {sum(x < .05 for x in c)/len(c):.0%}")
def spearman(x, y):
    def rk(v):
        o = sorted(range(len(v)), key=v.__getitem__); r = [0.0] * len(v)
        for i, j in enumerate(o): r[j] = i
        return r
    a, bb = rk(x), rk(y); ma, mb = st.mean(a), st.mean(bb)
    return sum((p - ma) * (q - mb) for p, q in zip(a, bb)) / math.sqrt(sum((p - ma) ** 2 for p in a) * sum((q - mb) ** 2 for q in bb))
tt = [(float(prov[r["genus"]]["timetree_crown"]), r) for r in allr if prov[r["genus"]]["timetree_crown"] and float(prov[r["genus"]]["timetree_crown"]) < r["stem"] and r["mrt_clean"] > 0]
print(f"\nAgainst TimeTree crown age ({len(tt)} genera): median root-to-tip, before -> after removing suspect units")
for lo, hi in [(0, 5), (5, 10), (10, 20), (20, 40), (40, 80), (80, 1000)]:
    s = [r for a, r in tt if lo <= a < hi]
    if s: print(f"   {lo:>3}-{hi:<5}{len(s):>4}  {st.median(r['mrt'] for r in s):.3f} -> {st.median(r['mrt_clean'] for r in s):.3f}   genera with suspect units: {sum(r['bad'] > 0 for r in s)}")
for name, key in (("before", "mrt"), ("after", "mrt_clean")):
    lr = [math.log2(r[key] / 0.0115 / a) for a, r in tt]
    print(f"   barcode minimum {name}: median estimate/TimeTree {2 ** st.median(lr):.2f}; rank correlation with TimeTree {spearman([r[key] for a, r in tt], [a for a, r in tt]):.2f}")
