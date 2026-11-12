"""How much does sample PD depend on genus crown ages? Random samples on the tree of a run,
with every genus crown at its lower bound, as built, at the midpoint of its bounds, and at its upper bound.

Usage: pd_sensitivity.py RUN_DIR   (the output folder of a run: timegraft.tre, provenance.tsv, subtrees/, genera/)
"""
import csv
import glob
import random
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

import dendropy

from timegraft.graft import rescale_crowns

csv.field_size_limit(10**9)
run = Path(sys.argv[1])


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


genus_units = {}
for path in glob.glob(str(run / "subtrees" / "shard_*.tsv")):
    for row in read(path):
        sub = dendropy.Tree.get(data=row["newick"], schema="newick", preserve_underscores=True)
        genus_units[row["genus"]] = [leaf.taxon.label for leaf in sub.leaf_node_iter()]
by_genus = {g: list(u) for g, u in genus_units.items()}
for row in read(run / "genera" / "singles.tsv"):
    by_genus[row["genus"]] = [row["unit_id"]]
prov = [r for r in read(run / "provenance.tsv") if r["status"] == "subtree"]
placed = {r["unit_id"] for r in read(run / "units_placed.tsv") if r["status"] in ("placed", "genus_stem")}
targets = {
    "lower bound": {r["genus"]: float(r["crown_lower"]) for r in prov},
    "as built": {r["genus"]: float(r["crown_age"]) for r in prov},
    "midpoint": {r["genus"]: (float(r["crown_lower"]) * float(r["crown_upper"])) ** 0.5 for r in prov},
    "upper bound": {r["genus"]: float(r["crown_upper"]) for r in prov},
}


def index(tree):
    parent, length, tip = {}, {}, {}
    for i, node in enumerate(tree.preorder_node_iter()):
        node._i = i
    for node in tree.preorder_node_iter():
        if node.parent_node is not None:
            parent[node._i] = node.parent_node._i
            length[node._i] = node.edge.length or 0.0
        if node.is_leaf():
            tip[node.taxon.label] = node._i
    return parent, length, tip


def faith_pd(tree, labels):
    parent, length, tip = tree
    seen, total = set(), 0.0
    for label in labels:
        i = tip[label]
        while i in parent and i not in seen:
            seen.add(i)
            total += length[i]
            i = parent[i]
    return total


def spearman(x, y):
    def rank(v):
        o = sorted(range(len(v)), key=v.__getitem__)
        r = [0.0] * len(v)
        for i, j in enumerate(o):
            r[j] = i
        return r
    a, b = rank(x), rank(y)
    ma, mb = st.mean(a), st.mean(b)
    return sum((p - ma) * (q - mb) for p, q in zip(a, b)) / (
        sum((p - ma) ** 2 for p in a) * sum((q - mb) ** 2 for q in b)) ** 0.5


trees = {}
for name, ages in targets.items():
    tree = dendropy.Tree.get(path=str(run / "timegraft.tre"), schema="newick", preserve_underscores=True,
                             rooting="force-rooted")
    rescale_crowns(tree, genus_units, ages, placed)
    total = sum(e.length or 0 for e in tree.postorder_edge_iter())
    crown = st.median(ages.values())
    print(f"{name:<12} median crown age {crown:6.1f} Ma; total branch length of the tree {total:,.0f} Ma")
    trees[name] = index(tree)
tips = trees["lower bound"][2]
by_genus = {g: [u for u in units if u in tips] for g, units in by_genus.items()}
genera = [g for g, units in by_genus.items() if units]

random.seed(1)
print("\nFaith's PD of random samples, relative to the tree with crowns at the lower bound (200 samples each)")
print(f"{'genera x units per genus':<28}{'as built':>10}{'midpoint':>10}{'upper bound':>13}{'rank corr. lower~built':>24}{'lower~upper':>13}")
for n_genera, per_genus in [(100, 1), (100, 2), (100, 5), (30, 3), (300, 3), (100, 20)]:
    pd = defaultdict(list)
    for _ in range(200):
        sample = []
        for g in random.sample(genera, n_genera):
            sample += random.sample(by_genus[g], min(per_genus, len(by_genus[g])))
        for name, tree in trees.items():
            pd[name].append(faith_pd(tree, sample))
    built = st.median(m / l for m, l in zip(pd["as built"], pd["lower bound"]))
    mid = st.median(m / l for m, l in zip(pd["midpoint"], pd["lower bound"]))
    up = st.median(u / l for u, l in zip(pd["upper bound"], pd["lower bound"]))
    print(f"{n_genera:>4} x {per_genus:<21}{built:>10.3f}{mid:>10.3f}{up:>13.3f}"
          f"{spearman(pd['lower bound'], pd['as built']):>24.3f}{spearman(pd['lower bound'], pd['upper bound']):>13.3f}")
