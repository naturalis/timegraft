"""Faith's PD of the same random samples on two trees of the same units (here: before and after screening).

Usage: sample_pd.py TREE_A TREE_B SCREENED   (SCREENED: units/screened.tsv of the run of TREE_B)
"""
import csv
import random
import statistics as st
import sys
from collections import defaultdict

import dendropy

csv.field_size_limit(10**9)
path_a, path_b, screened_path = sys.argv[1:4]


def index(path):
    tree = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True)
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


a, b = index(path_a), index(path_b)
rows = list(csv.DictReader(open(screened_path), delimiter="\t"))
aside = {r["unit_id"] for r in rows if r["decision"] in ("other_genus", "isolated")}
by_genus = defaultdict(list)
for r in rows:
    if r["unit_id"] in a[2] and r["unit_id"] in b[2]:
        by_genus[r["genus"]].append(r["unit_id"])
genera = list(by_genus)
units = [u for us in by_genus.values() for u in us]
random.seed(1)
print("Faith's PD after/before, same samples (200 per design)")
for label, draw in [
    ("150 units at random", lambda: random.sample(units, 150)),
    ("100 genera x 1 unit", lambda: [random.choice(by_genus[g]) for g in random.sample(genera, 100)]),
    ("100 genera x 5 units", lambda: [u for g in random.sample(genera, 100) for u in random.sample(by_genus[g], min(5, len(by_genus[g])))]),
]:
    pa, pb, share = [], [], []
    for _ in range(200):
        s = draw()
        pa.append(faith_pd(a, s))
        pb.append(faith_pd(b, s))
        share.append(sum(u in aside for u in s) / len(s))
    r = [y / x for x, y in zip(pa, pb)]
    print(f"   {label:<22} median {st.median(r):.3f} (range {min(r):.3f}-{max(r):.3f}); rank correlation {spearman(pa, pb):.3f}; "
          f"set-aside units in a sample: {st.mean(share):.1%}")
