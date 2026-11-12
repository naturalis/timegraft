"""Interpolated backbone node ages against held-out TimeTree ages (#16).

Usage: python holdout.py <run directory with backbone/> <richness.tsv> [out.tsv]

Reads the dated backbone and its node report, keeps the ages that came from
TimeTree directly, hides some of them and interpolates them again with each
method, then compares the result to the hidden TimeTree age. Two designs:
leave one out (each dated node hidden on its own, its neighbours known) and
repeated random halves (half of the dated nodes hidden at once, which gives
longer undated chains, closer to the poorly covered parts of the tree).

Methods: ``even`` (BLADJ), and Yule heights from the genera in the backbone
(``pb_sampled``), from the genera in the full clade of the synthetic tree
(``pb_genera``) and from its tips, mostly species (``pb_species``); the last
two come from ``richness.py``.
"""

import csv
import math
import random
import statistics
import sys
from pathlib import Path

import dendropy

from timegraft.backbone_build import interpolate

run = Path(sys.argv[1])
rich_path = Path(sys.argv[2])
out_path = Path(sys.argv[3]) if len(sys.argv) > 3 else None
tree = dendropy.Tree.get(path=run / "backbone" / "backbone.tre", schema="newick", rooting="force-rooted",
                         preserve_underscores=True)
rows = {r["node"]: r for r in csv.DictReader(open(run / "backbone" / "node_ages.tsv"), delimiter="\t")}
names = {f"n{i}": node for i, node in enumerate(tree.preorder_node_iter())}
# the report names nodes in the same preorder; check that it lines up
for name, row in rows.items():
    node = names[name]
    assert len(node.leaf_nodes()) == int(row["n_genera"]), name
root = tree.seed_node
rich_rows = list(csv.DictReader(open(rich_path), delimiter="\t"))
synth_genera = {names[r["node"]]: int(r["synth_genera"]) for r in rich_rows if r["synth_genera"]}
synth_tips = {names[r["node"]]: int(r["synth_tips"]) for r in rich_rows if r["synth_tips"]}
METHODS = {
    "even": dict(method="even"),
    "pb_sampled": dict(method="pure_birth"),
    "pb_genera": dict(method="pure_birth", richness=synth_genera),
    "pb_species": dict(method="pure_birth", richness=synth_tips),
}
anchors = {names[n]: float(r["age"]) for n, r in rows.items() if r["source"] in ("timetree", "timetree_higher_taxa")}
anchors[root] = float(rows["n0"]["age"])
truth = [names[n] for n, r in rows.items() if r["source"] == "timetree" and names[n] is not root]


def predict(hidden):
    known = {n: a for n, a in anchors.items() if n not in hidden}
    return {m: interpolate(tree, known, **kw) for m, kw in METHODS.items()}


records = []
for node in truth:
    pred = predict({node})
    records.append(("leave_one_out", 0, node, {m: pred[m][node] for m in METHODS}))
rng = random.Random(16)
for rep in range(50):
    hidden = set(rng.sample(truth, len(truth) // 2))
    pred = predict(hidden)
    for node in hidden:
        records.append(("random_half", rep, node, {m: pred[m][node] for m in METHODS}))

rev = {v: k for k, v in names.items()}
table = []
for design, rep, node, pred in records:
    t = anchors[node]
    n = len(node.leaf_nodes())
    row = {"design": design, "rep": rep, "node": rev[node], "n_genera": n, "timetree_age": f"{t:.3f}"}
    for m, p in pred.items():
        row[m] = f"{p:.3f}"
    table.append(row)
if out_path:
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0]), delimiter="\t")
        w.writeheader()
        w.writerows(table)


def summary(sel, label):
    if not sel:
        return
    parts = [f"{label:<22} n={len(sel):>5}"]
    for m in METHODS:
        err = [math.log(float(r[m]) / float(r["timetree_age"])) for r in sel]
        within2 = sum(abs(e) <= math.log(2) for e in err) / len(err)
        parts.append(f"{m:>10}: bias {statistics.median(err) / math.log(2):+.2f} "
                     f"|err| {statistics.median(abs(e) for e in err) / math.log(2):.2f} "
                     f"<2x {within2:.0%}")
    print("  ".join(parts))


for design in ("leave_one_out", "random_half"):
    sel = [r for r in table if r["design"] == design]
    print(design, "(bias: median log2 of interpolated/TimeTree age; |err|: median absolute log2; <2x: share within a factor 2)")
    summary(sel, "  all")
    for lo, hi in ((2, 2), (3, 5), (6, 20), (21, 100), (101, 10**6)):
        summary([r for r in sel if lo <= r["n_genera"] <= hi], f"  {lo}-{hi if hi < 10**6 else ''} genera")
    for lo, hi in ((0, 50), (50, 150), (150, 400), (400, 10**4)):
        summary([r for r in sel if lo <= float(r["timetree_age"]) < hi], f"  TimeTree age {lo}-{hi if hi < 10**4 else ''} Ma")
