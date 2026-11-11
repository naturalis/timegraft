"""Summaries of the chronos runs (chronos_bounds.R, chronos_fixed.R).

Usage: summarise.py DIR   (with trees.tsv, A_bounds.tsv, C_fixed.tsv including the newick column, C_big.tsv)
"""
import csv
import math
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

import dendropy

from timegraft.build import tree_length
from timegraft.dating import mpl_ultrametric

csv.field_size_limit(10**9)
out = Path(sys.argv[1])


def read(name):
    with open(out / name) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# A. crown bounded between half and twice the TimeTree age: where does chronos land?
rows = read("A_bounds.tsv")
print(f"A. bounded crown: {len(rows)} runs, {sum(r['ok'] != 'TRUE' for r in rows)} failed")
by = defaultdict(list)
for r in rows:
    if r["ok"] == "TRUE":
        by[(r["model"], r["genus"])].append(r)
for model in ("correlated", "relaxed", "clock"):
    spread, position, dpl, upper, lower, n = [], [], [], 0, 0, 0
    for (m, _), rs in by.items():
        if m != model or len(rs) < 3:
            continue
        ages = [float(r["root"]) for r in rs]
        lo, hi = float(rs[0]["lo"]), float(rs[0]["hi"])
        spread.append(max(ages) / min(ages))
        pl = [float(r["ploglik"]) for r in rs]
        dpl.append(max(pl) - min(pl))
        for a in ages:
            n += 1
            position.append(math.log(a / lo) / math.log(hi / lo))
            upper += a > 0.99 * hi
            lower += a < 1.01 * lo
    q = st.quantiles(spread, n=4)
    print(f"   {model:<11}{len(spread)} genera: oldest/youngest crown over 3 starts, median {st.median(spread):.2f} "
          f"(IQR {q[0]:.2f}-{q[2]:.2f}), above 1.25 in {sum(s > 1.25 for s in spread)}; at the upper bound "
          f"{upper / n:.0%}, at the lower {lower / n:.0%}; median difference in penalised log-likelihood "
          f"between starts {st.median(dpl):.2g}")

# C. crown fixed at 1: within-genus PD (tree length) against mean path length
trees = {r["genus"]: r for r in read("trees.tsv")}
rows = read("C_fixed.tsv")
mpl = {}
for genus in {r["genus"] for r in rows}:
    t = dendropy.Tree.get(data=trees[genus]["newick"], schema="newick", rooting="force-rooted", preserve_underscores=True)
    mpl[genus] = tree_length(mpl_ultrametric(t, 1.0))
by = defaultdict(list)
for r in rows:
    by[(r["model"], r["lambda"])].append(r)
print(f"\nC. crown fixed: {len(mpl)} genera of 4 to 60 units; within-genus PD relative to mean path length")
for (model, lam), rs in by.items():
    ok = [r for r in rs if r["ok"] == "TRUE"]
    ratio = [float(r["len"]) / mpl[r["genus"]] for r in ok]
    q = st.quantiles(ratio, n=20)
    cv = [float(r["rate_cv"]) for r in ok if r["rate_cv"] not in ("NA", "")]
    print(f"   {model:<11}lambda {lam:<4} ok {len(ok)}/{len(rs)}  median {st.median(ratio):.3f} "
          f"(5-95% {q[0]:.2f}-{q[-1]:.2f}); differs by more than 10% in {sum(abs(x - 1) > .1 for x in ratio)}, "
          f"25% in {sum(abs(x - 1) > .25 for x in ratio)}; rate CV median {st.median(cv) if cv else float('nan'):.2f}")
best = defaultdict(int)
phiic = defaultdict(dict)
for r in rows:
    if r["ok"] == "TRUE" and r["lambda"] == "1":
        phiic[r["genus"]][r["model"]] = float(r["phiic"])
for d in phiic.values():
    best[min(d, key=d.get)] += 1
print("   PHIIC prefers:", dict(best))
if (out / "C_big.tsv").exists():
    print("   run time, correlated model:", ", ".join(f"{r['genus']} ({r['n']} units) {float(r['secs']):.0f} s" for r in read("C_big.tsv")))
