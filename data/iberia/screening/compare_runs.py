"""Compare two runs that differ in unit screening: crown ages, agreement with TimeTree, placements.

Usage: compare_runs.py PROVENANCE_BEFORE PROVENANCE_AFTER UNITS_PLACED_AFTER SCREENED_AFTER [RATE]
Both provenance files need root_to_tip, stem_age, timetree_crown (as written since #36).
"""
import csv
import math
import statistics as st
import sys
from collections import Counter

csv.field_size_limit(10**9)
before_path, after_path, placed_path, screened_path = sys.argv[1:5]
RATE = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0115


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


before = {r["genus"]: r for r in read(before_path) if r["status"] == "subtree"}
after = {r["genus"]: r for r in read(after_path)}
screened = read(screened_path)
aside = Counter(r["genus"] for r in screened if r["decision"] in ("other_genus", "isolated"))

print("Screening decisions:", dict(Counter(r["decision"] for r in screened)))
print("Placement of set-aside units:", dict(Counter(r["status"] for r in read(placed_path) if r["reason"].startswith("screen_"))))
print("Genus status after screening, for genera with a subtree before:",
      dict(Counter(after[g]["status"] if g in after else "absent" for g in before)))


def spearman(x, y):
    def rank(v):
        o = sorted(range(len(v)), key=v.__getitem__)
        r = [0.0] * len(v)
        for i, j in enumerate(o):
            r[j] = i
        return r
    a, b = rank(x), rank(y)
    ma, mb = st.mean(a), st.mean(b)
    return sum((p - ma) * (q - mb) for p, q in zip(a, b)) / math.sqrt(
        sum((p - ma) ** 2 for p in a) * sum((q - mb) ** 2 for q in b))


both = [g for g in before if g in after and after[g]["status"] == "subtree"]
ratio = [float(after[g]["root_to_tip"]) / float(before[g]["root_to_tip"]) for g in both if float(before[g]["root_to_tip"]) > 0]
touched = [g for g in both if aside[g] and float(before[g]["root_to_tip"]) > 0]
tr = [float(after[g]["root_to_tip"]) / float(before[g]["root_to_tip"]) for g in touched]
q = st.quantiles(tr, n=4)
print(f"\nMean root-to-tip distance, genera with a subtree in both runs ({len(both)}): median ratio after/before {st.median(ratio):.3f}")
print(f"   genera that lost units ({len(touched)}): median {st.median(tr):.2f}, quartiles {q[0]:.2f}-{q[2]:.2f}; "
      f"down by more than a quarter {sum(x < .75 for x in tr)}, by more than half {sum(x < .5 for x in tr)}; up by more than a quarter {sum(x > 1.25 for x in tr)}")
cb = [float(before[g]["crown_age"]) for g in both]
ca = [float(after[g]["crown_age"]) for g in both]
print(f"   crown age median {st.median(cb):.2f} -> {st.median(ca):.2f} Ma; changed by more than a quarter in "
      f"{sum(abs(a / b - 1) > .25 for a, b in zip(ca, cb))} genera")


def holdout(label, prov):
    rows = []
    for g in both:
        p = prov[g]
        if not p["timetree_crown"] or float(p["root_to_tip"]) <= 0:
            continue
        tt, stem = float(p["timetree_crown"]), float(p["stem_age"])
        if tt >= stem:
            continue
        low = min(float(p["root_to_tip"]) / RATE, 0.9 * stem)
        rows.append((tt, low, math.sqrt(low * 0.9 * stem)))
    print(f"\n{label}: {len(rows)} genera with a TimeTree crown age below the stem age (TimeTree not used as a bound)")
    for name, k in (("barcode minimum", 1), ("midpoint", 2)):
        lr = [math.log2(r[k] / r[0]) for r in rows]
        print(f"   {name:<16} median estimate/TimeTree {2 ** st.median(lr):.2f}; within a factor 2 {sum(abs(x) <= 1 for x in lr) / len(lr):.0%}; "
              f"within a factor 3 {sum(abs(x) <= math.log2(3) for x in lr) / len(lr):.0%}; rank correlation {spearman([r[k] for r in rows], [r[0] for r in rows]):.2f}")


holdout("Before screening", before)
holdout("After screening", after)
