"""Hold-out comparison of crown age estimators against TimeTree crown ages.

Usage: holdout.py DIR PROVENANCE   (DIR with trees.tsv, taxonomy.tsv, outgroup_ages.tsv; the provenance.tsv
of the run before the change, provenance_before.tsv, whose crown_source says which TimeTree ages the earlier
rule accepted)
"""
import csv
import math
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

csv.field_size_limit(10**9)
out, prov_path = Path(sys.argv[1]), sys.argv[2]
RATE = 0.0115


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


trees = {r["genus"]: r for r in read(out / "trees.tsv")}
prov = {r["genus"]: r for r in read(prov_path) if r["status"] == "subtree"}
tax = {r["genus"]: r for r in read(out / "taxonomy.tsv")}
og = defaultdict(list)
for r in read(out / "outgroup_ages.tsv"):
    if r["div_age"]:
        og[r["genus"]].append((float(r["div_age"]), float(r["k2p"])))

# consistency with the archived run: rate-based crown ages should be reproduced
chk = [(float(prov[g]["crown_age"]), float(trees[g]["mrt"]) / RATE) for g in trees if prov[g]["crown_source"] == "rate"]
print("rebuild check: rate-based crown ages reproduced within 0.1%% in %d of %d genera" % (
    sum(abs(a - b) <= 1e-3 * max(a, 1e-9) + 1e-4 for a, b in chk), len(chk)))

# --- saturation: implied pairwise rate against backbone divergence age, over all genus/outgroup-unit pairs
pairs = [(t, d) for g in og for t, d in og[g] if t > 0]
print("\nIngroup-to-outgroup K2P against backbone divergence age (%d pairs)" % len(pairs))
print("age bin (Ma)      pairs  median K2P  implied rate per lineage (K2P / 2T)")
for lo, hi in [(0, 10), (10, 20), (20, 40), (40, 80), (80, 160), (160, 320), (320, 1000)]:
    sel = [(t, d) for t, d in pairs if lo <= t < hi]
    if sel:
        print(f"{lo:>4}-{hi:<5} {len(sel):>10}  {st.median(d for _, d in sel):>9.3f}  {st.median(d / (2 * t) for t, d in sel):.5f}")


# --- implied rate against TimeTree crown age (all TimeTree crown ages below the stem age)
tt = [(float(p["timetree_crown"]), float(trees[g]["mrt"])) for g, p in prov.items()
      if p["timetree_crown"] and g in trees and float(trees[g]["mrt"]) > 0 and float(p["timetree_crown"]) < float(p["stem_age"])]
print("\nMean root-to-tip K2P against TimeTree crown age (%d genera)" % len(tt))
print("crown age (Ma)   genera  median root-to-tip  implied rate per lineage")
for lo, hi in [(0, 5), (5, 10), (10, 20), (20, 40), (40, 80), (80, 1000)]:
    sel = [(a, d) for a, d in tt if lo <= a < hi]
    if sel:
        print(f"{lo:>4}-{hi:<5} {len(sel):>10}  {st.median(d for _, d in sel):>17.3f}  {st.median(d / a for a, d in sel):.4f}")


# --- estimators
def stem_rate(g):
    """Rate from the nearest sister clade: ingroup-to-outgroup K2P over twice the stem age."""
    stem = float(trees[g]["stem_age"])
    near = [d / (2 * t) for t, d in og[g] if abs(t - stem) < 1e-3 * stem]
    return st.median(near) if near else None


def relatives_rate(g, pool):
    """Median rate of calibrated genera in the lowest shared rank (leave-one-out)."""
    for rank in ("family", "order", "class", "phylum"):
        name = tax.get(g, {}).get(rank)
        peers = [r for h, r in pool.items() if h != g and name and tax.get(h, {}).get(rank) == name]
        if len(peers) >= 3:
            return st.median(peers), rank
    peers = [r for h, r in pool.items() if h != g]
    return st.median(peers), "all"


def report(label, truth):
    """truth: genus -> TimeTree crown age."""
    pool = {g: float(trees[g]["mrt"]) / a for g, a in truth.items() if float(trees[g]["mrt"]) > 0}
    est = defaultdict(dict)
    used_rank = defaultdict(int)
    for g, a in truth.items():
        mrt = float(trees[g]["mrt"])
        if mrt <= 0:
            continue
        stem = float(trees[g]["stem_age"])
        est["universal rate (current)"][g] = mrt / RATE
        est["one recalibrated rate (leave-one-out)"][g] = mrt / st.median(v for h, v in pool.items() if h != g)
        r, rank = relatives_rate(g, pool)
        used_rank[rank] += 1
        est["rate from calibrated relatives"][g] = mrt / r
        sr = stem_rate(g)
        if sr:
            est["rate from the stem (sister clade)"][g] = mrt / sr
        est["fraction 0.5 of stem (no barcodes)"][g] = 0.5 * stem
    print(f"\n{label}: {len(pool)} genera; pooled rate median {st.median(pool.values()):.4f}, "
          f"IQR {st.quantiles(pool.values(), n=4)[0]:.4f}-{st.quantiles(pool.values(), n=4)[2]:.4f}")
    print("   relatives' rate taken at rank:", dict(used_rank))
    print(f"   {'estimator':<36}{'n':>5}{'median est/TT':>15}{'median |log2 ratio|':>21}{'within 2x':>11}{'spearman':>10}")
    for name, e in est.items():
        ratios = [math.log2(min(e[g], 0.9 * float(trees[g]["stem_age"])) / truth[g]) for g in e]
        xs, ys = [e[g] for g in e], [truth[g] for g in e]
        print(f"   {name:<36}{len(e):>5}{2 ** st.median(ratios):>15.2f}{st.median(map(abs, ratios)):>21.2f}"
              f"{sum(abs(r) <= 1 for r in ratios) / len(ratios):>11.0%}{spearman(xs, ys):>10.2f}")
    return pool


def spearman(x, y):
    def rank(v):
        order = sorted(range(len(v)), key=v.__getitem__)
        r = [0.0] * len(v)
        for i, j in enumerate(order):
            r[j] = i
        return r
    rx, ry = rank(x), rank(y)
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    return num / math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))


accepted = {g: float(p["timetree_crown"]) for g, p in prov.items() if p["crown_source"].startswith("timetree") and g in trees}
below_stem = {g: float(p["timetree_crown"]) for g, p in prov.items()
              if p["timetree_crown"] and float(p["timetree_crown"]) < float(p["stem_age"]) and g in trees}
big = {g: a for g, a in accepted.items() if int(trees[g]["n"]) >= 4}
pool = report("Accepted TimeTree crown ages (selected by the 3x check, so favours the universal rate)", accepted)
report("Accepted, genera with at least 4 units", big)
report("All TimeTree crown ages below the stem age (includes inflated ones)", below_stem)

# rate by phylum and class among accepted genera
print("\nImplied rate (mean root-to-tip K2P / TimeTree crown age) of accepted genera, by class (n >= 5)")
by = defaultdict(list)
for g, r in pool.items():
    by[(tax.get(g, {}).get("phylum", ""), tax.get(g, {}).get("class", ""))].append(r)
for k, v in sorted(by.items(), key=lambda kv: -len(kv[1])):
    if len(v) >= 5:
        print(f"   {k[0]:<14}{k[1]:<18}{len(v):>4}  median {st.median(v):.4f}  range {min(v):.4f}-{max(v):.4f}")
