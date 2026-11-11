"""Crown age estimators that need no calibrated rate, against TimeTree crown ages (hold-out).

Usage: estimators.py DIR PROVENANCE_BEFORE   (DIR with trees.tsv of the dating benchmark)
"""
import csv
import math
import statistics as st
import sys
from pathlib import Path

csv.field_size_limit(10**9)
out, prov_path = Path(sys.argv[1]), sys.argv[2]
FAST = 0.0115  # per lineage; within the range of published COI rates from young splits


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


trees = {r["genus"]: r for r in read(out / "trees.tsv")}
prov = {r["genus"]: r for r in read(prov_path) if r["status"] == "subtree"}


def yule_fraction(n, q=0.5):
    """Quantile q of crown age over stem age for a pure-birth clade of n tips, birth rate at its ML value ln(n)/S.

    Given the stem age and n, the n-1 node ages are independent with distribution function
    F(t) = (1 - n^(-t/S)) / (1 - 1/n); the crown age is their maximum.
    """
    f = q ** (1.0 / (n - 1))
    return -math.log(1 - (1 - 1 / n) * f) / math.log(n)


rows = []
for g, p in prov.items():
    t = trees.get(g)
    if not t or not p["timetree_crown"] or float(t["mrt"]) <= 0:
        continue
    tt, stem = float(p["timetree_crown"]), float(p["stem_age"])
    if tt >= stem:
        continue
    rows.append({"g": g, "tt": tt, "stem": stem, "n": int(t["n"]), "mrt": float(t["mrt"]),
                 "clean": t["violations"] in ("0", ""), "accepted": p["crown_source"].startswith("timetree")})


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


def ols(X, y):
    """Least squares with intercept, by normal equations (few predictors)."""
    k = len(X[0]) + 1
    A = [[0.0] * k for _ in range(k)]
    b = [0.0] * k
    for xi, yi in zip(X, y):
        v = [1.0] + list(xi)
        for i in range(k):
            b[i] += v[i] * yi
            for j in range(k):
                A[i][j] += v[i] * v[j]
    for i in range(k):  # Gauss-Jordan
        piv = A[i][i]
        A[i] = [a / piv for a in A[i]]
        b[i] /= piv
        for j in range(k):
            if j != i:
                f = A[j][i]
                A[j] = [a - f * c for a, c in zip(A[j], A[i])]
                b[j] -= f * b[i]
    return b


def loo_regression(sel, features):
    est = {}
    for r in sel:
        rest = [o for o in sel if o is not r]
        coef = ols([[math.log(o[f]) for f in features] for o in rest], [math.log(o["tt"]) for o in rest])
        est[r["g"]] = math.exp(coef[0] + sum(c * math.log(r[f]) for c, f in zip(coef[1:], features)))
    return est


def evaluate(label, sel):
    print(f"\n{label}: {len(sel)} genera; TimeTree crown age median {st.median(r['tt'] for r in sel):.1f} Ma, "
          f"crown/stem median {st.median(r['tt'] / r['stem'] for r in sel):.2f}")
    loo_frac = {r["g"]: st.median(o["tt"] / o["stem"] for o in sel if o is not r) for r in sel}
    est = {
        "barcode minimum (0.0115)": {r["g"]: r["mrt"] / FAST for r in sel},
        "half the stem age": {r["g"]: 0.5 * r["stem"] for r in sel},
        "fitted fraction of the stem (LOO)": {r["g"]: loo_frac[r["g"]] * r["stem"] for r in sel},
        "birth model, median": {r["g"]: yule_fraction(r["n"]) * r["stem"] for r in sel},
        "geometric mean of minimum and 0.9 stem": {r["g"]: math.sqrt(min(r["mrt"] / FAST, .9 * r["stem"]) * .9 * r["stem"]) for r in sel},
        "regression on distance (LOO)": loo_regression(sel, ["mrt"]),
        "regression on stem (LOO)": loo_regression(sel, ["stem"]),
        "regression on distance and stem (LOO)": loo_regression(sel, ["mrt", "stem"]),
        "regression on distance, stem and units (LOO)": loo_regression(sel, ["mrt", "stem", "n"]),
    }
    print(f"   {'estimator':<46}{'median est/TT':>14}{'median |log2|':>15}{'within 2x':>11}{'within 3x':>11}{'spearman':>10}")
    for name, e in est.items():
        # every estimate is kept between the barcode minimum and 0.9 of the stem age, as the workflow would
        lr = []
        for r in sel:
            lo, hi = min(r["mrt"] / FAST, .9 * r["stem"]), .9 * r["stem"]
            v = e[r["g"]] if name.startswith("barcode") else min(max(e[r["g"]], lo), hi)
            lr.append(math.log2(v / r["tt"]))
        print(f"   {name:<46}{2 ** st.median(lr):>14.2f}{st.median(map(abs, lr)):>15.2f}"
              f"{sum(abs(x) <= 1 for x in lr) / len(lr):>11.0%}{sum(abs(x) <= math.log2(3) for x in lr) / len(lr):>11.0%}"
              f"{spearman([e[r['g']] for r in sel], [r['tt'] for r in sel]):>10.2f}")
    coef = ols([[math.log(r["mrt"]), math.log(r["stem"]), math.log(r["n"])] for r in sel], [math.log(r["tt"]) for r in sel])
    print("   log TT = %.2f + %.2f log distance + %.2f log stem + %.2f log units" % tuple(coef))


evaluate("All TimeTree crown ages below the stem age", rows)
evaluate("Genus a clade in the barcode tree (clean split)", [r for r in rows if r["clean"]])
evaluate("Clean split and at least 4 units", [r for r in rows if r["clean"] and r["n"] >= 4])
evaluate("Accepted by the earlier 3x rule (favours the barcode minimum)", [r for r in rows if r["accepted"]])
print("\nbirth model fraction by units:", {n: round(yule_fraction(n), 2) for n in (2, 3, 5, 10, 30, 100, 1000)})
print("\nLiterature ages (Ma): estimate at the barcode minimum, birth model median, 0.9 stem")
for g, lit in (("Gammarus", "61"), ("Niphargus", "87"), ("Daphnia", ">=145"), ("Alpheus", "")):
    t = trees[g]
    s, n = float(t["stem_age"]), int(t["n"])
    print(f"   {g:<10} literature {lit:<6} stem {s:6.1f}  units {n:5d}  minimum {float(t['mrt']) / FAST:5.1f}  birth {yule_fraction(n) * s:6.1f}  cap {0.9 * s:6.1f}")
