"""Precedence rule against the weighted fit, on held-out TimeTree ages (#52).

Usage: python holdout.py <backbone directory> [repeats] [out.tsv]

The backbone directory holds ``backbone.tre`` and ``node_ages.tsv`` of a run
with ``ages.method: fit`` (which reports every calibration with its
uncertainty). For each of ``repeats`` random halves of the calibrations
(the root always kept), the other half is hidden; the backbone is dated from
the visible half by each method, then interpolated as in the workflow
(pure birth on the clade sizes of the run), and the ages of the hidden nodes
are compared with TimeTree's.

Fitting the visible calibrations better says little by itself, since that is
what the fit minimises; this asks whether the ages it gives where it has no
calibration are closer to TimeTree's. TimeTree's own ages are noisy, so
the comparison is also made on the hidden nodes with at least four
estimates, where TimeTree is the better reference.
"""

import csv
import math
import random
import statistics
import sys
from pathlib import Path

import dendropy

from timegraft import agefit
from timegraft.backbone_build import _precedence, interpolate

folder = Path(sys.argv[1])
repeats = int(sys.argv[2]) if len(sys.argv) > 2 else 50
out = sys.argv[3] if len(sys.argv) > 3 else None
tree = dendropy.Tree.get(path=folder / "backbone.tre", schema="newick", rooting="force-rooted",
                         preserve_underscores=True)
rows = {r["node"]: r for r in csv.DictReader(open(folder / "node_ages.tsv"), delimiter="\t")}
names = {f"n{i}": n for i, n in enumerate(tree.preorder_node_iter())}
node_of = {n: k for k, n in names.items()}
root = tree.seed_node
richness = {names[k]: int(r["clade_genera"]) for k, r in rows.items() if r["clade_genera"]}
calibrations, source, n_estimates = {}, {}, {}
for key, r in rows.items():
    if r["calibration_sigma"] or r["source"] == "configured":
        node = names[key]
        calibrations[node] = agefit.Calibration(float(r["timetree_age"] or r["age"]),
                                                float(r["calibration_sigma"] or 1.0),
                                                r["source"].replace("_fit", "").replace("_outlier", ""),
                                                fixed=r["source"] == "configured")
        source[node] = calibrations[node].source
        n_estimates[node] = int(r["timetree_n_estimates"] or 0)
depth = {}
for node in tree.preorder_node_iter():
    depth[node] = 0 if node.parent_node is None else depth[node.parent_node] + 1


def date(visible):
    """Ages of all internal nodes from the visible calibrations, by each method."""
    out = {}
    # precedence
    age = {n: calibrations[n].age for n in visible}
    age.update({n: 0.0 for n in tree.leaf_node_iter()})
    src = {n: source[n] for n in visible}
    src.update({n: "tip" for n in tree.leaf_node_iter()})
    p_age, _, _, _ = _precedence(tree, age, src)
    p_age.update(interpolate(tree, p_age, method="pure_birth", richness=richness))
    out["precedence"] = p_age
    # fit
    parent, nearest = {}, {}
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            continue
        up = node.parent_node
        anc = None if up is None else (up if up in visible else nearest.get(up))
        nearest[node] = anc
        if node in visible:
            parent[node] = anc
    result = agefit.fit_ages(parent, depth, {n: calibrations[n] for n in visible})
    f_age = dict(result.age)
    f_age.update({n: 0.0 for n in tree.leaf_node_iter()})
    f_age.update(interpolate(tree, f_age, method="pure_birth", richness=richness))
    out["fit"] = f_age
    return out


rng = random.Random(52)
candidates = [n for n in calibrations if n is not root]
errors = {m: [] for m in ("precedence", "fit")}
strong = {m: [] for m in ("precedence", "fit")}
records = []
for rep in range(repeats):
    hidden = set(rng.sample(candidates, len(candidates) // 2))
    visible = set(calibrations) - hidden
    dated = date(visible)
    for node in hidden:
        for method, ages in dated.items():
            err = abs(math.log(ages[node] / calibrations[node].age))
            errors[method].append(err)
            if n_estimates[node] >= 4:
                strong[method].append(err)
            records.append({"repeat": rep, "node": node_of[node], "method": method,
                            "timetree_age": f"{calibrations[node].age:.4f}", "age": f"{ages[node]:.4f}",
                            "n_estimates": n_estimates[node]})


def summary(values):
    values = sorted(values)
    factor = [math.exp(v) for v in values]
    return (f"n={len(values)} median factor {statistics.median(factor):.3f}, "
            f"90% {factor[int(0.9 * len(factor))]:.3f}, off by more than 1.5: "
            f"{sum(f > 1.5 for f in factor) / len(factor):.1%}, more than 2: {sum(f > 2 for f in factor) / len(factor):.1%}")


print(f"{len(calibrations)} calibrations, {repeats} random halves hidden")
for method in errors:
    print(f"{method:>10}  all hidden: {summary(errors[method])}")
    print(f"{'':>10}  4+ estimates: {summary(strong[method])}")
if out:
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
