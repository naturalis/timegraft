"""Weighted global fit of backbone ages (#52)."""

import math
import random

import dendropy
import pytest

from timegraft import agefit
from timegraft.backbone_build import age_fit_summary, date_topology
from timegraft.timetree import TimeTree


def test_calibration_sigma():
    # TimeTree's interval, taken as 95%
    s = agefit.calibration_sigma(ci_low=80.0, ci_high=87.0, floor=0)
    assert s == pytest.approx((math.log(87) - math.log(80)) / 3.92)
    # the spread of estimates: standard error of their mean log
    s = agefit.calibration_sigma([50, 100], floor=0)
    assert s == pytest.approx(math.log(2) / math.sqrt(2) / math.sqrt(2))
    assert agefit.calibration_sigma([50]) == pytest.approx(math.sqrt(0.1 ** 2 + 0.5 ** 2))
    one = agefit.calibration_sigma([], floor=0.1, single=0.5)
    assert agefit.calibration_sigma([], higher_taxa=True) == pytest.approx(2 * one)
    assert agefit.calibration_sigma([], taxonomy_only=True) == pytest.approx(1.5 * one)
    # more studies, smaller uncertainty
    assert agefit.calibration_sigma([90, 100, 110, 95, 105, 100]) < agefit.calibration_sigma([90, 110])


def _brute(parent, values, weights, fixed=()):
    scipy = pytest.importorskip("scipy.optimize")
    import numpy as np

    nodes = list(values)
    index = {n: i for i, n in enumerate(nodes)}
    rows = []
    for n in nodes:
        if parent[n] is not None:
            row = np.zeros(len(nodes))
            row[index[n]], row[index[parent[n]]] = 1, -1
            rows.append(row)
    w = np.array([0.0 if n in fixed else weights[n] for n in nodes])
    y = np.array([values[n] for n in nodes])
    constraints = [scipy.LinearConstraint(np.array(rows), -np.inf, 0)]
    for n in fixed:
        row = np.zeros(len(nodes))
        row[index[n]] = 1
        constraints.append(scipy.LinearConstraint(row[None, :], values[n], values[n]))
    res = scipy.minimize(lambda x: (w * (x - y) ** 2).sum(), y.copy(), jac=lambda x: 2 * w * (x - y),
                         constraints=constraints, method="SLSQP", options={"ftol": 1e-14, "maxiter": 1000})
    return (w * (res.x - y) ** 2).sum()


def test_tree_isotonic_is_exact():
    rng = random.Random(7)
    for _ in range(100):
        n = rng.randint(2, 15)
        parent = {0: None, **{i: rng.randrange(i) for i in range(1, n)}}
        values = {i: rng.uniform(-3, 3) for i in range(n)}
        weights = {i: rng.uniform(0.1, 5) for i in range(n)}
        fit = agefit.tree_isotonic(parent, values, weights)
        assert all(fit[i] <= fit[parent[i]] + 1e-12 for i in range(1, n))
        objective = sum(weights[i] * (fit[i] - values[i]) ** 2 for i in range(n))
        assert objective <= _brute(parent, values, weights) + 1e-6


def test_tree_isotonic_with_ties():
    # a child closing against a tied parent that then merges upward (found by review)
    fit = agefit.tree_isotonic({"r": None, "z": "r", "zz": "z"}, {"r": 0, "z": 100, "zz": 100},
                               {"r": 1, "z": 1, "zz": 1})
    assert fit == pytest.approx({"r": 200 / 3, "z": 200 / 3, "zz": 200 / 3})
    # the same in every order of the names, which only break ties
    for names in (("a", "b", "c"), ("c", "b", "a"), ("b", "c", "a")):
        r, z, zz = names
        fit = agefit.tree_isotonic({r: None, z: r, zz: z}, {r: 0, z: 100, zz: 100}, {r: 1, z: 1, zz: 1})
        assert fit[zz] <= fit[z] + 1e-12 and fit[z] <= fit[r] + 1e-12


def test_tree_isotonic_is_exact_with_ties_and_fixed_nodes():
    rng = random.Random(11)
    for _ in range(150):
        n = rng.randint(2, 12)
        parent = {0: None, **{i: rng.randrange(i) for i in range(1, n)}}
        values = {i: float(rng.randint(0, 4)) for i in range(n)}  # many ties
        weights = {i: float(rng.randint(1, 3)) for i in range(n)}
        fixed = {0} if rng.random() < 0.5 else set()
        if fixed:
            values[0] = 5.0  # a fixed root above everything is always feasible
        fit = agefit.tree_isotonic(parent, values, weights, fixed)
        assert all(fit[i] <= fit[parent[i]] + 1e-12 for i in range(1, n))
        assert all(fit[i] == values[i] for i in fixed)
        objective = sum(weights[i] * (fit[i] - values[i]) ** 2 for i in range(n) if i not in fixed)
        assert objective <= _brute(parent, values, weights, fixed) + 1e-6


def test_fixed_values_hold_exactly():
    # a configured root at 100 against a very precise child at 200 (found by review): the root stays at 100
    parent, depth = {"r": None, "a": "r"}, {"r": 0, "a": 1}
    cal = {"r": agefit.Calibration(100, 1, "configured", fixed=True), "a": agefit.Calibration(200, 1e-9)}
    result = agefit.fit_ages(parent, depth, cal, min_gap=0.01)
    assert result.age["r"] == 100 and result.age["a"] == pytest.approx(100 / 1.01)
    with pytest.raises(ValueError):
        agefit.tree_isotonic({"r": None, "a": "r"}, {"r": 1.0, "a": 2.0}, {"r": 1, "a": 1}, {"r", "a"})


def test_consistent_calibrations_are_kept():
    parent, depth = {"r": None, "a": "r", "b": "a"}, {"r": 0, "a": 1, "b": 3}
    cal = {"r": agefit.Calibration(500, 0.2), "a": agefit.Calibration(300, 0.2), "b": agefit.Calibration(100, 0.2)}
    result = agefit.fit_ages(parent, depth, cal)
    assert result.age == pytest.approx({"r": 500, "a": 300, "b": 100})
    assert not result.outlier and all(abs(r) < 1e-12 for r in result.residual.values())


def test_weak_young_ancestor_gives_way_to_strong_old_descendants():
    # the crab case of the holistic analysis: an ancestor at 54.9 Ma from a query that did not span the
    # node, below descendants at 156.8 (six studies) and 195.5 Ma; the precedence rule kept 54.9
    parent = {"root": None, "crabs": "root", "eubrachyura": "crabs", "heterotremes": "eubrachyura"}
    depth = {"root": 0, "crabs": 1, "eubrachyura": 2, "heterotremes": 3}
    cal = {
        "root": agefit.Calibration(500, 0.15),
        "crabs": agefit.Calibration(54.9, agefit.calibration_sigma([50, 55, 60, 54])),
        "eubrachyura": agefit.Calibration(156.8, agefit.calibration_sigma([150, 160, 155, 158, 152, 165])),
        "heterotremes": agefit.Calibration(195.5, agefit.calibration_sigma([190, 201])),
    }
    result = agefit.fit_ages(parent, depth, cal)
    age = result.age
    assert age["root"] > age["crabs"] > age["eubrachyura"] > age["heterotremes"]  # hard constraint
    assert age["crabs"] > 150  # no longer 54.9
    assert "crabs" in result.outlier  # and reported as the odd one out
    assert abs(math.log(age["eubrachyura"] / 156.8)) < abs(math.log(age["crabs"] / 54.9))


def test_one_wild_calibration_does_not_drag_its_neighbours():
    # a chain of ten consistent nodes and one that is off by a factor of 20 (as Onchidella's Bilateria age was)
    parent = {i: (i - 1 if i else None) for i in range(10)}
    depth = {i: i for i in range(10)}
    truth = {i: 500 * 0.7 ** i for i in range(10)}
    cal = {i: agefit.Calibration(truth[i], 0.2) for i in range(10)}
    cal[6] = agefit.Calibration(truth[6] * 20, 0.2)
    result = agefit.fit_ages(parent, depth, cal)
    assert result.outlier == {6}
    for i in (0, 1, 2, 3, 4, 7, 8, 9):
        assert result.age[i] == pytest.approx(truth[i], rel=0.05)
    # without down-weighting it pulls the three nodes above it up with it
    plain = agefit.fit_ages(parent, depth, cal, huber_k=1e9)
    assert all(plain.age[i] > 1.2 * truth[i] for i in (3, 4, 5))


def test_min_gap_and_fixed_root():
    parent, depth = {"r": None, "a": "r"}, {"r": 0, "a": 2}
    cal = {"r": agefit.Calibration(100, 1, "configured", fixed=True), "a": agefit.Calibration(300, 0.2)}
    result = agefit.fit_ages(parent, depth, cal, min_gap=0.01)
    assert result.age["r"] == pytest.approx(100)
    assert result.age["a"] == pytest.approx(100 / 1.01 ** 2)  # two edges below the root, each at least 1%


def test_misfit():
    cal = {"a": agefit.Calibration(100, 0.1), "b": agefit.Calibration(10, 0.1), "c": agefit.Calibration(5, 0.1)}
    m = agefit.misfit({"a": 100, "b": 20, "c": 5}, cal)
    assert m["n"] == 3 and m["severe"] == 1 and m["max_factor"] == pytest.approx(2.0)
    assert m["median_factor"] == pytest.approx(1.0)


def _fake(ages, estimates=None):
    def fetch(url, payload=None):
        ids = [int(i) for i in url.split("/mrca/id/")[1].split("/")[0].split("+")]
        age = ages.get(frozenset(ids))
        est = (estimates or {}).get(frozenset(ids), [age] if age else [])
        return {"found_ids": ids, "precomputed_age": age,
                "time_estimates": "{" + ",".join(map(str, est)) + "}" if est else None}

    return fetch


def test_date_topology_fit_against_precedence():
    # (((A,B),C),D): (A,B,C) at 120 is older than the root at 100. The precedence rule keeps the root,
    # whatever the evidence; the fit weighs them, and here (A,B,C) has six studies and the root one.
    newick = "(((A,B),C),D);"
    ncbi = {"A": 1, "B": 2, "C": 3, "D": 4}
    ages = {frozenset({1, 2, 4}): 100.0, frozenset({1, 2, 3}): 120.0, frozenset({1, 2}): 30.0}
    estimates = {frozenset({1, 2, 3}): [110, 115, 120, 125, 118, 122]}
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted")
    _, report = date_topology(tree, ncbi, TimeTree(fetch=_fake(ages, estimates)), root_age=None)
    by = {r["example_genera"]: r for r in report}
    root, abc, ab = report[0], by["A;B;C"], by["A;B"]
    assert float(root["age"]) > float(abc["age"]) > float(ab["age"]) > 0  # monotonic
    assert abc["source"] == "timetree_fit" and abc["precedence_source"] == "interpolated_after_drop"
    assert root["source"] == "timetree_fit_outlier" or float(root["age"]) > 110
    assert float(abc["age"]) > 110  # close to its own, well-supported age
    assert abc["calibration_sigma"] and abc["calibration_z"] and abc["fitted_age"] == abc["age"]
    assert ab["source"] == "timetree_fit" and float(ab["age"]) == pytest.approx(30, rel=1e-6)
    summary = age_fit_summary(report)
    assert summary["final"]["n"] == summary["precedence"]["n"] == 3
    assert summary["final"]["weighted_rms"] < summary["precedence"]["weighted_rms"]
    # every edge has a positive length
    assert all(e.length is None or e.length > 0 for e in tree.preorder_edge_iter() if e.head_node.parent_node)


def test_conflict_metrics_use_the_same_weights_in_either_method():
    newick = "(((A,B),C),D);"
    ncbi = {"A": 1, "B": 2, "C": 3, "D": 4}
    ages = {frozenset({1, 2, 4}): 100.0, frozenset({1, 2, 3}): 120.0, frozenset({1, 2}): 30.0}
    estimates = {frozenset({1, 2, 3}): [110, 115, 120, 125, 118, 122]}
    reports = {}
    for method in ("fit", "precedence"):
        tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted")
        _, reports[method] = date_topology(tree, ncbi, TimeTree(fetch=_fake(ages, estimates)), method=method)
    sigma = {m: [r["calibration_sigma"] for r in rep] for m, rep in reports.items()}
    assert sigma["fit"] == sigma["precedence"] and any(sigma["fit"])
    assert len(set(s for s in sigma["fit"] if s)) > 1  # not one default for all
    assert age_fit_summary(reports["fit"])["precedence"] == age_fit_summary(reports["precedence"])["precedence"]
    assert all(r["fitted_age"] == "" for r in reports["precedence"])


def test_date_topology_unknown_method():
    tree = dendropy.Tree.get(data="(A,B);", schema="newick", rooting="force-rooted")
    with pytest.raises(ValueError):
        date_topology(tree, {"A": 1, "B": 2}, TimeTree(fetch=_fake({frozenset({1, 2}): 5.0})), method="vote")
