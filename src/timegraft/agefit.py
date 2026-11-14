"""Backbone node ages from all TimeTree calibrations at once (#52).

Before #52 a calibration that was not younger than its nearest dated
ancestor was dropped, by position in the tree: the ancestor won, however
weak its evidence (one estimate, or a query that did not span the node)
and however strong the descendant's. Here every calibration is a soft
constraint instead, weighted by its evidence, and all backbone ages are
fitted together under the hard constraint that ages decrease towards the
tips.

The fit is a weighted least-squares fit of log ages, so that a factor two
counts the same at 5 and at 500 Ma, made robust with Huber weights: a
calibration far out of line with the others (in units of its own
uncertainty) is down-weighted rather than allowed to drag its neighbours,
and is reported as an outlier. Nodes without a calibration are interpolated
between the fitted ages afterwards, as before.

The monotone fit is an isotonic regression on the tree of calibrated nodes,
solved exactly (:func:`tree_isotonic`). A minimum gap per edge keeps every
fitted age strictly younger than that of its calibrated ancestor, so no
branch has length zero.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Hashable, Mapping, Sequence

# Defaults for the uncertainty of a calibration, as a standard deviation of its log age.
SIGMA_FLOOR = 0.1  # even a well-studied node is not known to better than about 10%
SIGMA_SINGLE = 0.5  # a single study estimate: about a factor 1.65 either way
HIGHER_TAXA_FACTOR = 2.0  # ages from enclosing taxa date the node only indirectly
TAXONOMY_ONLY_FACTOR = 1.5  # nodes no source tree supports may not exist as dated
HUBER_K = 1.5
OUTLIER_Z = 3.0
MIN_GAP = 0.01  # each edge between calibrated nodes is at least 1% of the older age


@dataclass
class Calibration:
    age: float  # TimeTree's age, Ma
    sigma: float  # standard deviation of the log age
    source: str = "timetree"
    fixed: bool = False  # a configured root age: not fitted

    @property
    def weight(self) -> float:
        return 1.0 / (self.sigma * self.sigma)


def calibration_sigma(estimates: Sequence[float] = (), ci_low: float | None = None, ci_high: float | None = None,
                      higher_taxa: bool = False, taxonomy_only: bool = False, floor: float = SIGMA_FLOOR,
                      single: float = SIGMA_SINGLE, higher_factor: float = HIGHER_TAXA_FACTOR,
                      taxonomy_factor: float = TAXONOMY_ONLY_FACTOR) -> float:
    """Uncertainty of a TimeTree age as a standard deviation of its log, from the evidence behind it.

    From TimeTree's confidence interval where it gives one (taken as 95%),
    else from the spread of the study estimates (the standard error of their
    mean log), else ``single``; combined with ``floor`` and inflated for ages
    from enclosing taxa and for nodes that rest on the taxonomy alone.
    """
    estimates = [e for e in estimates if e and e > 0]
    if ci_low and ci_high and ci_high > ci_low > 0:
        spread = (math.log(ci_high) - math.log(ci_low)) / (2 * 1.96)
    elif len(estimates) >= 2:
        logs = [math.log(e) for e in estimates]
        mean = sum(logs) / len(logs)
        sd = math.sqrt(sum((x - mean) ** 2 for x in logs) / (len(logs) - 1))
        spread = sd / math.sqrt(len(logs))
    else:
        spread = single
    sigma = math.sqrt(floor * floor + spread * spread)
    if higher_taxa:
        sigma *= higher_factor
    if taxonomy_only:
        sigma *= taxonomy_factor
    return sigma


def tree_isotonic(parent: Mapping[Hashable, Hashable | None], values: Mapping[Hashable, float],
                  weights: Mapping[Hashable, float], fixed: set | frozenset = frozenset()) -> dict[Hashable, float]:
    """Weighted least-squares fit with every node at most its parent (the root is the maximum).

    ``parent`` maps each node to its parent (None for the root); every node
    needs a value and a positive weight. Nodes in ``fixed`` keep their value
    exactly. Exact, by Thompson's block-merging algorithm for tree orders:
    repeatedly take the open block with the largest value; if it is larger
    than the block holding its top node's parent, merge it into that block,
    otherwise close it.

    Why closed blocks stay valid: a block's value only goes down when it is
    merged upward, which happens only to the open block with the largest
    value; a block absorbing a child only goes up. So when a block is closed,
    the block above it is either closed already, the root's, or open with the
    same value, and of equal values the block nearest the root is taken first
    (merged upward or closed) before any block below it. Without that order,
    a tie lets a child close against a parent that then merges upward and
    drops below it.
    """
    nodes = list(values)
    depth: dict = {}

    def node_depth(n):
        chain = []
        while n not in depth and n is not None:
            chain.append(n)
            n = parent.get(n)
        d = -1 if n is None else depth[n]
        for m in reversed(chain):
            d += 1
            depth[m] = d
        return depth

    for n in nodes:
        node_depth(n)
    block_of = {n: n for n in nodes}  # union-find: representative
    total = {n: weights[n] * values[n] for n in nodes}
    mass = {n: float(weights[n]) for n in nodes}
    pinned = {n: (float(values[n]) if n in fixed else None) for n in nodes}  # a fixed value in the block
    top = {n: n for n in nodes}  # the node of the block nearest the root

    def find(n):
        while block_of[n] != n:
            block_of[n] = block_of[block_of[n]]
            n = block_of[n]
        return n

    def value(b):
        return pinned[b] if pinned[b] is not None else total[b] / mass[b]

    open_blocks = {n for n in nodes if parent.get(n) is not None}
    while open_blocks:
        b = max(open_blocks, key=lambda x: (value(x), -depth[top[x]], str(x)))
        open_blocks.discard(b)
        p = find(parent[top[b]])
        if value(b) <= value(p):
            continue
        if pinned[b] is not None and pinned[p] is not None:
            raise ValueError(f"fixed values conflict: {top[b]!r} ({pinned[b]}) is older than its ancestor "
                             f"{top[p]!r} ({pinned[p]})")
        # merge b into p; the merged block keeps p's top, and a fixed value of either
        block_of[b] = p
        total[p] += total[b]
        mass[p] += mass[b]
        if pinned[p] is None:
            pinned[p] = pinned[b]
    return {n: value(find(n)) for n in nodes}


@dataclass
class FitResult:
    age: dict = field(default_factory=dict)  # fitted age per calibrated node
    residual: dict = field(default_factory=dict)  # log(fitted / TimeTree)
    z: dict = field(default_factory=dict)  # residual over sigma
    weight: dict = field(default_factory=dict)  # final weight, after down-weighting
    outlier: set = field(default_factory=set)
    iterations: int = 0


def fit_ages(parent: Mapping[Hashable, Hashable | None], depth: Mapping[Hashable, int],
             calibrations: Mapping[Hashable, Calibration], min_gap: float = MIN_GAP, huber_k: float = HUBER_K,
             outlier_z: float = OUTLIER_Z, max_iterations: int = 50) -> FitResult:
    """Fit log ages to the calibrations of ``parent``'s nodes under monotonicity.

    ``parent`` is the tree of calibrated nodes (each node's nearest
    calibrated ancestor; None for the root) and ``depth`` the number of
    edges from the root of the full tree to each, used for the minimum gap.
    A fixed calibration keeps its age exactly; calibrations that would be
    older than it are fitted at its age less the gap.
    """
    gap = math.log1p(min_gap)
    # with z = log(age) + gap * depth, "age at most the ancestor's, less the gap per edge" is z <= z(parent)
    target = {n: math.log(c.age) + gap * depth[n] for n, c in calibrations.items()}
    fixed = {n for n, c in calibrations.items() if c.fixed}  # held exactly, not by a large weight
    base = {n: (1.0 if c.fixed else c.weight) for n, c in calibrations.items()}
    weights = dict(base)
    result = FitResult()
    previous = None
    for iteration in range(1, max_iterations + 1):
        z = tree_isotonic(parent, target, weights, fixed)
        result.iterations = iteration
        if previous is not None and max(abs(z[n] - previous[n]) for n in z) < 1e-9:
            break
        previous = z
        for n, c in calibrations.items():
            if c.fixed:
                continue
            standardised = abs(z[n] - target[n]) / c.sigma
            weights[n] = base[n] * (1.0 if standardised <= huber_k else huber_k / standardised)
    for n, c in calibrations.items():
        log_age = z[n] - gap * depth[n]
        result.age[n] = c.age if c.fixed else math.exp(log_age)  # exact, not through exp(log())
        result.residual[n] = log_age - math.log(c.age)
        result.z[n] = result.residual[n] / c.sigma
        result.weight[n] = weights[n] if not c.fixed else float("inf")
        if not c.fixed and abs(result.z[n]) > outlier_z:
            result.outlier.add(n)
    return result


def misfit(final: Mapping[Hashable, float], calibrations: Mapping[Hashable, Calibration],
           severe: float = 1.5) -> dict[str, float]:
    """How far the final ages are from TimeTree's, over all calibrations: the conflict metrics of #52.

    ``median_factor`` and ``p90_factor`` are factors (final over TimeTree or
    the reverse, whichever is larger); ``severe`` counts calibrations off by
    more than the factor ``severe``; ``weighted_rms`` is the root mean square
    of log residuals over sigma.
    """
    factors, zs = [], []
    for n, c in calibrations.items():
        if c.fixed or n not in final:
            continue
        r = math.log(final[n] / c.age)
        factors.append(math.exp(abs(r)))
        zs.append(r / c.sigma)
    if not factors:
        return {"n": 0}
    factors.sort()

    def quantile(q):
        return factors[min(len(factors) - 1, int(q * len(factors)))]

    return {
        "n": len(factors),
        "median_factor": round(quantile(0.5), 4),
        "p90_factor": round(quantile(0.9), 4),
        "max_factor": round(factors[-1], 4),
        "severe": sum(1 for f in factors if f > severe),
        "weighted_rms": round(math.sqrt(sum(z * z for z in zs) / len(zs)), 4),
    }
