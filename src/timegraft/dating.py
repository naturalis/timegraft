"""Dating a rooted genus subtree: the crown age and the ages of the nodes below it.

The crown age C of a genus lies between bounds: not younger than the age
from the barcodes (mean root-to-tip distance over a substitution rate), which
saturation makes a minimum, and not older than a fraction of the stem age S or
a TimeTree estimate for the genus, with a minimum and maximum from the user's
table where given. A point age in that table is used as it is. The age within
the bounds is the minimum itself, the midpoint of the range on a log scale, or
the TimeTree age where TimeTree gives one below the stem and the midpoint
otherwise (issue #42).
Bounds replace the earlier rule that accepted or rejected a TimeTree age at a
threshold (issue #36), and both bounds are recorded so that diversity can be
computed at either end.

Below the crown, node ages come from mean path length (MPL, Britton et al.
2002) or LSD2, both of which assume a clock within the genus. Penalised
likelihood (``chronos`` in ape) is available but is not the default: on the
Iberian genera it changed within-genus branch lengths little, and it cannot
estimate the crown age itself (docs/science/results/2026-11-dating.md).
"""

from __future__ import annotations

import math

import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import dendropy


@dataclass
class CrownAge:
    age: float
    # table | rate | fraction | midpoint | timetree; rate and fraction are followed by the
    # bound that was applied, if any: _capped (stem age), _timetree_max,
    # _table_max, _table_min
    source: str
    timetree: float | None = None  # the TimeTree estimate for the genus, if any
    barcode: float | None = None  # the age from the barcodes before any bound
    lower: float | None = None  # the youngest crown age the sources allow (None: no lower bound)
    upper: float | None = None  # the oldest
    upper_source: str = ""  # stem | timetree | table
    # TimeTree crown against the barcodes (#54): plausible | implausible | "" (not tested), the ratio of the
    # TimeTree crown to the age from the barcodes, and the largest ratio rate and saturation can explain
    plausibility: str = ""
    timetree_ratio: float | None = None
    ratio_allowed: float | None = None


def mean_root_to_tip(tree: dendropy.Tree) -> float:
    leaves = tree.leaf_nodes()
    if not leaves:
        return 0.0
    total = 0.0
    for leaf in leaves:
        node, path = leaf, 0.0
        while node.parent_node is not None:
            path += node.edge.length or 0.0
            node = node.parent_node
        total += path
    return total / len(leaves)


def _table_entry(value) -> dict[str, float | None]:
    """A table entry as {age, min, max}; a bare number is a point age."""
    if isinstance(value, Mapping):
        return {k: (None if value.get(k) in (None, "") else float(value[k])) for k in ("age", "min", "max")}
    return {"age": float(value), "min": None, "max": None}


def saturation_factor(distance: float, plateau: float) -> float:
    """How much larger the true distance can be than ``distance`` if distances saturate at ``plateau``.

    With saturation as d = plateau * (1 - exp(-D / plateau)), the true
    distance is D = -plateau * ln(1 - d / plateau); the factor is D / d. It
    is near 1 for small distances, where saturation cannot hide age, and
    unbounded as the distance nears the plateau.
    """
    if distance <= 0:
        return 1.0
    if distance >= plateau:
        return math.inf
    return -plateau * math.log(1 - distance / plateau) / distance


def timetree_plausibility(timetree_age: float, barcode_age: float, distance: float, rate_factor: float = 10.0,
                          plateau: float = 0.15, young_factor: float = 3.0) -> tuple[str, float, float]:
    """Can the barcodes be this young if the crown is as old as TimeTree says? (#54)

    The age from the barcodes (mean root-to-tip distance over the rate) is a
    minimum, too young for two reasons: the genus may evolve more slowly than
    the rate (up to ``rate_factor`` times; COI rates differ by about an
    order of magnitude between sponges or anthozoans and insects), and
    distances saturate (``saturation_factor``). A TimeTree crown older than
    both together allow is ``implausible``. The other way round, the barcode
    age can only be too old if the genus evolves faster than the rate; a
    TimeTree crown more than ``young_factor`` times younger than the barcode
    age is ``implausible_young``: TimeTree knew too few of the species
    queried, which then do not span the crown. Returns the verdict
    (``plausible`` otherwise), the ratio and the largest ratio allowed.
    """
    ratio = timetree_age / barcode_age if barcode_age > 0 else math.inf
    allowed = rate_factor * saturation_factor(distance, plateau)
    if ratio > allowed:
        return "implausible", ratio, allowed
    if young_factor and ratio < 1 / young_factor:
        return "implausible_young", ratio, allowed
    return "plausible", ratio, allowed


def crown_age(
    genus: str,
    stem_age: float,
    tree: dendropy.Tree,
    known: Mapping[str, object],
    method: str = "rate",
    rate: float | None = None,
    fraction: float | None = None,
    cap_fraction: float = 0.9,
    suggested: Mapping[str, float] | None = None,
    timetree_bound: str = "max",
    estimate: str = "minimum",
    plausibility: Mapping | None = None,
) -> CrownAge:
    """Choose the crown age for one genus, and the range the sources allow.

    ``known`` is the user's table: per genus a point age, or a mapping with
    ``age``, ``min`` and ``max``. A point age is used as it is (below the stem
    age). Otherwise the crown age lies between two bounds.

    The upper bound is the smallest of ``cap_fraction`` of the stem age, the
    TimeTree estimate in ``suggested`` (``timetree_bound`` "max"; "none"
    ignores it) and ``max`` of the table, which replaces TimeTree. The TimeTree
    estimate is the MRCA of the species TimeTree knows under the genus name:
    about right for a natural group, far too old for a genus that is not one,
    so it can only lower an age, never raise it.

    The lower bound (method ``rate``) is the age from the barcodes: the mean
    root-to-tip distance over ``rate``, in substitutions per site per million
    years per lineage. Distances saturate, also within genera, so with a rate
    from young calibrations this is a minimum, not an estimate. ``min`` of the
    table raises it. Method ``fraction`` has no lower bound of its own.

    ``estimate`` picks the age within the bounds:

    - ``minimum``: the age from the barcodes (or ``fraction`` of the stem
      age), moved to the nearest bound it violates;
    - ``midpoint`` (method ``rate``): the geometric mean of the two bounds,
      which is the age with the smallest worst-case relative error;
    - ``timetree_midpoint``: the TimeTree age where it is the upper bound
      (below ``cap_fraction`` of the stem age, and no table bounds), even if
      the barcode minimum is older; otherwise ``midpoint`` (method ``rate``)
      or ``minimum`` (method ``fraction``, which has no lower bound).
      TimeTree aggregates published timetrees and is the best estimate there
      is where it covers a genus; the midpoint is the best rate-free estimate
      where it does not.

    Where bounds conflict, the upper bound wins.

    ``plausibility`` (method ``rate`` only; ``enabled``, ``rate_factor``,
    ``saturation``) tests a TimeTree crown against the barcodes
    (:func:`timetree_plausibility`). An implausible one is still an upper
    bound, but with ``timetree_midpoint`` no longer the estimate: the age is
    the midpoint between the barcode age and it (``midpoint_timetree_demoted``).
    One far younger than the barcode age is not used at all (the bounds are
    as without TimeTree; ``plausibility`` ``implausible_young``).
    The holistic analysis found such crowns where the genus is not a natural
    group, or where the global genus is far older than the species sampled.
    """
    if timetree_bound not in ("max", "none"):
        raise ValueError(f"unknown timetree_bound {timetree_bound!r}; use 'max' or 'none'")
    if estimate not in ("minimum", "midpoint", "timetree_midpoint"):
        raise ValueError(f"unknown crown age estimate {estimate!r}; use 'minimum', 'midpoint' or 'timetree_midpoint'")
    cap = cap_fraction * stem_age
    tt = (suggested or {}).get(genus)
    entry = _table_entry(known[genus]) if genus in known else {"age": None, "min": None, "max": None}
    if entry["age"] is not None:
        age, source = (entry["age"], "table") if entry["age"] < stem_age else (cap, "table_capped")
        return CrownAge(age, source, tt, lower=age, upper=age, upper_source="table")

    if method == "rate":
        if not rate:
            raise ValueError("crown age method 'rate' needs a rate")
        own = mean_root_to_tip(tree) / rate
    elif method == "fraction":
        if fraction is None:
            raise ValueError("crown age method 'fraction' needs a fraction")
        if estimate == "midpoint":
            raise ValueError("crown age estimate 'midpoint' needs a lower bound, so method 'rate'")
        own = fraction * stem_age
    else:
        raise ValueError(f"unknown crown age method {method!r}")

    test = {}
    plausibility = dict(plausibility or {})
    if method == "rate" and tt is not None and plausibility.get("enabled", True) and plausibility:
        verdict, ratio, allowed = timetree_plausibility(tt, own, mean_root_to_tip(tree),
                                                        float(plausibility.get("rate_factor", 10.0)),
                                                        float(plausibility.get("saturation", 0.15)),
                                                        float(plausibility.get("young_factor", 3.0)))
        test = dict(plausibility=verdict, timetree_ratio=ratio, ratio_allowed=allowed)
    # a TimeTree crown far younger than the barcodes is no bound at all: its species did not span the crown
    tt_bound = None if test.get("plausibility") == "implausible_young" else tt

    upper, upper_source = cap, "stem"
    table_min = None
    if entry["min"] is not None or entry["max"] is not None:
        if entry["max"] is not None and entry["max"] < upper:
            upper, upper_source = entry["max"], "table"
        table_min = entry["min"]
    elif tt_bound is not None and timetree_bound == "max" and tt_bound < upper:
        upper, upper_source = tt_bound, "timetree"
    suffix = {"stem": "capped", "timetree": "timetree_max", "table": "table_max"}[upper_source]

    lower = own if method == "rate" else None
    if table_min is not None:
        lower = table_min if lower is None else max(lower, table_min)
    if lower is not None:
        lower = min(lower, upper)
    bounds = dict(lower=lower, upper=upper, upper_source=upper_source)

    if estimate == "timetree_midpoint":
        if upper_source == "timetree":
            if test.get("plausibility") == "implausible":
                return CrownAge((lower * upper) ** 0.5, "midpoint_timetree_demoted", tt, own, **bounds, **test)
            return CrownAge(upper, "timetree", tt, own, **bounds, **test)
        estimate = "midpoint" if method == "rate" else "minimum"
    if own > upper:
        return CrownAge(upper, f"{method}_{suffix}", tt, own, **bounds, **test)
    if table_min is not None and own < table_min:
        if table_min >= upper:
            return CrownAge(upper, f"{method}_{suffix}", tt, own, **bounds, **test)
        if estimate == "minimum":
            return CrownAge(table_min, f"{method}_table_min", tt, own, **bounds, **test)
    if estimate == "midpoint":
        return CrownAge((lower * upper) ** 0.5, "midpoint", tt, own, **bounds, **test)
    return CrownAge(own, method, tt, own, **bounds, **test)


def mpl_ultrametric(tree: dendropy.Tree, root_age: float) -> dendropy.Tree:
    """Make ``tree`` ultrametric with its root at ``root_age``, in place."""
    count: dict[dendropy.Node, int] = {}
    total: dict[dendropy.Node, float] = {}
    height: dict[dendropy.Node, int] = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            count[node], total[node], height[node] = 1, 0.0, 0
            continue
        count[node], total[node], height[node] = 0, 0.0, 0
        for child in node.child_node_iter():
            count[node] += count[child]
            total[node] += total[child] + count[child] * (child.edge.length or 0.0)
            height[node] = max(height[node], height[child] + 1)

    root = tree.seed_node
    mpl_root = total[root] / count[root]
    age: dict[dendropy.Node, float] = {root: root_age}
    for node in tree.preorder_node_iter():
        if node is root:
            continue
        parent_age = age[node.parent_node]
        if node.is_leaf():
            age[node] = 0.0
        else:
            if mpl_root > 0:
                guess = root_age * (total[node] / count[node]) / mpl_root
            else:
                guess = root_age * height[node] / height[root]
            # keep ages strictly decreasing towards the tips
            age[node] = min(guess, parent_age * (1.0 - 1e-6))
        node.edge.length = parent_age - age[node]
    root.edge.length = None
    return tree


def lsd2_ultrametric(
    tree: dendropy.Tree, root_age: float, seq_length: int = 650, executable: str = "lsd2"
) -> dendropy.Tree:
    """Date ``tree`` with LSD2 (To et al. 2016): tips at 0, the root fixed at ``root_age``.

    LSD2 fits a strict-clock least-squares model with variances that
    down-weight long branches. Trees with fewer than three tips are dated
    with mean path length, which gives the same answer there.
    """
    leaves = tree.leaf_nodes()
    if len(leaves) < 3:
        return mpl_ultrametric(tree, root_age)
    exe = shutil.which(executable) or (executable if Path(executable).exists() else None)
    if exe is None:
        raise RuntimeError(f"{executable} not found; build it with the build_lsd2 rule or use dating method 'mpl'")
    names = {leaf.taxon.label: f"t{i}" for i, leaf in enumerate(leaves)}
    back = {v: k for k, v in names.items()}
    work = tree.clone(depth=2)
    for leaf in work.leaf_node_iter():
        leaf.taxon.label = names[leaf.taxon.label]
    for edge in work.postorder_edge_iter():
        if edge.length is not None and edge.length < 0:
            edge.length = 0.0
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "in.nwk").write_text(work.as_string(schema="newick", suppress_rooting=True).strip() + "\n")
        (tmp / "dates.txt").write_text(f"{len(names)}\n" + "".join(f"{n} 0\n" for n in names.values()))
        subprocess.run(
            [exe, "-i", "in.nwk", "-d", "dates.txt", "-a", str(-root_age), "-s", str(seq_length), "-o", "out"],
            cwd=tmp, check=True, capture_output=True, text=True,
        )
        nexus = (tmp / "out.date.nexus").read_text()
    line = next(l for l in nexus.splitlines() if l.strip().startswith("tree "))
    newick = re.sub(r"\[&[^\]]*\]", "", line.split("=", 1)[1]).strip()
    dated = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True)
    for leaf in dated.leaf_node_iter():
        leaf.taxon.label = back[leaf.taxon.label]
    dated.seed_node.edge.length = None
    return dated


_CHRONOS_R = """
suppressMessages(library(ape))
a <- commandArgs(TRUE)
t <- read.tree(a[1])
floor <- as.numeric(a[5])
t$edge.length[is.na(t$edge.length) | t$edge.length < floor] <- floor
if (!is.binary(t)) { t <- multi2di(t, random = FALSE); t$edge.length[t$edge.length < floor] <- floor }
set.seed(1)
ch <- suppressWarnings(chronos(t, lambda = as.numeric(a[4]), model = a[3], quiet = TRUE,
                               calibration = makeChronosCalib(t, node = "root", age.min = 1, age.max = 1)))
write.tree(ch, a[2])
"""


def chronos_ultrametric(
    tree: dendropy.Tree,
    root_age: float,
    model: str = "correlated",
    lam: float = 1.0,
    floor: float = 1e-5,
    executable: str = "Rscript",
) -> dendropy.Tree:
    """Date ``tree`` by penalised likelihood with ``chronos`` (ape; Paradis 2013), the root fixed at ``root_age``.

    An option, not the default. With all tips at the present, penalised
    likelihood has no information on the root age (only on the ages of the
    nodes below it), so the root is fixed here as for the other methods.
    Branches shorter than ``floor`` (zero or negative BioNJ lengths) are set to
    it, as ``chronos`` cannot handle them. Trees with fewer than four tips are
    dated with mean path length.
    """
    leaves = tree.leaf_nodes()
    if len(leaves) < 4:
        return mpl_ultrametric(tree, root_age)
    if model not in ("correlated", "relaxed", "discrete", "clock"):
        raise ValueError(f"unknown chronos model {model!r}")
    exe = shutil.which(executable)
    if exe is None:
        raise RuntimeError(f"{executable} not found; install R with the ape package or use dating method 'mpl'")
    names = {leaf.taxon.label: f"t{i}" for i, leaf in enumerate(leaves)}
    back = {v: k for k, v in names.items()}
    work = tree.clone(depth=2)
    for leaf in work.leaf_node_iter():
        leaf.taxon.label = names[leaf.taxon.label]
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "in.nwk").write_text(work.as_string(schema="newick", suppress_rooting=True).strip() + "\n")
        (tmp / "chronos.R").write_text(_CHRONOS_R)
        result = subprocess.run(
            [exe, "chronos.R", "in.nwk", "out.nwk", model, str(lam), str(floor)],
            cwd=tmp, capture_output=True, text=True,
        )
        if result.returncode != 0 or not (tmp / "out.nwk").exists():
            raise RuntimeError(f"chronos failed: {result.stderr.strip()[-500:]}")
        newick = (tmp / "out.nwk").read_text()
    dated = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True)
    for leaf in dated.leaf_node_iter():
        leaf.taxon.label = back[leaf.taxon.label]
    # chronos returns ages relative to a root at 1, and tips within rounding of the present
    height = max(leaf.distance_from_root() for leaf in dated.leaf_node_iter())
    for edge in dated.postorder_edge_iter():
        if edge.length is not None:
            edge.length *= root_age / height
    dated.seed_node.edge.length = None
    for leaf in dated.leaf_node_iter():  # close the rounding gap so that the tree is exactly ultrametric
        leaf.edge.length += root_age - leaf.distance_from_root()
    return dated


def is_ultrametric(tree: dendropy.Tree, tolerance: float = 1e-6) -> bool:
    depths = []
    for leaf in tree.leaf_node_iter():
        node, path = leaf, 0.0
        while node.parent_node is not None:
            path += node.edge.length or 0.0
            node = node.parent_node
        depths.append(path)
    return max(depths) - min(depths) <= tolerance * max(1.0, max(depths))


def root_age(tree: dendropy.Tree) -> float:
    leaf = next(tree.leaf_node_iter())
    node, path = leaf, 0.0
    while node.parent_node is not None:
        path += node.edge.length or 0.0
        node = node.parent_node
    return path
