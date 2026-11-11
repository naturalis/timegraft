"""Dating a rooted genus subtree: the crown age and the ages of the nodes below it.

The crown age C of a genus is an age from the barcodes (mean root-to-tip
distance over a substitution rate, or a stated fraction of the stem age S)
kept within bounds from outside sources: never older than a fraction of the
stem age, nor than a TimeTree estimate for the genus, and within a minimum
and maximum from the user's table. A point age in that table is used as it is.
Bounds replace the earlier rule that accepted or rejected a TimeTree age at a
threshold (issue #36).

Below the crown, node ages come from mean path length (MPL, Britton et al.
2002) or LSD2, both of which assume a clock within the genus. Penalised
likelihood (``chronos`` in ape) is available but is not the default: on the
Iberian genera it changed within-genus branch lengths little, and it cannot
estimate the crown age itself (docs/science/results/2026-11-dating.md).
"""

from __future__ import annotations

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
    # table | rate | fraction, followed for the last two by the bound that was
    # applied, if any: _capped (stem age), _timetree_max, _table_max, _table_min
    source: str
    timetree: float | None = None  # the TimeTree estimate for the genus, if any
    barcode: float | None = None  # the age from the barcodes before any bound


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
) -> CrownAge:
    """Choose the crown age for one genus.

    ``known`` is the user's table: per genus a point age, or a mapping with
    ``age``, ``min`` and ``max``. A point age is used as it is (below the stem
    age). Otherwise the age from the barcodes (``rate``, in substitutions per
    site per million years per lineage, or ``fraction`` of the stem age) is
    moved to the nearest bound it violates:

    - at most ``cap_fraction`` of the stem age;
    - at most the TimeTree estimate in ``suggested`` (``timetree_bound`` "max";
      "none" ignores it). That estimate is the MRCA of the species TimeTree
      knows under the genus name: about right for a natural group, far too old
      for a genus that is not one, so it can only lower an age, never raise it;
    - within ``min`` and ``max`` of the table, which then replace TimeTree.

    Where bounds conflict, the stem age wins over the table minimum.
    """
    if timetree_bound not in ("max", "none"):
        raise ValueError(f"unknown timetree_bound {timetree_bound!r}; use 'max' or 'none'")
    cap = cap_fraction * stem_age
    tt = (suggested or {}).get(genus)
    entry = _table_entry(known[genus]) if genus in known else {"age": None, "min": None, "max": None}
    if entry["age"] is not None:
        age = entry["age"]
        return CrownAge(age, "table", tt) if age < stem_age else CrownAge(cap, "table_capped", tt)

    if method == "rate":
        if not rate:
            raise ValueError("crown age method 'rate' needs a rate")
        own = mean_root_to_tip(tree) / rate
    elif method == "fraction":
        if fraction is None:
            raise ValueError("crown age method 'fraction' needs a fraction")
        own = fraction * stem_age
    else:
        raise ValueError(f"unknown crown age method {method!r}")

    upper, upper_source = cap, "capped"
    lower = None
    if entry["min"] is not None or entry["max"] is not None:
        if entry["max"] is not None and entry["max"] < upper:
            upper, upper_source = entry["max"], "table_max"
        if entry["min"] is not None:
            lower = min(entry["min"], upper)
    elif tt is not None and timetree_bound == "max" and tt < upper:
        upper, upper_source = tt, "timetree_max"

    if own > upper:
        return CrownAge(upper, f"{method}_{upper_source}", tt, own)
    if lower is not None and own < lower:
        return CrownAge(lower, f"{method}_table_min" if lower < upper else f"{method}_{upper_source}", tt, own)
    return CrownAge(own, method, tt, own)


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
