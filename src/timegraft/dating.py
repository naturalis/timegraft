"""Dating a rooted genus subtree: the 'crown age' plug-in and rate smoothing.

The crown age C of a genus comes, in order of preference, from a table of
known ages (e.g. TimeTree), from a substitution-rate prior, or from a stated
fraction of the genus stem age S. It is always kept below S.

Smoothing uses mean path length (MPL, Britton et al. 2002): a node's age is
proportional to its mean distance to the tips below it. This is a simple,
dependency-free default; treePL can be plugged in later for larger subtrees.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import dendropy


@dataclass
class CrownAge:
    age: float
    source: str  # table | timetree | rate | fraction, each possibly *_capped; or *_timetree_rejected
    timetree: float | None = None  # the TimeTree estimate that was considered, if any


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


def crown_age(
    genus: str,
    stem_age: float,
    tree: dendropy.Tree,
    known: dict[str, float],
    method: str = "rate",
    rate: float | None = None,
    fraction: float | None = None,
    cap_fraction: float = 0.9,
    suggested: dict[str, float] | None = None,
    max_ratio: float | None = 3.0,
) -> CrownAge:
    """Choose the crown age for one genus.

    In order: the user's table (``known``); a TimeTree estimate
    (``suggested``) if it passes two checks; otherwise the marker's method.
    ``rate`` is in substitutions per site per million years, per lineage.

    A TimeTree crown age is the MRCA of every species TimeTree knows under the
    genus name, which is far too old when the genus is not monophyletic or
    holds misplaced species. It is therefore rejected when it is not younger
    than the stem age, or, with the rate method, when it exceeds ``max_ratio``
    times the age implied by the barcodes and the rate.
    """
    cap = cap_fraction * stem_age

    def capped(age: float, source: str, tt: float | None = None) -> CrownAge:
        return CrownAge(age, source, tt) if age <= cap else CrownAge(cap, f"{source}_capped", tt)

    if genus in known:
        age = known[genus]
        return CrownAge(age, "table") if age < stem_age else CrownAge(cap, "table_capped")

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

    tt = (suggested or {}).get(genus)
    if tt is not None:
        too_old = tt >= stem_age or (method == "rate" and max_ratio and own > 0 and tt > max_ratio * own)
        if not too_old:
            return capped(tt, "timetree", tt)
        return capped(own, f"{method}_timetree_rejected", tt)
    return capped(own, method)


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
