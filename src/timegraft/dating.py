"""Dating a rooted genus subtree: the 'crown age' plug-in and rate smoothing.

The crown age C of a genus comes, in order of preference, from a table of
known ages (e.g. TimeTree), from a substitution-rate prior, or from a stated
fraction of the genus stem age S. It is always kept below S.

Smoothing uses mean path length (MPL, Britton et al. 2002): a node's age is
proportional to its mean distance to the tips below it. This is a simple,
dependency-free default; treePL can be plugged in later for larger subtrees.
"""

from __future__ import annotations

from dataclasses import dataclass

import dendropy


@dataclass
class CrownAge:
    age: float
    source: str  # table | rate | rate_capped | fraction | fraction_capped


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
) -> CrownAge:
    """Choose the crown age for one genus.

    ``rate`` is in substitutions per site per million years, per lineage.
    """
    cap = cap_fraction * stem_age
    if genus in known:
        age = known[genus]
        if age < stem_age:
            return CrownAge(age, "table")
        return CrownAge(cap, "table_capped")
    if method == "rate":
        if not rate:
            raise ValueError("crown age method 'rate' needs a rate")
        age = mean_root_to_tip(tree) / rate
        return CrownAge(age, "rate") if age <= cap else CrownAge(cap, "rate_capped")
    if method == "fraction":
        if fraction is None:
            raise ValueError("crown age method 'fraction' needs a fraction")
        age = fraction * stem_age
        return CrownAge(age, "fraction") if age <= cap else CrownAge(cap, "fraction_capped")
    raise ValueError(f"unknown crown age method {method!r}")


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
