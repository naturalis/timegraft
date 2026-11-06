"""Dating a genus-level topology with TimeTree node ages, interpolating the rest.

For every internal node, a few genera from each child clade are sent to
TimeTree, and the age of the node that joins them is taken as a calibration.
Ages that would make a node older than a dated ancestor are discarded. Nodes
without an age are spaced evenly between their nearest dated ancestor and
their oldest dated descendant, as in BLADJ.
"""

from __future__ import annotations

import dendropy

from .timetree import TimeTree

REPORT_COLUMNS = ("node", "n_genera", "age", "source", "example_genera")


def _side(node: dendropy.Node, ncbi: dict[str, int | None], per_side: int) -> list[int]:
    genera = sorted(l.taxon.label for l in node.leaf_iter() if ncbi.get(l.taxon.label))
    # spread the picks over the clade rather than taking neighbours
    if len(genera) > per_side:
        step = len(genera) / per_side
        genera = [genera[int(i * step)] for i in range(per_side)]
    return [ncbi[g] for g in genera]


def date_topology(
    tree: dendropy.Tree,
    ncbi: dict[str, int | None],
    timetree: TimeTree,
    per_side: int = 2,
    root_age: float | None = None,
) -> tuple[dendropy.Tree, list[dict]]:
    """Set branch lengths in Ma on a genus-level topology, in place."""
    age: dict[dendropy.Node, float] = {}
    source: dict[dendropy.Node, str] = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            age[node], source[node] = 0.0, "tip"
            continue
        sides = [_side(child, ncbi, per_side) for child in node.child_node_iter()]
        value = timetree.split_age([s for s in sides if s])
        if value is not None:
            age[node], source[node] = value, "timetree"

    root = tree.seed_node
    if root not in age:
        if root_age is None:
            raise ValueError("TimeTree has no age for the root; set backbone_root_age")
        age[root], source[root] = root_age, "configured"

    # drop calibrations that are not younger than a dated ancestor
    def prune(node: dendropy.Node, ceiling: float) -> None:
        for child in node.child_node_iter():
            if child.is_leaf():
                continue
            next_ceiling = ceiling
            if child in age:
                if age[child] >= ceiling:
                    del age[child]
                    source[child] = "dropped_inconsistent"
                else:
                    next_ceiling = age[child]
            prune(child, next_ceiling)

    prune(root, age[root])

    # oldest dated descendant and the length of the undated chain leading to it
    below: dict[dendropy.Node, tuple[float, int]] = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            below[node] = (0.0, 0)
            continue
        best = (0.0, 0)
        for child in node.child_node_iter():
            if child in age:
                candidate = (age[child], 0)
            else:
                candidate = below[child]
            if candidate[0] > best[0] or (candidate[0] == best[0] and candidate[1] > best[1]):
                best = candidate
        below[node] = (best[0], best[1] + 1)

    for node in tree.preorder_node_iter():
        if node is root or node.is_leaf() or node in age:
            continue
        parent_age = age[node.parent_node]
        floor, steps = below[node]
        age[node] = parent_age - (parent_age - floor) / (steps + 1)
        source[node] = "interpolated" if source.get(node) != "dropped_inconsistent" else "interpolated_after_drop"

    for node in tree.preorder_node_iter():
        if node is root:
            node.edge.length = None
        else:
            node.edge.length = age[node.parent_node] - age[node]

    report = []
    for i, node in enumerate(tree.preorder_node_iter()):
        if node.is_leaf():
            continue
        genera = sorted(l.taxon.label for l in node.leaf_iter())
        report.append(
            {
                "node": f"n{i}",
                "n_genera": len(genera),
                "age": f"{age[node]:.4f}",
                "source": source[node],
                "example_genera": ";".join(genera[:3]),
            }
        )
    return tree, report
