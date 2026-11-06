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

REPORT_COLUMNS = (
    "node", "parent", "n_genera", "age", "source", "timetree_age", "conflict_ratio", "query_ids", "example_genera",
)


def _ids(value) -> list[int]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [int(v) for v in value if v]
    return [int(value)]


def _side(
    node: dendropy.Node, ncbi: dict[str, object], per_side: int, attempt: int = 0, per_genus: int = 3
) -> list[int]:
    """NCBI ids standing in for one child clade.

    ``ncbi`` maps a genus to its NCBI id, or to a list whose first element is
    the genus id and the rest species ids, since TimeTree knows many species
    but not every genus. Each attempt shifts to other genera of the clade and
    to the next ``per_genus`` species of each genus.
    """
    genera = sorted(l.taxon.label for l in node.leaf_iter() if _ids(ncbi.get(l.taxon.label)))
    if not genera:
        return []
    # spread the picks over the clade rather than taking neighbours
    if len(genera) > per_side:
        step = len(genera) / per_side
        genera = sorted({genera[(int(i * step) + attempt) % len(genera)] for i in range(per_side)})
    out = []
    for genus in genera:
        value = ncbi[genus]
        if isinstance(value, (list, tuple)):
            head = [int(value[0])] if value and value[0] else []
            species = [int(v) for v in value[1:] if v]
        else:
            head, species = _ids(value), []
        window = species[attempt * per_genus : (attempt + 1) * per_genus]
        if attempt and not window and len(genera) <= per_side:
            continue
        out.extend(head + window)
    return out


def _higher_sides(node: dendropy.Node, higher: dict, limit: int) -> list[list[int]]:
    children = list(node.child_node_iter())
    genera = [[l.taxon.label for l in c.leaf_iter()] for c in children]
    lineages = [{uid for g in side for uid, _ in higher.get(g, [])} for side in genera]
    sides = []
    for i, side in enumerate(genera):
        others = set().union(*(lineages[j] for j in range(len(children)) if j != i))
        ids = []
        for genus in side:
            for uid, ncbi in higher.get(genus, []):
                if uid in others:
                    break  # from here up the taxon also holds genera of another side
                if ncbi and ncbi not in ids:
                    ids.append(ncbi)
        if ids:
            sides.append(ids[:limit])
    return sides


def date_topology(
    tree: dendropy.Tree,
    ncbi: dict[str, object],
    timetree: TimeTree,
    per_side: int = 2,
    root_age: float | None = None,
    attempts: int = 4,
    per_genus: int = 3,
    higher: dict[str, list[tuple[int, int | None]]] | None = None,
) -> tuple[dendropy.Tree, list[dict]]:
    """Set branch lengths in Ma on a genus-level topology, in place.

    When TimeTree knows none of the picked taxa on one side of a node, up to
    ``attempts`` different picks are tried. As a last resort, ``higher`` (per
    genus, its enclosing taxa as (OTT id, NCBI id) pairs, nearest first) is
    used: each side is represented by the enclosing taxa of its genera that
    contain no genus of the other sides, so the MRCA stays the same node.
    Such ages are marked ``timetree_higher_taxa``.
    """
    age: dict[dendropy.Node, float] = {}
    source: dict[dendropy.Node, str] = {}
    raw: dict[dendropy.Node, float] = {}  # TimeTree's answer, kept even when it is later dropped
    used: dict[dendropy.Node, list[int]] = {}
    # first picks for every node, fetched concurrently; later attempts are fetched as needed
    if hasattr(timetree, "prefetch"):
        timetree.prefetch(
            [i for s in (_side(c, ncbi, per_side, 0, per_genus) for c in node.child_node_iter()) for i in s]
            for node in tree.postorder_node_iter()
            if not node.is_leaf()
        )
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            age[node], source[node] = 0.0, "tip"
            continue
        for attempt in range(attempts):
            sides = [_side(child, ncbi, per_side, attempt, per_genus) for child in node.child_node_iter()]
            sides = [s for s in sides if s]
            if len(sides) < 2:
                break
            value = timetree.split_age(sides)
            if value is not None:
                age[node], source[node] = value, "timetree"
                raw[node], used[node] = value, sorted({i for side in sides for i in side})
                break
        if node not in age and higher:
            sides = _higher_sides(node, higher, per_side * per_genus)
            value = timetree.split_age(sides)
            if value is not None:
                age[node], source[node] = value, "timetree_higher_taxa"
                raw[node], used[node] = value, sorted({i for side in sides for i in side})

    root = tree.seed_node
    if root not in age:
        if root_age is None:
            raise ValueError("TimeTree has no age for the root; set backbone_build.root_age")
        age[root], source[root] = root_age, "configured"

    # Resolve calibrations that are not younger than their nearest dated
    # ancestor. Equal ages mean TimeTree does not resolve the node: the node
    # is left to interpolation. Otherwise a direct TimeTree age outranks one
    # from enclosing taxa; between equals, the ancestor is kept.
    def dated_ancestor(node: dendropy.Node) -> dendropy.Node | None:
        up = node.parent_node
        while up is not None and up not in age:
            up = up.parent_node
        return up

    conflict: dict[dendropy.Node, float] = {}
    changed = True
    while changed:
        changed = False
        for node in tree.preorder_node_iter():
            if node.is_leaf() or node is root or node not in age:
                continue
            anc = dated_ancestor(node)
            if anc is None or age[node] < age[anc]:
                continue
            if abs(age[node] - age[anc]) <= 1e-9 * max(1.0, age[anc]):
                del age[node]
                source[node] = "timetree_same_as_ancestor"
            elif source[anc] == "timetree_higher_taxa" and source[node] == "timetree" and anc is not root:
                conflict[anc] = age[node] / age[anc]
                del age[anc]
                source[anc] = "dropped_inconsistent"
            else:
                conflict[node] = age[node] / age[anc]
                del age[node]
                source[node] = "dropped_inconsistent"
            changed = True
            break

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
        previous = source.get(node)
        source[node] = {
            "dropped_inconsistent": "interpolated_after_drop",
            "timetree_same_as_ancestor": "interpolated_unresolved_in_timetree",
        }.get(previous, "interpolated")

    for node in tree.preorder_node_iter():
        if node is root:
            node.edge.length = None
        else:
            node.edge.length = age[node.parent_node] - age[node]

    report = []
    names = {node: f"n{i}" for i, node in enumerate(tree.preorder_node_iter())}
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            continue
        genera = sorted(l.taxon.label for l in node.leaf_iter())
        report.append(
            {
                "node": names[node],
                "parent": names.get(node.parent_node, ""),
                "n_genera": len(genera),
                "age": f"{age[node]:.4f}",
                "source": source[node],
                "timetree_age": f"{raw[node]:.4f}" if node in raw else "",
                # TimeTree age over that of the dated ancestor it conflicted with; large values
                # often point to a genus misplaced in the synthetic tree
                "conflict_ratio": f"{conflict[node]:.3f}" if node in conflict else "",
                "query_ids": "+".join(map(str, used.get(node, []))),
                "example_genera": ";".join(genera[:3]),
            }
        )
    return tree, report
