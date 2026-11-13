"""Dating a genus-level topology with TimeTree node ages, interpolating the rest.

For every internal node, a few genera from each child clade are sent to
TimeTree, and the age of the node that joins them is taken as a calibration.
Conflicting calibrations are resolved by a weighted fit of all of them under
the constraint that ages decrease towards the tips (#52, :mod:`timegraft.agefit`),
or, with the precedence rule used before, by discarding ages that would make
a node older than a dated ancestor. Nodes
without an age are placed between their nearest dated ancestor and their
oldest dated descendant, by default by their expected height under a
pure-birth model, or evenly as in BLADJ (#16).
"""

from __future__ import annotations

import logging
import math

import dendropy

from .col import GENUS_COLUMNS as COL_GENUS_COLUMNS
from .col import compare_ids
from . import agefit
from . import requery as rq
from .requery import COLUMNS as REQUERY_COLUMNS
from .timetree import EVIDENCE_COLUMNS, TimeTree

REPORT_COLUMNS = (
    "node", "parent", "n_genera", "clade_genera", "age", "source", "timetree_age", "conflict_ratio", "query_ids",
    "example_genera",
    # evidence and provenance for the consistency layer (#46, #47); none of these changes an age
    "timetree_attempt", "conflict_with",
    *(f"timetree_{c}" for c in EVIDENCE_COLUMNS),
    "synth_node", *(f"synth_{k}" for k in ("supported_by", "conflicts_with", "resolves", "partial_path_of",
                                           "terminal")),
    # weighted fit of all calibrations (#52): the calibration's uncertainty and weight, the fitted age, the
    # residual log(fitted / TimeTree) and over sigma, whether it is an outlier; and, for comparison, the age
    # and source the precedence rule before #52 gives the node
    "calibration_sigma", "calibration_weight", "fitted_age", "calibration_residual", "calibration_z",
    "calibration_outlier", "precedence_age", "precedence_source",
    # re-query of conflicting nodes with disjoint picks, and what it says about the cause (#53)
    *REQUERY_COLUMNS,
)
# columns reserved for later issues, written empty until those land
RESERVED_COLUMNS = ()
AGE_METHODS = ("fit", "precedence")


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


INTERPOLATION_METHODS = ("pure_birth", "even")


def yule_height(n: int) -> float:
    """Expected crown age of a pure-birth (Yule) clade of ``n`` tips, in units of 1/lambda.

    With k lineages the waiting time to the next split has mean 1/(k lambda),
    so the crown of n tips is expected sum_{k=2}^{n} 1/k = H_n - 1 time
    units before the present; a single tip has height 0.
    """
    if n <= 1:
        return 0.0
    if n <= 1000:
        return sum(1.0 / k for k in range(2, n + 1))
    # asymptotic expansion of the harmonic number, exact to well below 1e-12 here
    return math.log(n) + 0.5772156649015329 + 1 / (2 * n) - 1 / (12 * n * n) - 1.0


def interpolate(
    tree: dendropy.Tree,
    age: dict[dendropy.Node, float],
    method: str = "pure_birth",
    richness: dict[dendropy.Node, int] | None = None,
) -> dict[dendropy.Node, float]:
    """Ages for the internal nodes missing from ``age``; the root must be dated.

    Each undated node is placed, from the root down, between the age of its
    parent and that of its oldest dated descendant (a tip, at 0, if there is
    none), so that ages decrease strictly towards the tips and every dated
    age is kept.

    ``even`` spaces the nodes of an undated chain evenly, as BLADJ does.
    ``pure_birth`` places them by their expected height under a pure-birth
    model (:func:`yule_height`): a node of richness ``n`` under a parent of
    richness ``n_p`` and above a dated floor of richness ``n_f`` is put at the
    fraction (h(n) - h(n_f)) / (h(n_p) - h(n_f)) of the way from the floor to
    the parent, with h = 0 for a tip. This amounts to a constant speciation
    rate between each pair of dated nodes, fitted to that pair.

    ``richness`` is the number of species (or other terminal taxa) in the
    full clade of each internal node, for example from the synthetic tree.
    Without it the genera in ``tree`` are counted, which with sparse sampling
    makes shallow nodes far too young (see ``data/iberia/backbone-interpolation``).
    Where the heights do not increase strictly towards the root, the node is
    spaced evenly instead.
    """
    if method not in INTERPOLATION_METHODS:
        raise ValueError(f"unknown interpolation method {method!r}; use one of {', '.join(INTERPOLATION_METHODS)}")
    root = tree.seed_node
    if root not in age:
        raise ValueError("the root needs an age before the other nodes can be interpolated")
    tips: dict[dendropy.Node, int] = {}
    for node in tree.postorder_node_iter():
        tips[node] = 1 if node.is_leaf() else sum(tips[c] for c in node.child_node_iter())

    def height(node: dendropy.Node) -> float:
        if node.is_leaf():
            return 0.0
        n = richness.get(node) if richness else None
        return yule_height(n if n else tips[node])

    heights = {node: height(node) for node in tree.preorder_node_iter()} if method == "pure_birth" else {}
    # oldest dated descendant: its age, the length of the undated chain to it, the node itself
    below: dict[dendropy.Node, tuple[float, int, dendropy.Node]] = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            below[node] = (0.0, 0, node)
            continue
        best = None
        for child in node.child_node_iter():
            candidate = (age[child], 0, child) if child in age else below[child]
            if best is None or candidate[0] > best[0] or (candidate[0] == best[0] and candidate[1] > best[1]):
                best = candidate
        below[node] = (best[0], best[1] + 1, best[2])

    out: dict[dendropy.Node, float] = {}
    for node in tree.preorder_node_iter():
        if node is root or node.is_leaf() or node in age:
            continue
        parent = node.parent_node
        parent_age = age[parent] if parent in age else out[parent]
        floor, steps, floor_node = below[node]
        fraction = 1.0 / (steps + 1)
        if method == "pure_birth":
            h_parent, h_node, h_floor = heights[parent], heights[node], heights[floor_node]
            if h_parent > h_node > h_floor:
                fraction = (h_parent - h_node) / (h_parent - h_floor)
        out[node] = parent_age - (parent_age - floor) * fraction
    return out


def date_topology(
    tree: dendropy.Tree,
    ncbi: dict[str, object],
    timetree: TimeTree,
    per_side: int = 2,
    root_age: float | None = None,
    attempts: int = 4,
    per_genus: int = 3,
    higher: dict[str, list[tuple[int, int | None]]] | None = None,
    interpolation: str = "pure_birth",
    richness: dict[dendropy.Node, int] | None = None,
    synth_nodes: dict[dendropy.Node, str] | None = None,
    support: dict[str, dict[str, int]] | None = None,
    method: str = "fit",
    fit: dict | None = None,
    requery: dict | None = None,
) -> tuple[dendropy.Tree, list[dict]]:
    """Set branch lengths in Ma on a genus-level topology, in place.

    Nodes left without an age are filled in by :func:`interpolate`, with
    ``interpolation`` as its method. ``pure_birth`` needs ``richness``, the
    number of genera in the full clade of each internal node (see
    :func:`timegraft.otol_static.clade_richness`); without it the nodes are
    spaced evenly, since a pure-birth model of the sampled genera alone makes
    shallow nodes far too young (#16).

    When TimeTree knows none of the picked taxa on one side of a node, up to
    ``attempts`` different picks are tried. As a last resort, ``higher`` (per
    genus, its enclosing taxa as (OTT id, NCBI id) pairs, nearest first) is
    used: each side is represented by the enclosing taxa of its genera that
    contain no genus of the other sides, so the MRCA stays the same node.
    Such ages are marked ``timetree_higher_taxa``.

    The report has one row per internal node (see ``REPORT_COLUMNS``). Besides
    the age and its source it records the evidence behind it, none of which
    changes an age (#47): which attempt gave TimeTree's answer, the node it
    conflicted with, what TimeTree said (see :meth:`TimeTree.evidence`), and,
    given ``synth_nodes`` (the label of the synthesis node behind each node)
    and ``support`` (per label, the number of source trees per kind of
    annotation, see :func:`timegraft.otol_static.synthesis_support`), how
    well the synthetic tree supports the node.

    How conflicting calibrations are resolved depends on ``method``. With
    ``fit`` (#52) all TimeTree ages are fitted together as soft constraints,
    weighted by their evidence, under the hard constraint that ages decrease
    towards the tips (:mod:`timegraft.agefit`; ``fit`` holds its settings:
    ``sigma_floor``, ``sigma_single``, ``higher_taxa_factor``,
    ``taxonomy_only_factor``, ``huber_k``, ``outlier_z``, ``min_gap``). With
    ``precedence`` an age that is not younger than its nearest dated
    ancestor is dropped, as before #52. Either way the report also gives the
    ages and sources of the precedence rule, for comparison.

    ``requery`` (#53; settings ``mode``: ``classify``, ``apply`` or ``off``;
    ``z_trigger``, ``ratio_trigger``, ``move``) re-asks TimeTree about
    conflicting nodes with picks disjoint from the first query, and
    classifies the conflict (:mod:`timegraft.requery`). With ``apply`` (fit
    only) the re-query age replaces the first for nodes classed as
    ``under_spanning`` or ``pick_specific``, and the ages are fitted again.
    """
    if method not in AGE_METHODS:
        raise ValueError(f"unknown age method {method!r}; use one of {', '.join(AGE_METHODS)}")
    fit = dict(fit or {})
    age: dict[dendropy.Node, float] = {}
    source: dict[dendropy.Node, str] = {}
    raw: dict[dendropy.Node, float] = {}  # TimeTree's answer, kept even when it is later dropped
    used: dict[dendropy.Node, list[int]] = {}
    attempt_of: dict[dendropy.Node, str] = {}  # which picks gave the age: 0, 1, ... or "higher_taxa"
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
                attempt_of[node] = str(attempt)
                break
        if node not in age and higher:
            sides = _higher_sides(node, higher, per_side * per_genus)
            value = timetree.split_age(sides)
            if value is not None:
                age[node], source[node] = value, "timetree_higher_taxa"
                raw[node], used[node] = value, sorted({i for side in sides for i in side})
                attempt_of[node] = "higher_taxa"

    root = tree.seed_node
    if root not in age:
        if root_age is None:
            raise ValueError("TimeTree has no age for the root; set backbone_build.root_age")
        age[root], source[root] = root_age, "configured"
    if interpolation == "pure_birth" and not richness:
        logging.getLogger(__name__).warning(
            "no clade richness for a pure-birth interpolation (it needs the static Open Tree files); "
            "spacing undated nodes evenly instead"
        )
        interpolation = "even"

    # the precedence rule (before #52), kept for comparison and as a fallback
    p_age, p_source, conflict, conflict_with = _precedence(tree, dict(age), dict(source))
    for node, value in interpolate(tree, p_age, method=interpolation, richness=richness).items():
        p_age[node] = value
        p_source[node] = {
            "dropped_inconsistent": "interpolated_after_drop",
            "timetree_same_as_ancestor": "interpolated_unresolved_in_timetree",
        }.get(p_source.get(node), "interpolated")

    # the calibrations and their uncertainty, whichever method, so that the conflict metrics of both
    # methods are computed with the same weights
    calibrations, cal_parent, unresolved = _calibrations(tree, age, source, used, timetree, synth_nodes, support,
                                                         fit)
    calibration_rows = _calibration_rows(calibrations)
    requery = dict(requery or {})
    rq_mode = {False: "off", True: "classify"}.get(requery.get("mode"), requery.get("mode", "classify"))
    if rq_mode not in ("classify", "apply", "off"):
        raise ValueError(f"requery mode must be classify, apply or off, not {rq_mode!r}")
    requery_rows: dict[dendropy.Node, dict[str, str]] = {}
    if method == "precedence":
        if rq_mode == "apply":
            raise ValueError("requery mode 'apply' requires age method 'fit'")
        age, source = p_age, p_source
        if rq_mode != "off":
            requery_rows, _ = _requery(tree, calibrations, {}, set(), conflict, conflict_with, p_age, used, ncbi,
                                       timetree, per_side, per_genus, requery, fit, synth_nodes, support)
    else:
        age, source, fitted_rows, result = _fit(tree, calibrations, cal_parent, unresolved, source, fit)
        if rq_mode != "off":
            requery_rows, replace = _requery(tree, calibrations, result.z, result.outlier, conflict, conflict_with,
                                             result.age, used, ncbi, timetree, per_side, per_genus, requery,
                                             fit, synth_nodes, support)
            if rq_mode == "apply" and replace:
                for node, (value, sigma) in replace.items():
                    calibrations[node] = agefit.Calibration(value, sigma, calibrations[node].source)
                    requery_rows[node]["requery_applied"] = "yes"
                # calibration_sigma keeps describing the first age (timetree_age); the replacement's is
                # requery_sigma, so that every reported metric pairs an age with its own uncertainty
                age, source, fitted_rows, result = _fit(tree, calibrations, cal_parent, unresolved, source, fit)
        for node, row in fitted_rows.items():
            calibration_rows[node].update(row)
        for node, value in interpolate(tree, age, method=interpolation, richness=richness).items():
            age[node] = value
            source[node] = ("interpolated_unresolved_in_timetree" if source.get(node) == "timetree_same_as_ancestor"
                            else "interpolated")

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
                "clade_genera": richness.get(node, "") if richness else "",
                "age": f"{age[node]:.4f}",
                "source": source[node],
                "timetree_age": f"{raw[node]:.4f}" if node in raw else "",
                # TimeTree age over that of the dated ancestor it conflicted with; large values
                # often point to a genus misplaced in the synthetic tree
                "conflict_ratio": f"{conflict[node]:.3f}" if node in conflict else "",
                "query_ids": "+".join(map(str, used.get(node, []))),
                "example_genera": ";".join(genera[:3]),
                "timetree_attempt": attempt_of.get(node, ""),
                "conflict_with": names[conflict_with[node]] if node in conflict_with else "",
            }
        )
        row = report[-1]
        evidence = timetree.evidence(used[node]) if node in used else {c: "" for c in EVIDENCE_COLUMNS}
        row.update({f"timetree_{c}": v for c, v in evidence.items()})
        label = (synth_nodes or {}).get(node, "")
        row["synth_node"] = label
        # a node without annotations rests on the taxonomy alone: no source tree either way
        counts = support.get(label, {}) if label and support is not None else None
        for kind in ("supported_by", "conflicts_with", "resolves", "partial_path_of", "terminal"):
            row[f"synth_{kind}"] = "" if counts is None else str(counts.get(kind, 0))
        row.update({c: "" for c in ("calibration_sigma", "calibration_weight", "fitted_age", "calibration_residual",
                                    "calibration_z", "calibration_outlier")})
        row.update(calibration_rows.get(node, {}))
        row.update(precedence_age=f"{p_age[node]:.4f}", precedence_source=p_source[node])
        row.update({c: "" for c in REQUERY_COLUMNS})
        row.update(requery_rows.get(node, {}))
        row.update({c: "" for c in RESERVED_COLUMNS})
    return tree, report


def _precedence(tree, age, source):
    """The conflict rule before #52, on copies of ``age`` and ``source``: drops calibrations in place.

    A calibration that is not younger than its nearest dated ancestor is
    dropped; equal ages mean TimeTree does not resolve the node. A direct
    TimeTree age outranks one from enclosing taxa; between equals, the
    ancestor is kept.
    """
    root = tree.seed_node

    def dated_ancestor(node):
        up = node.parent_node
        while up is not None and up not in age:
            up = up.parent_node
        return up

    conflict, conflict_with = {}, {}
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
                conflict_with[anc] = node
                del age[anc]
                source[anc] = "dropped_inconsistent"
            else:
                conflict[node] = age[node] / age[anc]
                conflict_with[node] = anc
                del age[node]
                source[node] = "dropped_inconsistent"
            changed = True
            break
    return age, source, conflict, conflict_with


def _cache_entry(timetree, ids) -> dict:
    entry = timetree.cache.get(timetree._key(ids)) if hasattr(timetree, "cache") else None
    return entry if isinstance(entry, dict) else {}


def _valid_estimates(entry: dict) -> list[float]:
    """The study estimates of a cache entry that are usable ages (TimeTree can return sentinels such as -1)."""
    return [e for e in entry.get("estimates") or () if e and e > 0]


def _sigma(entry: dict, node, higher_taxa: bool, synth_nodes, support, settings) -> float:
    """Uncertainty of a TimeTree age (#52), with the configured settings and the node's synthesis support."""
    label = (synth_nodes or {}).get(node, "")
    taxonomy_only = bool(label) and support is not None and not (support.get(label) or {}).get("supported_by")
    return agefit.calibration_sigma(
        _valid_estimates(entry), entry.get("ci_low"), entry.get("ci_high"),
        higher_taxa=higher_taxa, taxonomy_only=taxonomy_only,
        floor=float(settings.get("sigma_floor", agefit.SIGMA_FLOOR)),
        single=float(settings.get("sigma_single", agefit.SIGMA_SINGLE)),
        higher_factor=float(settings.get("higher_taxa_factor", agefit.HIGHER_TAXA_FACTOR)),
        taxonomy_factor=float(settings.get("taxonomy_only_factor", agefit.TAXONOMY_ONLY_FACTOR)))


def _calibrations(tree, age, source, used, timetree, synth_nodes, support, settings):
    """Every usable calibration with its uncertainty, whichever age method is used (#52).

    Returns the calibrations, each one's nearest calibrated ancestor (None
    for the root), and the nodes left out because their TimeTree age equals
    that of their nearest calibrated ancestor (TimeTree does not resolve
    them).
    """
    root = tree.seed_node
    calibrations: dict[dendropy.Node, agefit.Calibration] = {}
    parent: dict[dendropy.Node, dendropy.Node | None] = {}
    unresolved: set[dendropy.Node] = set()
    nearest: dict[dendropy.Node, dendropy.Node | None] = {}
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            continue
        up = node.parent_node
        anc = None if up is None else (up if up in calibrations else nearest.get(up))
        nearest[node] = anc
        if node not in age:
            continue
        if node is not root and anc is not None and abs(age[node] - calibrations[anc].age) <= 1e-9 * max(
                1.0, calibrations[anc].age):
            unresolved.add(node)
            continue
        if source[node] == "configured":
            calibrations[node] = agefit.Calibration(age[node], 1.0, "configured", fixed=True)
        else:
            sigma = _sigma(_cache_entry(timetree, used[node]), node, source[node] == "timetree_higher_taxa",
                           synth_nodes, support, settings)
            calibrations[node] = agefit.Calibration(age[node], sigma, source[node])
        parent[node] = anc
    return calibrations, parent, unresolved


def _calibration_rows(calibrations) -> dict:
    """The uncertainty columns of node_ages.tsv, which do not depend on the age method."""
    return {node: {"calibration_sigma": f"{c.sigma:.4f}", "calibration_weight": f"{c.weight:.4f}"}
            for node, c in calibrations.items() if not c.fixed}


def _fit(tree, calibrations, parent, unresolved, source, settings):
    """Ages of the calibrated nodes from a weighted fit of all calibrations (#52); the rest is left undated."""
    depth: dict[dendropy.Node, int] = {}
    for node in tree.preorder_node_iter():
        depth[node] = 0 if node.parent_node is None else depth[node.parent_node] + 1
    new_source = {n: s for n, s in source.items() if n.is_leaf()}
    new_source.update({n: "timetree_same_as_ancestor" for n in unresolved})
    result = agefit.fit_ages(parent, depth, calibrations, min_gap=float(settings.get("min_gap", agefit.MIN_GAP)),
                             huber_k=float(settings.get("huber_k", agefit.HUBER_K)),
                             outlier_z=float(settings.get("outlier_z", agefit.OUTLIER_Z)))
    fitted = {n: 0.0 for n in tree.leaf_node_iter()}
    rows = {}
    for node, c in calibrations.items():
        fitted[node] = result.age[node]
        if c.fixed:
            new_source[node] = "configured"
            continue
        suffix = "_outlier" if node in result.outlier else ""
        new_source[node] = f"{c.source}_fit{suffix}"
        rows[node] = {
            "calibration_weight": f"{result.weight[node]:.4f}",
            "fitted_age": f"{result.age[node]:.4f}",
            "calibration_residual": f"{result.residual[node]:.4f}",
            "calibration_z": f"{result.z[node]:.3f}",
            "calibration_outlier": "yes" if node in result.outlier else "",
        }
    return fitted, new_source, rows, result


def _bracket(node, reference) -> tuple[float | None, float | None]:
    """The ages the rest of the tree allows a node: the oldest of its nearest dated descendants (tips at 0)
    and its nearest dated ancestor, from ``reference`` (ages of the other nodes; the node's own is ignored)."""
    up = node.parent_node
    while up is not None and up not in reference:
        up = up.parent_node
    stack, low = list(node.child_node_iter()), 0.0
    while stack:
        n = stack.pop()
        if n.is_leaf():
            continue
        if n in reference:
            low = max(low, reference[n])
        else:
            stack.extend(n.child_node_iter())
    return low, (reference[up] if up is not None else None)


def _requery(tree, calibrations, z, outliers, conflict, conflict_with, reference, used, ncbi, timetree, per_side,
             per_genus, settings, fit_settings=None, synth_nodes=None, support=None):
    """Re-query triggered nodes with disjoint picks and classify them (#53).

    Returns the report columns per triggered node, and for the nodes whose
    class says the first age was the picks' doing (``under_spanning``,
    ``pick_specific``) the re-query age with its uncertainty, for ``apply``.
    """
    todo = rq.triggers(z, outliers, conflict, conflict_with, float(settings.get("z_trigger", rq.Z_TRIGGER)),
                       float(settings.get("ratio_trigger", rq.RATIO_TRIGGER)))
    rows, replace = {}, {}
    for node in tree.preorder_node_iter():
        if node not in todo or node not in calibrations or calibrations[node].fixed or node not in used:
            continue
        first = calibrations[node].age
        exclude = set(used[node])
        sides = [rq.disjoint_side([l.taxon.label for l in child.leaf_iter()], ncbi, exclude, per_side, per_genus)
                 for child in node.child_node_iter()]
        row = {c: "" for c in rq.COLUMNS}
        row["requery_trigger"] = todo[node]
        if sum(1 for s in sides if s) < 2:
            row.update(conflict_class="unresolved", conflict_detail="no_disjoint_picks")
            rows[node] = row
            continue
        # up to `attempts` disjoint pick sets, until TimeTree finds taxa on two sides
        again, ids = None, sorted({i for side in sides for i in side})
        for attempt in range(int(settings.get("attempts", 3))):
            if attempt:
                sides = [rq.disjoint_side([l.taxon.label for l in child.leaf_iter()], ncbi, exclude, per_side,
                                          per_genus, attempt) for child in node.child_node_iter()]
                if sum(1 for s in sides if s) < 2:
                    break
            ids = sorted({i for side in sides for i in side})
            again = timetree.split_age([s for s in sides if s])
            if again is not None:
                break
        evidence = timetree.evidence(ids) if hasattr(timetree, "evidence") else {}
        # only usable estimates count towards confidence, as towards the uncertainty
        n_first = len(_valid_estimates(_cache_entry(timetree, used[node])))
        n_again = len(_valid_estimates(_cache_entry(timetree, ids))) if again is not None else 0
        cls, detail, confidence = rq.classify(first, again, reference.get(node, first), n_first, n_again,
                                              float(settings.get("move", rq.MOVE)),
                                              bracket=_bracket(node, reference))
        row.update(requery_ids="+".join(map(str, ids)), requery_age="" if again is None else f"{again:.4f}",
                   requery_mrca_ncbi=evidence.get("mrca_ncbi", "") if again is not None else "",
                   requery_n_estimates=str(n_again) if again is not None else "",
                   conflict_class=cls, conflict_detail=detail, conflict_confidence=confidence)
        rows[node] = row
        if cls in ("under_spanning", "pick_specific") and again:
            # the same rules as for the first age: configured floor and single-estimate default, and the
            # node's synthesis support; disjoint picks are genera and species, never enclosing taxa
            sigma = _sigma(_cache_entry(timetree, ids), node, False, synth_nodes, support, fit_settings or {})
            row["requery_sigma"] = f"{sigma:.4f}"
            replace[node] = (again, sigma)
    return rows, replace


def age_fit_summary(report: list[dict], severe: float = 1.5) -> dict:
    """Conflict metrics over all calibrations, for the final ages and for the precedence rule (#52)."""
    calibrations, final, precedence = {}, {}, {}
    for row in report:
        # the calibrations of the fit, which carry their uncertainty in either method
        if not row.get("timetree_age") or not row.get("calibration_sigma"):
            continue
        calibrations[row["node"]] = agefit.Calibration(float(row["timetree_age"]), float(row["calibration_sigma"]))
        final[row["node"]] = float(row["age"])
        precedence[row["node"]] = float(row["precedence_age"])
    # with requery.mode apply, some calibrations were replaced by their re-query age (#53)
    effective = {n: agefit.Calibration(float(r["requery_age"]), float(r["requery_sigma"]))
                 for r in report for n in [r["node"]] if n in calibrations and r.get("requery_applied") == "yes"}
    return {
        "final": agefit.misfit(final, calibrations, severe),
        "final_against_requery": agefit.misfit(final, {**calibrations, **effective}, severe) if effective else None,
        "precedence": agefit.misfit(precedence, calibrations, severe),
        "outliers": sum(1 for r in report if r.get("calibration_outlier") == "yes"),
    }


# genera.tsv (#47). The first four columns are as before; the rest is evidence for the
# consistency layer (#46) and does not change how a genus is matched or placed.
GENERA_COLUMNS = (
    "genus", "status", "ott_id", "ncbi_id",
    # how the name was matched in OTT (static source only); ott_resolution is "col" where COL chose
    # between OTT homonyms (#48)
    "ott_resolution", "ott_candidates", "ott_candidates_in_context", "ott_candidate_ids", "ott_flags", "ott_sourceinfo",
    "ott_ncbi_ids", "ott_worms_ids", "ott_gbif_ids",
    # OTT lineage of the matched genus
    "ott_family", "ott_order", "ott_class", "ott_phylum",
    # where it sits in the synthetic tree
    "anchor_node", "anchor_clade_genera", "n_species_ott", "n_species_in_synthesis",
    # how a genus that is not a node of the synthesis was anchored on its species (#51)
    "anchor_rule", "anchor_fallback", "anchor_basis", "anchor_species", "anchor_species_inside", "anchor_outside_share",
    "anchor_outside_species", "anchor_mrca_node", "anchor_mrca_clade_genera",
    # lineage from the reference records: the most common name per rank, and its share
    "records_n", "records_family", "records_family_share", "records_order", "records_order_share",
    "records_class", "records_class_share", "records_phylum", "records_phylum_share",
    # COL XR as identifier hub (#48): the COL usage, its lineage, and the NCBI and OTT ids derived from it
    *COL_GENUS_COLUMNS,
    # the NCBI id the genus is sent to TimeTree with, and where it came from (#48)
    "join_ncbi_id", "join_ncbi_source",
    # placement concordance (#49; details in concordance.tsv): rank level at which the OTT placement and the
    # reference family part, and the merged NCBI ids part; share of barcode neighbours in the phyla of the
    # backbone sisters; flag, review, empty or not_evaluable
    "lineage_distance", "merged_id_distance", "neighbour_concordance", "concordance_flag",
    # re-anchoring of flagged genera (#50): what was done, from where to where, why, on whose authority
    "placement_action", "original_status", "original_anchor", "reanchor_target", "reanchor_how",
    "reanchor_trigger", "reanchor_source", "reanchor_note",
)
GENERA_RESERVED = ()
LINEAGE_RANKS = ("family", "order", "class", "phylum")


def records_lineage(records) -> dict[str, dict[str, str]]:
    """Per genus, the most common family, order, class and phylum among its records, with its share.

    The share is over the records of the genus that name that rank. This is
    the reference lineage that placement concordance compares OTT with (#49).
    """
    from collections import Counter, defaultdict

    counts: dict[str, dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
    n: Counter = Counter()
    for record in records:
        genus = (record.get("genus") or "").strip()
        if not genus:
            continue
        n[genus] += 1
        for rank in LINEAGE_RANKS:
            name = (record.get(rank) or "").strip()
            if name:
                counts[genus][rank][name] += 1
    out = {}
    for genus in n:
        row = {"records_n": str(n[genus])}
        for rank in LINEAGE_RANKS:
            c = counts[genus][rank]
            if c:
                name, votes = c.most_common(1)[0]
                row[f"records_{rank}"], row[f"records_{rank}_share"] = name, f"{votes / sum(c.values()):.3f}"
            else:
                row[f"records_{rank}"], row[f"records_{rank}_share"] = "", ""
        out[genus] = row
    return out


def genera_rows(matches: dict[str, dict], *extra: dict[str, dict]) -> list[dict]:
    """Rows of genera.tsv: the match of each genus, plus any per-genus columns in ``extra``."""
    rows = []
    for genus, match in sorted(matches.items()):
        row = {c: "" for c in GENERA_COLUMNS}
        row.update({k: ("" if v is None else v) for k, v in match.items() if k in GENERA_COLUMNS})
        for table in extra:
            row.update({k: v for k, v in (table.get(genus) or {}).items() if k in GENERA_COLUMNS})
        row["genus"] = genus
        rows.append(row)
    return rows


def timetree_summary(report: list[dict], crowns: list[dict] | None = None) -> dict:
    """How much TimeTree evidence the outputs carry, for metadata.json."""
    from collections import Counter

    evidence, versions = Counter(), Counter()
    for row in list(report) + list(crowns or []):
        if row.get("query_ids"):
            evidence[row.get("timetree_evidence") or "none"] += 1
            if row.get("timetree_api_version"):
                versions[row["timetree_api_version"]] += 1
    return {"queries_by_evidence": dict(evidence), "api_versions": dict(versions)}


JOIN_MODES = ("ott", "col")
JOIN_CONFIDENCE = ("high", "medium")


def join_ncbi(ncbi: dict[str, object], matches: dict[str, dict], col_rows: dict[str, dict],
              info: dict[str, dict] | None = None, mode: str = "col") -> dict[str, dict[str, str]]:
    """Choose the NCBI id each genus is sent to TimeTree with; changes ``ncbi`` in place (#48).

    With ``mode`` ``ott`` the id is OTT's, as before. With ``col`` a COL-derived
    id of ``high`` or ``medium`` confidence replaces it where OTT has none, or
    where it is another of the ids OTT merged into the genus (OTT lists them in
    no particular order, and the first is not always right: for Onchidella it
    is the sand dollar Peronella). A COL id that OTT does not hold at all is
    not used, only reported (``col_ncbi_vs_ott`` = ``different``), since that
    is a disagreement for placement concordance (#49) to weigh, not a join.
    Species ids are left alone. Returns the ``join_ncbi_*`` and
    ``col_ncbi_vs_ott`` columns per genus.
    """
    if mode not in JOIN_MODES:
        raise ValueError(f"unknown col.join {mode!r}; use one of {', '.join(JOIN_MODES)}")
    out = {}
    for genus, match in matches.items():
        ott_first = "" if match.get("ncbi_id") in (None, "") else str(match["ncbi_id"])
        ott_ids = (info or {}).get(genus, {}).get("ott_ncbi_ids", "")
        col = col_rows.get(genus) or {}
        relation = compare_ids(col.get("col_ncbi_id", ""), ott_ids, ott_first)
        row = {"join_ncbi_id": ott_first, "join_ncbi_source": "ott" if ott_first else "", "col_ncbi_vs_ott": relation}
        usable = (mode == "col" and col.get("col_confidence") in JOIN_CONFIDENCE and col.get("col_ncbi_id")
                  and match.get("ott_id") not in (None, ""))
        if usable and (relation == "merged" or not ott_first):
            new = int(col["col_ncbi_id"])
            value = ncbi.get(genus)
            if isinstance(value, (list, tuple)):
                ncbi[genus] = [new] + list(value[1:])
            else:
                ncbi[genus] = new
            row.update(join_ncbi_id=str(new), join_ncbi_source="col")
        out[genus] = row
    return out
