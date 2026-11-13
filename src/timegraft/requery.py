"""Re-querying conflicting backbone nodes with disjoint picks, to tell causes apart (#53).

A TimeTree age that conflicts with the rest of the tree has two common
causes, and they leave different traces. A query whose picks do not
straddle the intended node in TimeTree's own tree dates a shallower node:
the age is too young, never too old (*under-spanning*). A taxon placed
differently by OToL and TimeTree (a misplaced genus, a genus that is not a
natural group, a wrong id) makes TimeTree date a deeper split than the one
in the backbone, or a shallower one, however it is asked (*topology*).

Asking again with picks that share no taxon with the first query separates
them: an under-spanned node comes back older, closer to what the rest of
the tree says; a topological conflict comes back with the same age; and
an old age that does not come back points to one of the first picks
(*pick_specific*). Everything is cached, so the classification is
reproducible on a fixed cache.
"""

from __future__ import annotations

import math
from typing import Mapping, Sequence

CLASSES = ("under_spanning", "pick_specific", "topology", "consistent", "unresolved")
COLUMNS = ("requery_trigger", "requery_ids", "requery_age", "requery_mrca_ncbi", "requery_n_estimates",
           "conflict_class", "conflict_detail", "conflict_confidence", "requery_applied")
Z_TRIGGER = 2.0
RATIO_TRIGGER = 1.5
MOVE = 1.25  # a re-query age within this factor of the first counts as the same
STRONG = 1.5  # a move by this factor or more, or a match within 1.1, counts towards high confidence


def disjoint_side(leaf_genera: Sequence[str], ncbi: Mapping[str, object], exclude: set[int], per_side: int = 2,
                  per_genus: int = 3, attempt: int = 0) -> list[int]:
    """NCBI ids for one side of a node that share none with ``exclude`` (the first query).

    Genera not used before come first, spread over the clade; where every
    genus was used, other species of the same genera are taken. Each
    ``attempt`` shifts to other genera, and to the next species of each.
    """
    genera = sorted(g for g in leaf_genera if ncbi.get(g))
    if not genera:
        return []

    def ids(genus):
        value = ncbi[genus]
        if isinstance(value, (list, tuple)):
            return [int(v) for v in value if v]
        return [int(value)] if value else []

    fresh = [g for g in genera if not set(ids(g)) & exclude]
    pool = fresh or genera
    if len(pool) > per_side:
        step = len(pool) / per_side
        pool = sorted({pool[(int(i * step) + attempt) % len(pool)] for i in range(per_side)})
        shift = 0
    else:
        shift = attempt  # no other genera to turn to: other species instead
    out = []
    for genus in pool:
        unused = [i for i in ids(genus) if i not in exclude]
        out.extend(unused[shift * (1 + per_genus): (shift + 1) * (1 + per_genus)])
    return out


def triggers(z: Mapping[object, float], outliers: set, conflict: Mapping[object, float],
             conflict_with: Mapping[object, object], z_trigger: float = Z_TRIGGER,
             ratio_trigger: float = RATIO_TRIGGER) -> dict[object, str]:
    """Nodes to re-query and why: a fit residual of ``z_trigger`` sigma or more (``z``), an outlier of the
    fit (``outlier``), or a conflict under the precedence rule of ``ratio_trigger`` or more, which takes in
    both nodes of the conflict (``conflict``)."""
    out: dict[object, list[str]] = {}
    for node, value in z.items():
        if abs(value) >= z_trigger:
            out.setdefault(node, []).append("z")
    for node in outliers:
        out.setdefault(node, []).append("outlier")
    for node, ratio in conflict.items():
        if ratio >= ratio_trigger:
            for n in (node, conflict_with.get(node)):
                if n is not None:
                    out.setdefault(n, []).append("conflict")
    return {n: ";".join(dict.fromkeys(reasons)) for n, reasons in out.items()}


def classify(first: float, again: float | None, reference: float, n_first: int = 0, n_again: int = 0,
             move: float = MOVE, strong: float = STRONG) -> tuple[str, str, str]:
    """Class, detail and confidence for one re-queried node.

    ``first`` is TimeTree's first age, ``again`` the re-query's (None if it
    gave none), ``reference`` what the rest of the tree says (the fitted age,
    or the precedence age); ``n_first`` and ``n_again`` the estimates behind
    each. ``detail`` is ``moved_older``, ``moved_younger``, ``stable_young``,
    ``stable_old``, ``close_to_reference`` (class ``consistent``: the first
    age is within ``move`` of the reference and the re-query agrees, which
    happens to the partner of a conflict whose other side was at fault),
    ``contrary`` (moved away from the reference) or ``no_answer``;
    confidence is ``high``, ``medium`` or ``low``.
    """
    if again is None or again <= 0:
        return "unresolved", "no_answer", ""
    shift = math.log(again / first)
    young = reference > first  # the first age is younger than the rest of the tree says
    if abs(math.log(reference / first)) < math.log(move) and abs(shift) < math.log(move):
        # within the factor of what the rest of the tree says, and asked again the same: no conflict to explain
        return "consistent", "close_to_reference", "high" if n_first >= 2 and n_again >= 2 else "medium"
    if abs(shift) < math.log(move):
        cls, detail = "topology", "stable_young" if young else "stable_old"
        decisive = abs(shift) < math.log(1.1)
    elif shift > 0 and young and abs(math.log(again / reference)) < abs(math.log(first / reference)):
        cls, detail = "under_spanning", "moved_older"
        decisive = shift >= math.log(strong)
    elif shift < 0 and not young and abs(math.log(again / reference)) < abs(math.log(first / reference)):
        cls, detail = "pick_specific", "moved_younger"
        decisive = -shift >= math.log(strong)
    else:
        return "unresolved", "contrary", ""
    supported = n_first >= 2 and n_again >= 2
    confidence = "high" if decisive and supported else "medium" if decisive or supported else "low"
    return cls, detail, confidence
