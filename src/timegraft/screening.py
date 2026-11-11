"""Screening reference units for misidentifications and outliers (issue #38).

Within-genus barcode distances set the lower bound of a genus's crown age
(issue #36). A unit that does not belong to its genus, through
misidentification, contamination or a genus that is not a natural group,
inflates those distances. Before genus subtrees are built, every unit
representative is compared with all others, and a unit is set aside when:

- ``other_genus``: it is at least ``margin`` more similar (in identity) to a
  unit of another genus than to any unit of its own genus; or
- ``isolated``: no unit of its own genus reaches ``min_identity``, while the
  genus has units that do (a core to be isolated from).

A genus never loses all its units: when every unit of a genus would be set
aside, all are kept and the decision says so. Genera with one unit have
nothing to compare with and are not screened.

Set-aside units are not discarded. They are placed on the tree like units
without a usable genus, so they remain in diversity, but they no longer
shape the genus subtree or its crown age.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Mapping

from .placement import vsearch_hits

COLUMNS = (
    "unit_id", "genus", "genus_units", "n_records", "decision", "best_congener", "best_congener_identity",
    "best_other", "best_other_genus", "best_other_identity",
)


def screen(
    genus_of: Mapping[str, str],
    hits: Mapping[str, Iterable[tuple[str, float]]],
    margin: float = 0.02,
    min_identity: float = 0.80,
) -> list[dict]:
    """Decide per unit: ``keep``, ``other_genus``, ``isolated`` or ``kept_whole_genus``.

    ``genus_of`` maps every unit to its genus; ``hits`` maps a unit to
    (unit, identity) pairs from a search of all units against each other
    (self hits are ignored), with identity in 0..1.
    """
    members = defaultdict(list)
    for unit, genus in genus_of.items():
        members[genus].append(unit)
    best_same: dict[str, tuple[str, float]] = {}
    best_other: dict[str, tuple[str, float]] = {}
    for unit in genus_of:
        for target, ident in hits.get(unit, ()):
            if target == unit or target not in genus_of:
                continue
            table = best_same if genus_of[target] == genus_of[unit] else best_other
            # a hit is evidence for both units: heuristic searches do not always find it both ways
            for a, b in ((unit, target), (target, unit)):
                if ident > table.get(a, ("", 0.0))[1]:
                    table[a] = (b, ident)

    decision = {}
    for genus, units in members.items():
        if len(units) < 2:
            for unit in units:
                decision[unit] = "keep"
            continue
        has_core = any(best_same.get(u, ("", 0.0))[1] >= min_identity for u in units)
        for unit in units:
            same = best_same.get(unit, ("", 0.0))[1]
            other = best_other.get(unit, ("", 0.0))[1]
            if other and other >= same + margin:
                decision[unit] = "other_genus"
            elif same < min_identity and has_core:
                decision[unit] = "isolated"
            else:
                decision[unit] = "keep"
        if all(decision[u] != "keep" for u in units):
            for unit in units:
                decision[unit] = "kept_whole_genus"

    rows = []
    for unit, genus in genus_of.items():
        same_unit, same = best_same.get(unit, ("", 0.0))
        other_unit, other = best_other.get(unit, ("", 0.0))
        rows.append({
            "unit_id": unit, "genus": genus, "genus_units": len(members[genus]), "decision": decision[unit],
            "best_congener": same_unit, "best_congener_identity": f"{same:.4f}" if same_unit else "",
            "best_other": other_unit, "best_other_genus": genus_of.get(other_unit, "") if other_unit else "",
            "best_other_identity": f"{other:.4f}" if other_unit else "",
        })
    return rows


def search_units(
    reps: list[tuple[str, str]],
    min_identity: float,
    query_cov: float = 0.7,
    maxaccepts: int = 10,
    threads: int = 1,
) -> dict[str, list[tuple[str, float]]]:
    """All representatives against each other with VSEARCH; self hits removed."""
    hits = vsearch_hits(reps, reps, min_identity=min_identity, query_cov=query_cov,
                        maxaccepts=maxaccepts + 1, threads=threads)
    return {q: [(t, i) for t, i in h if t != q] for q, h in hits.items()}


def set_aside(rows: Iterable[Mapping[str, str]]) -> set[str]:
    """Units the screening takes out of their genus subtree."""
    return {r["unit_id"] for r in rows if r["decision"] in ("other_genus", "isolated")}
