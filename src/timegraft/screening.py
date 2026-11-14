"""Screening reference units for misidentifications and outliers (issues #38, #55).

Within-genus barcode distances set the lower bound of a genus's crown age
(issue #36). A unit that does not belong to its genus, through
misidentification, contamination or a genus that is not a natural group,
inflates those distances. Before genus subtrees are built, every unit
representative is compared with all others, and a unit is set aside when:

- ``other_genus``: it is at least ``margin`` more similar (in identity) to a
  unit of another genus than to any unit of its own genus; or
- ``isolated``: no unit of its own genus reaches ``min_identity``, while the
  genus has units that do (a core to be isolated from); or
- ``numt_like``: only when ``translation.action`` is ``set_aside``, its codons
  differ from those of its nearest relative in a way that real mitochondrial
  sequences rarely do (see ``compare_codons``).

A genus never loses all its units: when every unit of a genus would be set
aside, all are kept and the decision says so. Genera with one unit have
nothing to compare with and are not screened.

Issue #55 widened this in four ways, all recorded per unit in the output:

- **Reference scope.** Besides the units of the run, sequences from an
  external reference (``reference.path``) can be searched, restricted to the
  families (or another rank) present in the run. They add evidence (a closer
  congener, or a closer unit of another genus) but are never decided on. Each
  best hit says whether it came from the ``run`` or the ``external`` reference.
- **Thresholds per class.** ``margin`` and ``min_identity`` can differ per
  taxon at a chosen rank (``classes``, by default the class); the values used
  for each unit are written next to its decision. ``calibrate`` summarises the
  observed identities per class, to set those thresholds from the data.
- **Codon diagnostics.** For coding markers, the unit is compared codon by
  codon with its nearest relative, in its own genetic code: how many codons
  differ, how many of those change the amino acid, and whether the alignment
  is out of frame (an indel whose length is not a multiple of three). Many
  amino-acid changes point to a nuclear copy (NUMT); a frameshift to a poor
  sequence in either of the two.
- **Reasons.** Every decision comes with the comparison that made it, in words.

Set-aside units are not discarded. They are placed on the tree like units
without a usable genus, so they remain in diversity, but they no longer
shape the genus subtree or its crown age.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Iterable, Mapping

from . import geneticcode
from .placement import vsearch_search

COLUMNS = (
    "unit_id", "genus", "genus_units", "n_records", "decision",
    "best_congener", "best_congener_identity", "best_congener_scope",
    "best_other", "best_other_genus", "best_other_identity", "best_other_scope",
    "reference_scope", "threshold_class", "margin", "min_identity", "reason",
    "translation", "translation_ref", "translation_reason", "codons_compared", "codons_changed", "aa_changes",
    "nonsynonymous_share", "frameshift_stretches",
)

SET_ASIDE = ("other_genus", "isolated", "numt_like")

DEFAULT_TRANSLATION = {"min_changed_codons": 8, "max_nonsynonymous_share": 0.5, "action": "flag"}


# Thresholds ------------------------------------------------------------------


def resolve_thresholds(
    class_of: Mapping[str, str] | None,
    units: Iterable[str],
    margin: float,
    min_identity: float,
    thresholds: Mapping[str, Mapping[str, float]] | None = None,
) -> dict[str, tuple[str, float, float]]:
    """Per unit: (class whose thresholds apply, or ``default``; margin; min_identity)."""
    thresholds = thresholds or {}
    out = {}
    for unit in units:
        name = (class_of or {}).get(unit, "")
        if name and name in thresholds:
            entry = thresholds[name]
            out[unit] = (name, float(entry.get("margin", margin)), float(entry.get("min_identity", min_identity)))
        else:
            out[unit] = ("default", float(margin), float(min_identity))
    return out


def read_threshold_table(path: str) -> dict[str, dict[str, float]]:
    """A table with columns ``class``, ``margin`` and ``min_identity``; empty cells use the default."""
    from .seqio import read_tsv

    out = {}
    for row in read_tsv(path):
        entry = {k: float(row[k]) for k in ("margin", "min_identity") if (row.get(k) or "").strip()}
        if entry:
            out[row["class"]] = entry
    return out


# Decisions -------------------------------------------------------------------


def _fmt(value: float) -> str:
    return f"{value:.4f}"


def screen(
    genus_of: Mapping[str, str],
    hits: Mapping[str, Iterable[tuple[str, float]]],
    margin: float = 0.02,
    min_identity: float = 0.80,
    external: Mapping[str, str] | None = None,
    class_of: Mapping[str, str] | None = None,
    thresholds: Mapping[str, Mapping[str, float]] | None = None,
    reference_scope: str = "run",
    mutual_pairs: str = "both",
) -> list[dict]:
    """Decide per unit: ``keep``, ``other_genus``, ``isolated`` or ``kept_whole_genus``.

    ``genus_of`` maps every unit of the run to its genus; ``hits`` maps a unit
    to (unit, identity) pairs from a search of all units against each other
    and against the external reference (self hits are ignored), with identity
    in 0..1. ``external`` maps the identifiers of external reference sequences
    to their genus: they count as evidence but are not decided on.
    ``class_of`` and ``thresholds`` give per-class ``margin`` and
    ``min_identity`` (see ``resolve_thresholds``); ``reference_scope`` is
    written to every row (e.g. ``run`` or ``run+external:family``).

    ``mutual_pairs`` decides about two units of different genera that are
    each other's nearest unit of another genus, and both ``other_genus``: a
    pair cannot tell which of the two carries the wrong name. ``both`` sets
    both aside; ``further`` keeps the one that is closer to its own genus
    (the smaller excess of the other-genus identity over its best congener),
    and sets aside the other. Equal excesses set aside both.
    """
    if mutual_pairs not in ("both", "further"):
        raise ValueError(f"mutual_pairs must be both or further, not {mutual_pairs!r}")
    external = external or {}
    genus_any = {**external, **genus_of}
    members = defaultdict(list)
    for unit, genus in genus_of.items():
        members[genus].append(unit)
    best_same: dict[str, tuple[str, float]] = {}
    best_other: dict[str, tuple[str, float]] = {}
    for unit in genus_of:
        for target, ident in hits.get(unit, ()):
            if target == unit or target not in genus_any:
                continue
            table = best_same if genus_any[target] == genus_of[unit] else best_other
            # a hit is evidence for both units: heuristic searches do not always find it both ways
            for a, b in ((unit, target), (target, unit)):
                if a in genus_of and ident > table.get(a, ("", 0.0))[1]:
                    table[a] = (b, ident)

    limits = resolve_thresholds(class_of, genus_of, margin, min_identity, thresholds)
    decision, reason = {}, {}
    for genus, units in members.items():
        if len(units) < 2:
            for unit in units:
                decision[unit] = "keep"
                reason[unit] = "only unit of its genus in the run: not screened"
            continue
        has_core = any(best_same.get(u, ("", 0.0))[1] >= limits[u][2] for u in units)
        for unit in units:
            _, unit_margin, unit_min = limits[unit]
            same_unit, same = best_same.get(unit, ("", 0.0))
            other_unit, other = best_other.get(unit, ("", 0.0))
            other_text = f"{genus_any[other_unit]} at {_fmt(other)}" if other_unit else "none"
            same_text = _fmt(same) if same_unit else "none"
            if other and other >= same + unit_margin:
                decision[unit] = "other_genus"
                reason[unit] = (f"closer to another genus ({other_text}) than to any congener ({same_text}) "
                                f"by {_fmt(other - same)} >= margin {unit_margin:g}")
            elif same < unit_min and has_core:
                decision[unit] = "isolated"
                found = f"best congener {same_text} below" if same_unit else "no congener at"
                reason[unit] = (f"{found} min_identity {unit_min:g}, while other units of the genus reach it")
            else:
                decision[unit] = "keep"
                if same < unit_min:
                    reason[unit] = f"no congener at {unit_min:g}, but neither has any unit of the genus"
                else:
                    reason[unit] = f"best congener {same_text}; nearest other genus {other_text}"
    if mutual_pairs == "further":
        excess = {u: best_other[u][1] - best_same.get(u, ("", 0.0))[1] for u in decision if u in best_other}
        for unit in list(decision):
            other_unit = best_other.get(unit, ("", 0.0))[0]
            if (decision[unit] == "other_genus" and other_unit in genus_of and decision.get(other_unit) == "other_genus"
                    and best_other.get(other_unit, ("", 0.0))[0] == unit and excess[unit] < excess[other_unit]):
                decision[unit] = "keep"
                reason[unit] = (f"nearest to {other_unit} of {genus_of[other_unit]}, and it to this unit; "
                                f"kept as the closer of the two to its own genus "
                                f"(excess {_fmt(excess[unit])} against {_fmt(excess[other_unit])})")
    _guard_genera(members, decision, reason)

    rows = []
    for unit, genus in genus_of.items():
        same_unit, same = best_same.get(unit, ("", 0.0))
        other_unit, other = best_other.get(unit, ("", 0.0))
        name, unit_margin, unit_min = limits[unit]
        rows.append({
            "unit_id": unit, "genus": genus, "genus_units": len(members[genus]), "decision": decision[unit],
            "best_congener": same_unit, "best_congener_identity": _fmt(same) if same_unit else "",
            "best_congener_scope": _scope(same_unit, external),
            "best_other": other_unit, "best_other_genus": genus_any.get(other_unit, "") if other_unit else "",
            "best_other_identity": _fmt(other) if other_unit else "",
            "best_other_scope": _scope(other_unit, external),
            "reference_scope": reference_scope, "threshold_class": name,
            "margin": f"{unit_margin:g}", "min_identity": f"{unit_min:g}", "reason": reason[unit],
        })
    return rows


def _scope(unit: str, external: Mapping[str, str]) -> str:
    if not unit:
        return ""
    return "external" if unit in external else "run"


def _guard_genera(members: Mapping[str, list[str]], decision: dict, reason: dict) -> None:
    """A genus never loses all its units: if all would be set aside, all are kept."""
    for genus, units in members.items():
        if len(units) > 1 and all(decision[u] in SET_ASIDE for u in units):
            for unit in units:
                reason[unit] = (f"would be {decision[unit]} ({reason[unit]}), but so would every unit "
                                f"of the genus: all kept")
                decision[unit] = "kept_whole_genus"


def apply_translation(rows: list[dict], action: str = "flag") -> None:
    """Set aside ``numt_like`` units when ``action`` is ``set_aside``; otherwise only flag them.

    Only units that would otherwise be kept change, and the whole-genus rule
    still holds. Changes ``rows`` in place.
    """
    if action not in ("flag", "set_aside"):
        raise ValueError(f"translation.action must be flag or set_aside, not {action!r}")
    if action != "set_aside":
        return
    by_unit = {r["unit_id"]: r for r in rows}
    members = defaultdict(list)
    for r in rows:
        members[r["genus"]].append(r["unit_id"])
    changed = set()
    for r in rows:
        if r.get("translation") == "numt_like" and r["decision"] == "keep" and len(members[r["genus"]]) > 1:
            r["decision"] = "numt_like"
            r["reason"] = f"codons unlike those of {r['translation_ref']}: {r['translation_reason']}"
            changed.add(r["genus"])
    for genus in changed:
        units = members[genus]
        if all(by_unit[u]["decision"] in SET_ASIDE for u in units):
            decision = {u: by_unit[u]["decision"] for u in units}
            reason = {u: by_unit[u]["reason"] for u in units}
            _guard_genera({genus: units}, decision, reason)
            for u in units:
                by_unit[u]["decision"], by_unit[u]["reason"] = decision[u], reason[u]


def set_aside(rows: Iterable[Mapping[str, str]]) -> set[str]:
    """Units the screening takes out of their genus subtree."""
    return {r["unit_id"] for r in rows if r["decision"] in SET_ASIDE}


# Search ----------------------------------------------------------------------


def search_units(
    reps: list[tuple[str, str]],
    min_identity: float,
    query_cov: float = 0.7,
    maxaccepts: int = 10,
    threads: int = 1,
    external: list[tuple[str, str]] | None = None,
    genus_of: Mapping[str, str] | None = None,
    alignments: bool | None = None,
    deep: int = 100,
) -> dict[str, list[tuple[str, float]]] | tuple[dict, dict]:
    """All representatives against each other (and the external reference) with VSEARCH.

    Self hits are removed. VSEARCH keeps the first ``maxaccepts`` hits of a
    query, so in a large genus its congeners can crowd out every unit of
    another genus, or the reverse, and the evidence from that side is lost.
    With ``genus_of`` (genus per unit and per external sequence), queries
    whose hits fill ``maxaccepts`` while one side (own genus or other genera)
    is empty are searched again with ``deep`` accepts; the extra hits are
    added. With ``alignments`` (by default, whenever ``genus_of`` is given),
    the alignment of each query with its best hit on each side is kept too,
    and ``(hits, alignments)`` is returned; ``alignments`` maps (query,
    target) to (query start, query row, target row).
    """
    if alignments is None:
        alignments = genus_of is not None
    if alignments and genus_of is None:
        raise ValueError("alignments need genus_of")
    refs = list(reps) + list(external or [])
    hits: dict[str, dict[str, float]] = defaultdict(dict)
    rows: dict[tuple[str, str], tuple[int, str, str]] = {}
    best: dict[tuple[str, bool], tuple[str, float]] = {}
    fields = ("qilo", "qrow", "trow") if alignments else ()

    def collect(queries, accepts):
        for query, target, ident, *extra in vsearch_search(queries, refs, min_identity, query_cov, accepts,
                                                           threads, fields):
            if target == query:
                continue
            hits[query][target] = max(ident, hits[query].get(target, 0.0))
            if alignments:
                key = (query, genus_of.get(target) == genus_of.get(query))
                if ident > best.get(key, ("", -1.0))[1]:
                    old = best.get(key)
                    if old:
                        rows.pop((query, old[0]), None)
                    best[key] = (target, ident)
                    rows[(query, target)] = (int(extra[0]), extra[1], extra[2])

    collect(reps, maxaccepts + 1)
    if genus_of is not None and deep > maxaccepts:
        def one_sided(query):
            sides = {genus_of.get(t) == genus_of.get(query) for t in hits.get(query, {})}
            return len(hits.get(query, {})) >= maxaccepts and len(sides) == 1
        again = [(q, seq) for q, seq in reps if one_sided(q)]
        if again:
            collect(again, deep + 1)
    out = {q: sorted(h.items(), key=lambda kv: -kv[1]) for q, h in hits.items()}
    if not alignments:
        return out
    return out, rows


def in_scope(
    reference: Iterable[Mapping[str, str]],
    run_taxa: Mapping[str, set[str]],
    scope: str = "family",
    run_units: Iterable[str] = (),
) -> list[Mapping[str, str]]:
    """External reference rows within the scope of the run.

    ``scope`` is a rank (``genus``, ``family``, ``order``, ``class``,
    ``phylum``) or ``all``: a row is kept when its name at that rank occurs in
    the run (``run_taxa[rank]``). Rows for units that are in the run are
    dropped, so that the reference adds evidence without repeating it.
    """
    run_units = set(run_units)
    out = []
    for row in reference:
        if row.get("unit_id") and row["unit_id"] in run_units:
            continue
        if scope != "all" and (row.get(scope) or "") not in run_taxa.get(scope, set()):
            continue
        out.append(row)
    return out


# Codon diagnostics -----------------------------------------------------------


def compare_codons(qilo: int, qrow: str, trow: str, query: str, code: int) -> dict:
    """Compare a query with an aligned target codon by codon, in the query's reading frame.

    ``qilo`` is the 1-based start of the alignment on the (ungapped) query,
    ``qrow`` and ``trow`` the aligned rows, ``query`` the whole query.
    Codons are compared where all three bases of the query align to bases of
    the target without a gap between them, and both are unambiguous.

    The frame is the forward frame of the query with the fewest stop codons.
    Where several frames have equally few (common in codes with only two
    stops, such as the echinoderm code), the one with the fewest amino-acid
    changes is taken: of the candidate frames, the reading frame is the one in
    which the differences are most often silent, so this never makes a pair
    look more like a nuclear copy than it is.

    Returns codons compared and changed, amino-acid changes, the frame, and
    the number of stretches where the alignment is out of frame (see
    ``_frameshifts``).
    """
    query = geneticcode._degap(query)
    counts = geneticcode.frame_stops(query, code)
    table = geneticcode._MAPS[geneticcode._check(code)]
    column = {}
    pos = qilo - 1
    for col, base in enumerate(qrow):
        if base != "-":
            column[pos] = col
            pos += 1
    best = None
    for frame in (f for f in range(3) if counts[f] == min(counts)):
        compared = changed = aa = 0
        for start in range(frame, len(query) - 2, 3):
            cols = [column.get(start + k) for k in range(3)]
            if None in cols or cols[2] - cols[0] != 2:
                continue
            qc = "".join(qrow[c] for c in cols).upper()
            tc = "".join(trow[c] for c in cols).upper()
            if qc not in table or tc not in table:
                continue
            compared += 1
            if qc != tc:
                changed += 1
                aa += table[qc] != table[tc]
        if best is None or (aa, -compared) < (best["aa_changes"], -best["codons_compared"]):
            best = {"codons_compared": compared, "codons_changed": changed, "aa_changes": aa, "frame": frame}
    best["frameshift_stretches"] = _frameshifts(qrow, trow)
    return best


def _frameshifts(qrow: str, trow: str, window: int = 15) -> int:
    """Stretches of the alignment that are out of frame between query and target.

    The offset is the running difference between gaps in the two rows. An
    indel shifts the frame when the offset is not a multiple of three for more
    than ``window`` columns: two nearby indels that cancel each other out are
    an artefact of the alignment, not a frameshift.
    """
    count, offset, since = 0, 0, None
    for col, (q, t) in enumerate(zip(qrow, trow)):
        offset += (t == "-") - (q == "-")
        if offset % 3:
            if since is None:
                since = col
            elif col - since == window:
                count += 1
        else:
            since = None
    return count


def judge_codons(stats: Mapping[str, int], settings: Mapping | None = None) -> tuple[str, str]:
    """``numt_like``, ``frameshift`` or ``ok``, and why.

    ``numt_like``: at least ``min_changed_codons`` codons differ and more than
    ``max_nonsynonymous_share`` of them change the amino acid. Diverging
    mitochondrial copies change mostly third positions, which rarely change
    the amino acid; a nuclear copy is free of that constraint.
    ``frameshift``: the alignment is out of frame over a stretch. A pair cannot
    tell which of the two sequences has the indel, so this is reported only.
    """
    cfg = {**DEFAULT_TRANSLATION, **(settings or {})}
    changed, aa = stats["codons_changed"], stats["aa_changes"]
    share = aa / changed if changed else 0.0
    text = f"{aa} of {changed} changed codons change the amino acid"
    if changed >= cfg["min_changed_codons"] and share > cfg["max_nonsynonymous_share"]:
        return "numt_like", f"{text} (share {share:.2f} > {cfg['max_nonsynonymous_share']:g})"
    if stats["frameshift_stretches"]:
        return "frameshift", f"{text}; out of frame over {stats['frameshift_stretches']} stretch(es)"
    return "ok", text


def assess_translation(
    rows: list[dict],
    alignments: Mapping[tuple[str, str], tuple[int, str, str]],
    sequences: Mapping[str, str],
    code_of: Mapping[str, int],
    settings: Mapping | None = None,
) -> None:
    """Fill the codon columns of ``rows`` in place, against the best congener or else the best other hit.

    The alignment of a pair may have been kept from either side; from the
    relative's side, the comparison runs in the relative's frame and genetic
    code (the unit's code if the relative has none, as for an external
    sequence). Units without an alignment
    are ``not_assessed``.
    """
    for row in rows:
        unit = row["unit_id"]
        ref = row.get("best_congener") or row.get("best_other") or ""
        row.update({"translation": "not_assessed", "translation_ref": ref, "translation_reason": ""})
        for key in ("codons_compared", "codons_changed", "aa_changes", "nonsynonymous_share", "frameshift_stretches"):
            row[key] = ""
        if not ref:
            continue
        if (unit, ref) in alignments:
            qilo, qrow, trow = alignments[(unit, ref)]
            stats = compare_codons(qilo, qrow, trow, sequences[unit], code_of.get(unit, 5))
        elif (ref, unit) in alignments and ref in sequences:
            # the alignment as seen from the relative: read in its frame and its genetic code
            qilo, qrow, trow = alignments[(ref, unit)]
            stats = compare_codons(qilo, qrow, trow, sequences[ref], code_of.get(ref, code_of.get(unit, 5)))
        else:
            continue
        verdict, why = judge_codons(stats, settings)
        changed = stats["codons_changed"]
        row.update({k: str(v) for k, v in stats.items()})
        row["nonsynonymous_share"] = f"{stats['aa_changes'] / changed:.3f}" if changed else ""
        row["translation"], row["translation_reason"] = verdict, why


def not_coding(rows: list[dict]) -> None:
    """Mark the codon columns for markers without a reading frame."""
    for row in rows:
        row["translation"] = "not_coding"


# Calibration -----------------------------------------------------------------

CALIBRATION_COLUMNS = (
    "class", "n_units", "n_genera", "screened_units", "congener_p01", "congener_p05", "congener_p50",
    "other_closer", "other_closer_share", "set_aside", "set_aside_share", "nonsynonymous_p95", "numt_like",
    "suggested_min_identity",
)


def _quantile(values: list[float], q: float) -> float:
    values = sorted(values)
    if not values:
        return math.nan
    pos = q * (len(values) - 1)
    lo, hi = math.floor(pos), math.ceil(pos)
    return values[lo] + (values[hi] - values[lo]) * (pos - lo)


def calibrate(
    rows: Iterable[Mapping[str, str]],
    class_of: Mapping[str, str],
    floor: float = 0.80,
    ceiling: float = 0.95,
    quantile: float = 0.01,
    min_units: int = 50,
) -> list[dict]:
    """Observed identities per class, and a ``min_identity`` that follows from them.

    For units in screened genera (two or more units) with a congener hit, the
    distribution of best-congener identity is summarised. The suggested
    ``min_identity`` is its ``quantile`` (1% by default) among kept units,
    rounded down to a percent and held between ``floor`` (the search's
    identity cut-off, below which nothing is seen) and ``ceiling``. Classes
    with fewer than ``min_units`` such units get no suggestion: the default
    applies. ``other_closer`` counts units nearer to another genus by any
    amount: their share bounds what any margin can set aside.
    """
    groups = defaultdict(list)
    for row in rows:
        groups[class_of.get(row["unit_id"], "") or "unknown"].append(row)
    out = []
    for name in sorted(groups):
        group = groups[name]
        screened = [r for r in group if int(r.get("genus_units") or 0) > 1]
        congener = [float(r["best_congener_identity"]) for r in screened
                    if r.get("best_congener_identity") and r["decision"] in ("keep", "kept_whole_genus")]
        closer = sum(1 for r in screened if r.get("best_other_identity")
                     and float(r["best_other_identity"]) > float(r.get("best_congener_identity") or 0))
        aside = sum(1 for r in screened if r["decision"] in SET_ASIDE)
        shares = [float(r["nonsynonymous_share"]) for r in group if r.get("nonsynonymous_share")]
        suggestion = ""
        if len(congener) >= min_units:
            value = math.floor(_quantile(congener, quantile) * 100) / 100
            suggestion = f"{min(max(value, floor), ceiling):.2f}"
        out.append({
            "class": name, "n_units": len(group), "n_genera": len({r["genus"] for r in group}),
            "screened_units": len(screened),
            "congener_p01": f"{_quantile(congener, 0.01):.4f}" if congener else "",
            "congener_p05": f"{_quantile(congener, 0.05):.4f}" if congener else "",
            "congener_p50": f"{_quantile(congener, 0.50):.4f}" if congener else "",
            "other_closer": closer, "other_closer_share": f"{closer / len(screened):.4f}" if screened else "",
            "set_aside": aside, "set_aside_share": f"{aside / len(screened):.4f}" if screened else "",
            "nonsynonymous_p95": f"{_quantile(shares, 0.95):.3f}" if shares else "",
            "numt_like": sum(1 for r in group if r.get("translation") == "numt_like"),
            "suggested_min_identity": suggestion,
        })
    return out
