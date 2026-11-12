"""Re-anchoring genera that placement concordance flags, before dating (#50).

A genus flagged by placement concordance (#49) sits in the wrong clade of
the synthetic tree, mostly because OTT merged the wrong records into it (as
for Onchidella). Its stem age, the ages above it, the outgroups of its
neighbours and its own outgroup all follow from that placement, so the
correction is made to the topology, before any of those are derived.

The genus is moved to the most specific taxon of its reference lineage
(COL, else its records) that OTT holds and that agrees with the rest of
that lineage, and is attached in the synthetic tree at that taxon's node;
where the taxon is not a node of the synthesis, at the MRCA of its members
there, and failing that at its nearest ancestor in the synthesis. The
genus's own synthesis node and species are left out of that search, since
they are what is misplaced. TimeTree queries and outgroups are then built
from the corrected backbone as for any other genus.

A curated table can force a correction (``reanchor``) or prevent one
(``keep``). Every decision is recorded in ``genera.tsv`` and
``concordance.tsv``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np

from .seqio import read_tsv

MODES = ("apply", "report", "off")
OVERRIDE_ACTIONS = ("reanchor", "keep")
OVERRIDE_COLUMNS = ("genus", "action", "target", "source", "note")
# genera.tsv columns (#50)
COLUMNS = ("placement_action", "original_status", "original_anchor", "reanchor_target", "reanchor_how",
           "reanchor_trigger", "reanchor_source", "reanchor_note")


def read_overrides(path: str | Path | None) -> dict[str, dict[str, str]]:
    """The curated table: per genus an ``action`` (``reanchor`` or ``keep``), for ``reanchor`` a ``target``
    (``rank:name``, e.g. ``family:Onchidiidae``, or ``ott:<id>``), and who decided and why (``source``,
    ``note``)."""
    if not path:
        return {}
    out = {}
    for row in read_tsv(path):
        genus = (row.get("genus") or "").strip()
        if not genus:
            continue
        action = (row.get("action") or "").strip()
        if action not in OVERRIDE_ACTIONS:
            raise ValueError(f"{path}: action for {genus} must be one of {OVERRIDE_ACTIONS}, not {action!r}")
        if action == "reanchor" and not (row.get("target") or "").strip():
            raise ValueError(f"{path}: a reanchor of {genus} needs a target (rank:name or ott:<id>)")
        out[genus] = {k: (row.get(k) or "").strip() for k in OVERRIDE_COLUMNS}
    return out


def _lineage_names(taxonomy, uid: int, limit: int = 200) -> set[str]:
    names, node = set(), int(uid)
    while node >= 0 and limit:
        rank_name = taxonomy.rank_name.get(node)
        if rank_name:
            names.add(rank_name[1])
        node = int(taxonomy.parent[node]) if node < len(taxonomy.parent) else -1
        limit -= 1
    return names


def find_target(taxonomy, names: Mapping[tuple[str, str], list[int]], reference: Mapping[str, str],
                target: str = "") -> tuple[int | None, str, str]:
    """The OTT taxon to move a genus to: its id, ``rank:name``, and how it was chosen.

    With ``target`` (from the curated table) that taxon. Otherwise the
    reference family, else order, else class, if OTT holds exactly one taxon
    of that name whose lineage contains another name of the reference
    lineage (so a homonym family elsewhere is not taken), or exactly one at
    all where none does.
    """
    if target:
        if target.startswith("ott:"):
            uid = int(target[4:])
            rank_name = taxonomy.rank_name.get(uid)
            return uid, f"{rank_name[0]}:{rank_name[1]}" if rank_name else target, "override"
        rank, _, name = target.partition(":")
        found = names.get((rank, name), [])
        return (found[0], target, "override") if len(found) == 1 else (None, target, "override_not_found")
    ranks = ("family", "order", "class")
    for i, rank in enumerate(ranks):
        name = (reference.get(rank) or "").strip()
        candidates = names.get((rank, name), []) if name else []
        if not candidates:
            continue
        higher = {(reference.get(r) or "").strip() for r in ranks[i + 1:] + ("phylum",)} - {""}
        agreeing = [uid for uid in candidates if higher & _lineage_names(taxonomy, uid)]
        if len(agreeing) == 1:
            return agreeing[0], f"{rank}:{name}", "reference_lineage"
        if not agreeing and len(candidates) == 1:
            return candidates[0], f"{rank}:{name}", "reference_name"
        return None, f"{rank}:{name}", "ambiguous_target"
    return None, "", "no_target"


def synthesis_node(synth, taxonomy, uid: int, exclude: int | None = None, max_steps: int = 15) -> tuple[int | None, str]:
    """Node of the synthetic tree standing for OTT taxon ``uid``, and how it was found.

    ``taxon_in_synthesis`` where the taxon is a node; ``mrca_of_members``,
    the MRCA of the synthesis nodes of OTT taxa within it (leaving out those
    within ``exclude``, the misplaced genus); ``ancestor_in_synthesis``, the
    nearest enclosing taxon that is a node.
    """
    if f"ott{uid}" in synth.index:
        return synth.index[f"ott{uid}"], "taxon_in_synthesis"
    nodes, uids = _synth_ott_nodes(synth)
    inside = np.zeros(len(uids), dtype=bool)
    excluded = np.zeros(len(uids), dtype=bool)
    current = uids.copy()
    parent = taxonomy.parent
    for _ in range(max_steps):
        inside |= current == uid
        if exclude is not None:
            excluded |= current == exclude
        valid = (current >= 0) & (current < len(parent))
        current = np.where(valid, parent[np.clip(current, 0, len(parent) - 1)], -1)
    members = nodes[inside & ~excluded].tolist()
    if members:
        return (synth.mrca(members) if len(members) > 1 else int(synth.parent[members[0]])), "mrca_of_members"
    up = int(taxonomy.parent[uid])
    while up >= 0 and f"ott{up}" not in synth.index:
        up = int(taxonomy.parent[up])
    return (synth.index[f"ott{up}"], "ancestor_in_synthesis") if up >= 0 else (None, "not_in_synthesis")


def _synth_ott_nodes(synth) -> tuple[np.ndarray, np.ndarray]:
    if getattr(synth, "_ott_nodes", None) is None:
        pairs = [(i, int(label[3:])) for i, label in enumerate(synth.label)
                 if label.startswith("ott") and label[3:].isdigit()]
        synth._ott_nodes = (np.array([p[0] for p in pairs], dtype=np.int64),
                            np.array([p[1] for p in pairs], dtype=np.int64))
    return synth._ott_nodes


def decide(concordance: list[dict], overrides: Mapping[str, dict], min_severity: int,
           genera_in_ott: Mapping[str, int]) -> dict[str, dict[str, str]]:
    """Which genera to move, and why: flagged by concordance with ``severity`` at least ``min_severity``,
    or listed for ``reanchor`` in the curated table; ``keep`` in the table wins over a flag."""
    out = {}
    by_genus = {r["genus"]: r for r in concordance}
    for genus in genera_in_ott:
        row = by_genus.get(genus) or {}
        override = overrides.get(genus)
        flagged = row.get("concordance_flag") == "flag" and int(row.get("severity") or 0) >= min_severity
        trigger = (f"flags={row.get('n_flags', '')}/{row.get('n_evaluable', '')};severity={row.get('severity', '')};"
                   f"lineage={row.get('lineage_distance', '')};merged={row.get('merged_id_distance', '')};"
                   f"neighbours={row.get('neighbour_concordance', '')}") if row else ""
        if override and override["action"] == "keep":
            if flagged:
                out[genus] = {"decision": "kept_override", "trigger": trigger, "override": override}
        elif override and override["action"] == "reanchor":
            out[genus] = {"decision": "reanchor_override", "trigger": trigger, "override": override}
        elif flagged:
            out[genus] = {"decision": "reanchor_auto", "trigger": trigger, "override": None}
    return out


def lineage_from(taxonomy, uid: int) -> list[tuple[int, int | None]]:
    """(OTT id, NCBI id) of a taxon and every taxon above it, as :func:`timegraft.otol_static.lineage`
    gives for a genus (which starts above the genus)."""
    out, node = [], int(uid)
    while node >= 0 and len(out) < 100:
        out.append((node, taxonomy.ncbi_of.get(node)))
        node = int(taxonomy.parent[node])
    return out


def plan(synth, taxonomy, names, matches: Mapping[str, dict], anchors: Mapping[str, int], concordance: list[dict],
         reference: Mapping[str, Mapping[str, str]], overrides: Mapping[str, dict], mode: str,
         min_severity: int) -> tuple[dict[str, dict], dict[str, tuple[int, int]]]:
    """Decide and locate every correction.

    Returns the genera.tsv columns (``COLUMNS``) per genus decided on, and,
    with ``mode`` ``apply``, the genera to move with their new synthesis node
    and target OTT id. ``placement_action`` is ``reanchored``,
    ``would_reanchor`` (mode ``report``), ``kept_override`` or
    ``flagged_no_target`` (no taxon of the reference lineage could be found).
    """
    if mode not in MODES:
        raise ValueError(f"reanchor.mode must be one of {MODES}, not {mode!r}")
    if mode == "off":
        return {}, {}
    decisions = decide(concordance, overrides, min_severity,
                       {g: m["ott_id"] for g, m in matches.items() if m.get("ott_id")})
    provenance, moved = {}, {}
    for genus, d in sorted(decisions.items()):
        row = {c: "" for c in COLUMNS}
        row.update(original_status=matches[genus]["status"], reanchor_trigger=d["trigger"],
                   original_anchor=synth.label[anchors[genus]] if genus in anchors else "")
        if d["override"]:
            row.update(reanchor_source=f"override:{d['override'].get('source') or 'curated'}",
                       reanchor_note=d["override"].get("note", ""))
        else:
            row["reanchor_source"] = "auto"
        if d["decision"] == "kept_override":
            row["placement_action"] = "kept_override"
            provenance[genus] = row
            continue
        target = d["override"]["target"] if d["override"] else ""
        uid, taxon, how = find_target(taxonomy, names, reference.get(genus) or {}, target)
        node, where = (synthesis_node(synth, taxonomy, uid, exclude=matches[genus]["ott_id"])
                       if uid is not None else (None, how))
        row.update(reanchor_target=f"{taxon} (ott{uid})" if uid is not None else taxon,
                   reanchor_how=f"{how};{where}" if uid is not None else how)
        if node is None:
            row["placement_action"] = "flagged_no_target"
        elif mode == "report":
            row["placement_action"] = "would_reanchor"
        else:
            row["placement_action"] = "reanchored"
            moved[genus] = (node, uid)
        provenance[genus] = row
    return provenance, moved
