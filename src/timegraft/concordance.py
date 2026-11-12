"""Placement concordance: is a genus where OToL puts it? Report only (#49).

Three checks, each independent of the others and of OToL's placement:

- **Reference lineage.** Where in OTT is the family of the genus according
  to a reference classification (COL where it resolved the genus, otherwise
  the most common family among its records)? The distance is how far up
  OTT one has to go from the genus to reach that family.
- **Merged ids.** Where OTT merged more than one NCBI id into the genus, how
  far apart are those ids in NCBI Taxonomy?
- **Barcode neighbours.** Do the genera that the genus's units come closest
  to (outside the genus) share a phylum with its sisters in the backbone?

A distance is the rank level at which two lineages part (``RANK_LEVEL``):
0 where they agree, 1 when they part at family level, up to 5 for
different phyla and 6 above that. Each check alone has false positives
(classification differences, synonyms, slow mitochondrial evolution);
a genus is flagged when two of the checks that could be made agree. The
holistic analysis found that combination isolates the one real
misplacement of the Iberian run (Onchidella) with nothing else.

Nothing here changes the backbone; correction is #50.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Iterable, Mapping, Sequence

RANK_LEVEL = {
    "tribe": 1, "subtribe": 1, "supertribe": 1, "subfamily": 1, "family": 1,
    "superfamily": 2,
    "parvorder": 3, "infraorder": 3, "suborder": 3, "order": 3, "superorder": 3, "cohort": 3, "subcohort": 3,
    "infraclass": 4, "subclass": 4, "class": 4, "superclass": 4,
    "infraphylum": 5, "subphylum": 5, "phylum": 5, "division": 5, "subdivision": 5, "superphylum": 5,
    "infrakingdom": 6, "subkingdom": 6, "kingdom": 6, "superkingdom": 6, "domain": 6, "realm": 6,
}
LEVEL_NAME = {0: "", 1: "family", 2: "superfamily", 3: "order", 4: "class", 5: "phylum", 6: "kingdom"}
LEVELS = {v: k for k, v in LEVEL_NAME.items() if v}

COLUMNS = (
    "rank", "genus", "status", "ott_id", "ott_family", "ott_class", "ott_phylum",
    # reference lineage
    "reference_source", "reference_taxon", "lineage_mrca", "lineage_split", "lineage_distance", "lineage_flag",
    # merged ids
    "merged_ncbi_ids", "merged_mrca", "merged_split", "merged_id_distance", "merged_flag",
    # barcode neighbours
    "neighbour_hits", "neighbour_phyla", "sister_phyla", "neighbour_concordance", "neighbour_flag",
    # aggregate
    "n_evaluable", "n_flags", "severity", "concordance_flag",
)


def split_level(a: Sequence[tuple[object, str]], b: Sequence[tuple[object, str]]) -> tuple[int, object]:
    """Level at which two lineages part, and the key of their last shared taxon.

    Each lineage is a list of (key, rank) from the taxon itself up to the
    root; keys identify taxa (ids, or names within one classification). The
    level is the highest ``RANK_LEVEL`` among the taxa below the shared part
    on either side, so a split between two phyla at an unranked node is still
    level 5. 0 if one lineage contains the other's taxon.
    """
    keys_a, keys_b = {k for k, _ in a}, {k for k, _ in b}
    if a and a[0][0] in keys_b:
        return 0, a[0][0]
    if b and b[0][0] in keys_a:
        return 0, b[0][0]
    shared = next((k for k, _ in a if k in keys_b), None)
    below = []
    for lineage in (a, b):
        for key, rank in lineage:
            if key == shared:
                break
            below.append(rank)
    return max((RANK_LEVEL.get(r, 0) for r in below), default=0), shared


def ott_lineage(taxonomy, uid: int, limit: int = 200) -> list[tuple[int, str]]:
    """(OTT id, rank) of a taxon and all that enclose it, nearest first."""
    out, node = [], int(uid)
    while node >= 0 and len(out) < limit:
        out.append((node, (taxonomy.rank_name.get(node) or ("", ""))[0]))
        node = int(taxonomy.parent[node]) if node < len(taxonomy.parent) else -1
    return out


def _ott_name(taxonomy, uid) -> str:
    if uid is None:
        return ""
    rank_name = taxonomy.rank_name.get(int(uid))
    return rank_name[1] if rank_name else f"ott{uid}"


def name_index(taxonomy, ranks: Iterable[str] = ("family", "order", "class", "phylum")) -> dict[tuple[str, str], list[int]]:
    """OTT ids of the higher taxa of the given ranks, by (rank, name), within the matching context."""
    wanted = set(ranks)
    out: dict[tuple[str, str], list[int]] = defaultdict(list)
    for uid, (rank, name) in taxonomy.rank_name.items():
        if rank in wanted and taxonomy._within_context(uid):
            out[(rank, name)].append(uid)
    return out


def lineage_distance(taxonomy, names: Mapping[tuple[str, str], list[int]], genus_uid: int,
                     reference: Mapping[str, str]) -> dict[str, str]:
    """How far the OTT placement of a genus is from the family (else order, class) a reference gives it.

    Where the reference name is held by more than one OTT taxon, the closest
    is used, so homonyms never make a genus look misplaced.
    """
    out = {"reference_taxon": "", "lineage_mrca": "", "lineage_split": "", "lineage_distance": ""}
    genus = ott_lineage(taxonomy, genus_uid)
    for rank in ("family", "order", "class"):
        name = (reference.get(rank) or "").strip()
        if not name or (rank, name) not in names:
            continue
        best = min((split_level(genus, ott_lineage(taxonomy, uid)) for uid in names[(rank, name)]),
                   key=lambda x: x[0])
        out.update(reference_taxon=f"{rank}:{name}", lineage_mrca=_ott_name(taxonomy, best[1]),
                   lineage_split=LEVEL_NAME[best[0]], lineage_distance=str(best[0]))
        break
    return out


def classification_lineage(classification: Sequence[Mapping]) -> list[tuple[str, str]]:
    """(name, rank) from a ChecklistBank classification (root first), nearest first."""
    return [(c.get("name"), c.get("rank") or "") for c in reversed(list(classification)) if c.get("name")]


def merged_id_distance(classifications: Mapping[str, Sequence[Mapping]], ids: Sequence[str]) -> dict[str, str]:
    """How far apart in NCBI Taxonomy the NCBI ids OTT merged into one genus are (the widest pair)."""
    out = {"merged_mrca": "", "merged_split": "", "merged_id_distance": ""}
    known = [i for i in ids if classifications.get(i)]
    if len(ids) < 2 or len(known) < 2:
        return out
    worst = (-1, None)
    for x in range(len(known)):
        for y in range(x + 1, len(known)):
            a = [(f"id:{known[x]}", "genus")] + classification_lineage(classifications[known[x]])
            b = [(f"id:{known[y]}", "genus")] + classification_lineage(classifications[known[y]])
            level, shared = split_level(a, b)
            if level > worst[0]:
                worst = (level, shared)
    out.update(merged_mrca=worst[1] or "", merged_split=LEVEL_NAME[worst[0]], merged_id_distance=str(worst[0]))
    return out


def backbone_sisters(tree) -> dict[str, list[str]]:
    """Per genus tip of the backbone, the genera of its sister clade(s)."""
    out = {}
    for leaf in tree.leaf_node_iter():
        parent = leaf.parent_node
        if parent is None:
            continue
        out[leaf.taxon.label] = [l.taxon.label for c in parent.child_node_iter() if c is not leaf
                                 for l in c.leaf_iter()]
    return out


def neighbour_concordance(genus: str, hits: Sequence[str], sisters: Sequence[str], phylum: Mapping[str, str],
                          min_hits: int = 5) -> dict[str, str]:
    """Share of the genus's cross-genus barcode hits whose phylum is among those of its backbone sisters.

    ``hits`` are the nearest genera outside the genus, one per unit that has
    one. Not evaluable with fewer than ``min_hits`` hits of known phylum, or
    no sister of known phylum.
    """
    counts = Counter(phylum.get(g, "") for g in hits)
    counts.pop("", None)
    sister_phyla = sorted({phylum[g] for g in sisters if phylum.get(g)})
    out = {"neighbour_hits": str(sum(counts.values())),
           "neighbour_phyla": ";".join(f"{p}:{n}" for p, n in counts.most_common()),
           "sister_phyla": ";".join(sister_phyla[:10]) + (";..." if len(sister_phyla) > 10 else ""),
           "neighbour_concordance": ""}
    total = sum(counts.values())
    if total >= min_hits and sister_phyla:
        out["neighbour_concordance"] = f"{sum(n for p, n in counts.items() if p in sister_phyla) / total:.3f}"
    return out


def aggregate(row: dict, flag_level: int, max_share: float, min_flags: int = 2) -> dict:
    """Component flags, the number of checks made and failed, severity, and the overall flag.

    ``concordance_flag`` is ``flag`` when at least ``min_flags`` checks fail,
    ``review`` when one does, empty when none does, and ``not_evaluable``
    when no check could be made.
    """
    flags = {
        "lineage_flag": row.get("lineage_distance", "") != "" and int(row["lineage_distance"]) >= flag_level,
        "merged_flag": row.get("merged_id_distance", "") != "" and int(row["merged_id_distance"]) >= flag_level,
        "neighbour_flag": row.get("neighbour_concordance", "") != "" and float(row["neighbour_concordance"]) <= max_share,
    }
    made = sum(1 for k in ("lineage_distance", "merged_id_distance", "neighbour_concordance") if row.get(k, "") != "")
    failed = sum(flags.values())
    levels = [int(row[k]) for k in ("lineage_distance", "merged_id_distance") if row.get(k, "") != ""]
    row.update({k: "yes" if v else "" for k, v in flags.items()})
    row.update(n_evaluable=str(made), n_flags=str(failed), severity=str(max(levels, default=0)),
               concordance_flag="not_evaluable" if not made else
               "flag" if failed >= min_flags else "review" if failed else "")
    return row


def rank_rows(rows: list[dict]) -> list[dict]:
    """Sort for manual review: most failed checks first, then the widest split; number them."""
    rows = sorted(rows, key=lambda r: (-int(r["n_flags"] or 0), -int(r["severity"] or 0),
                                       r["neighbour_concordance"] or "9", r["genus"]))
    for i, row in enumerate(rows, 1):
        row["rank"] = str(i)
    return rows
