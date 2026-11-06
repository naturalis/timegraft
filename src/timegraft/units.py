"""From records to one representative sequence per reference unit (BIN or SH)."""

from __future__ import annotations

from collections import Counter, defaultdict

from . import distances

UNASSIGNED_COLUMNS = ("unit_id", "reason", "genera", "n_records", "record_id", "sequence")
REPRESENTATIVE_COLUMNS = ("unit_id", "genus", "n_records", "genus_share", "record_id", "sequence")


def medoid(records: list[dict], k: int = 8) -> dict:
    """Record with the smallest summed k-mer distance to the others.

    k-mer distance is used regardless of the marker so that medoid choice
    never depends on an alignment.
    """
    if len(records) <= 2:
        return records[0]
    seqs = [r["sequence"] for r in records]
    dist = distances.matrix(seqs, method="kmer", k=k)
    return records[int(dist.sum(axis=1).argmin())]


def representatives(
    records: list[dict],
    max_records: int = 50,
    genus_agreement: float = 0.8,
    k: int = 8,
) -> tuple[list[dict], list[dict]]:
    """Pick a medoid per unit and decide the unit's genus.

    A unit is assigned to a genus when at least ``genus_agreement`` of its
    records that name a genus agree. Units without a genus, or with
    conflicting genera, are returned separately for later placement.
    """
    by_unit: dict[str, list[dict]] = defaultdict(list)
    for record in records:
        if record["unit_id"]:
            by_unit[record["unit_id"]].append(record)

    reps, unassigned = [], []
    for unit_id in sorted(by_unit):
        group = sorted(by_unit[unit_id], key=lambda r: r["record_id"])
        genera = Counter(r["genus"] for r in group if r["genus"])
        if not genera:
            chosen = medoid(group[:max_records], k=k)
            unassigned.append({"unit_id": unit_id, "reason": "no_genus", "genera": "", "n_records": len(group),
                               "record_id": chosen["record_id"], "sequence": chosen["sequence"]})
            continue
        genus, votes = genera.most_common(1)[0]
        share = votes / sum(genera.values())
        if share < genus_agreement:
            chosen = medoid(group[:max_records], k=k)
            unassigned.append(
                {
                    "unit_id": unit_id,
                    "reason": "genus_conflict",
                    "genera": ";".join(f"{g}:{n}" for g, n in sorted(genera.items())),
                    "n_records": len(group),
                    "record_id": chosen["record_id"],
                    "sequence": chosen["sequence"],
                }
            )
            continue
        pool = [r for r in group if r["genus"] == genus][:max_records]
        chosen = medoid(pool, k=k)
        reps.append(
            {
                "unit_id": unit_id,
                "genus": genus,
                "n_records": len(group),
                "genus_share": f"{share:.3f}",
                "record_id": chosen["record_id"],
                "sequence": chosen["sequence"],
            }
        )
    return reps, unassigned
