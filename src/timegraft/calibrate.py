"""Calibrating the substitution rate per taxon from genera of known crown age (issue #18).

For a genus with a known crown age, the rate implied by its barcodes is the
mean root-to-tip distance of its rooted subtree divided by that age. Known ages
are point ages from the user's table (literature) and, failing that, TimeTree
crown ages below the stem age. A TimeTree crown age is the MRCA of the species
TimeTree knows under the genus name and is far too old for genera that are not
natural groups, so the implied rates have a tail of values that are much too
low. The rate of a taxon is therefore a trimmed median: the median of its
genera, recomputed without those more than ``window`` times above or below it,
until nothing changes.

    python -m timegraft.calibrate --provenance results/x/provenance.tsv \\
        --taxonomy results/x/genera/taxonomy.tsv --out rates.tsv

The result can be given to the workflow as ``crown_age.rate_table``.
"""

from __future__ import annotations

import argparse
import statistics
from typing import Iterable, Mapping

from .rates import RANKS
from .seqio import read_tsv, write_tsv

COLUMNS = ("rank", "name", "rate", "n", "n_trimmed", "q1", "q3", "untrimmed_median")


def implied_rates(
    provenance: Iterable[Mapping[str, str]],
    table: Mapping[str, float] | None = None,
    min_units: int = 2,
) -> dict[str, tuple[float, str]]:
    """Implied rate and its source (table or timetree) per genus with a usable known age."""
    out = {}
    for row in provenance:
        if row.get("status", "subtree") != "subtree" or not row.get("root_to_tip"):
            continue
        distance, stem = float(row["root_to_tip"]), float(row["stem_age"])
        if distance <= 0 or int(row["n_units"]) < min_units:
            continue
        genus = row["genus"]
        if table and genus in table:
            age, source = float(table[genus]), "table"
        elif row.get("timetree_crown"):
            age, source = float(row["timetree_crown"]), "timetree"
        else:
            continue
        if 0 < age < stem:
            out[genus] = (distance / age, source)
    return out


def trimmed_median(values: list[float], window: float = 3.0, max_iter: int = 100) -> tuple[float, list[float]]:
    """Median of ``values`` after dropping, repeatedly, those beyond ``window`` times the median either way."""
    kept = list(values)
    centre = statistics.median(kept)
    for _ in range(max_iter):
        inside = [v for v in values if centre / window <= v <= centre * window]
        if not inside or inside == kept:
            break
        kept, centre = inside, statistics.median(inside)
    return centre, kept


def calibrate(
    rates: Mapping[str, tuple[float, str]],
    taxonomy: Mapping[str, Mapping[str, str]],
    min_genera: int = 10,
    window: float = 3.0,
    ranks: Iterable[str] = ("class",),
) -> list[dict]:
    """Rate per taxon of the given ``ranks`` with at least ``min_genera`` genera left after trimming.

    One pooled row (rank ``all``) comes first; it is the value to use as the
    marker's default rate. Lower ranks than class have few genera each and
    mostly add noise, hence the default.
    """
    ranks = tuple(ranks)
    unknown = set(ranks) - set(RANKS[1:])
    if unknown:
        raise ValueError(f"cannot calibrate at ranks {sorted(unknown)}; use {list(RANKS[1:])}")
    groups: dict[tuple[str, str], list[float]] = {("all", "all"): [r for r, _ in rates.values()]}
    for genus, (rate, _) in rates.items():
        for rank in ranks:
            name = (taxonomy.get(genus, {}).get(rank) or "").strip()
            if name:
                groups.setdefault((rank, name), []).append(rate)
    rows = []
    for (rank, name), values in groups.items():
        if len(values) < min_genera:
            continue
        centre, kept = trimmed_median(values, window)
        if len(kept) < min_genera:
            continue
        q = statistics.quantiles(kept, n=4)
        rows.append({
            "rank": rank, "name": name, "rate": f"{centre:.5f}", "n": len(kept), "n_trimmed": len(values) - len(kept),
            "q1": f"{q[0]:.5f}", "q3": f"{q[2]:.5f}", "untrimmed_median": f"{statistics.median(values):.5f}",
        })
    order = {rank: i for i, rank in enumerate(("all",) + RANKS[::-1])}
    rows.sort(key=lambda r: (order[r["rank"]], -r["n"], r["name"]))
    return rows


def read_crown_table(path) -> dict[str, dict[str, float | None]]:
    """The user's crown age table: ``genus`` and any of ``crown_age``, ``crown_age_min``, ``crown_age_max``."""
    out = {}
    for row in read_tsv(path):
        entry = {key: (float(row[col]) if (row.get(col) or "").strip() else None)
                 for key, col in (("age", "crown_age"), ("min", "crown_age_min"), ("max", "crown_age_max"))}
        if any(v is not None for v in entry.values()):
            out[row["genus"].strip()] = entry
    return out


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--provenance", required=True, help="provenance.tsv of a run (needs the root_to_tip column)")
    parser.add_argument("--taxonomy", required=True, help="genera/taxonomy.tsv of the run: genus, family, order, class, phylum")
    parser.add_argument("--crown-ages", help="table of known crown ages (genus, crown_age); used before TimeTree")
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-genera", type=int, default=10)
    parser.add_argument("--window", type=float, default=3.0)
    parser.add_argument("--min-units", type=int, default=2, help="leave out genera with fewer units")
    parser.add_argument("--ranks", nargs="+", default=["class"], help="ranks to calibrate at (default: class)")
    args = parser.parse_args(argv)
    table = {}
    if args.crown_ages:
        table = {g: e["age"] for g, e in read_crown_table(args.crown_ages).items() if e["age"] is not None}
    rates = implied_rates(read_tsv(args.provenance), table, args.min_units)
    taxonomy = {row["genus"]: row for row in read_tsv(args.taxonomy)}
    rows = calibrate(rates, taxonomy, args.min_genera, args.window, args.ranks)
    write_tsv(args.out, rows, COLUMNS)
    print(f"{len(rates)} genera with a known crown age; {len(rows)} rates written to {args.out}")


if __name__ == "__main__":
    main()
