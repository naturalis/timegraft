"""Turn a BOLD data package (BCDM TSV) into a timegraft record table.

Streams the file, so a full BOLD dump can be filtered without loading it.
Usable from the workflow (rule ``bcdm_records``) or on its own:

    python -m timegraft.bcdm BOLD_Public.tsv records.tsv \\
        --taxon phylum=Annelida,Mollusca,Arthropoda --country Portugal,Spain
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

RANKS = ("phylum", "class", "order", "family", "subfamily", "genus", "species")
OUT_COLUMNS = ("record_id", "unit_id", "genus", "family", "order", "class", "phylum", "species", "country", "sequence")

# BCDM column names, with variants seen in older exports
COLUMN_OPTIONS = {
    "record_id": ["processid"],
    "unit_id": ["bin_uri", "bin"],
    "sequence": ["nuc", "nucleotides"],
    "marker": ["marker_code", "markercode"],
    "country": ["country/ocean", "country", "country_ocean"],
}


def _resolve(fieldnames: list[str]) -> dict[str, str | None]:
    lower = {name.lower(): name for name in fieldnames}
    out = {}
    for key, options in COLUMN_OPTIONS.items():
        out[key] = next((lower[o] for o in options if o in lower), None)
    for rank in RANKS:
        out[rank] = lower.get(rank)
    for key in ("record_id", "unit_id", "sequence"):
        if out[key] is None:
            raise ValueError(f"column for {key} not found; looked for {COLUMN_OPTIONS[key]}")
    return out


def _clean(value: str | None) -> str:
    value = (value or "").strip()
    return "" if value.lower() in ("none", "na", "null") else value


def convert(
    source: str | Path,
    target: str | Path,
    marker: str = "COI-5P",
    taxa: dict[str, set[str]] | None = None,
    countries: set[str] | None = None,
) -> dict[str, int]:
    """Filter ``source`` and write ``target``; returns counts of kept and skipped records."""
    taxa = taxa or {}
    counts = {"read": 0, "kept": 0, "other_marker": 0, "no_unit": 0, "no_sequence": 0, "taxon": 0, "country": 0}
    csv.field_size_limit(sys.maxsize)
    Path(target).parent.mkdir(parents=True, exist_ok=True)
    with open(source, newline="", encoding="utf-8", errors="replace") as src, open(target, "w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t", quoting=csv.QUOTE_NONE)
        cols = _resolve(reader.fieldnames or [])
        writer = csv.DictWriter(dst, fieldnames=OUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            counts["read"] += 1
            if marker and cols["marker"] and _clean(row.get(cols["marker"])) != marker:
                counts["other_marker"] += 1
                continue
            unit = _clean(row.get(cols["unit_id"]))
            if not unit:
                counts["no_unit"] += 1
                continue
            seq = _clean(row.get(cols["sequence"])).upper()
            if not seq:
                counts["no_sequence"] += 1
                continue
            values = {rank: _clean(row.get(cols[rank])) if cols[rank] else "" for rank in RANKS}
            if any(values.get(rank) not in allowed for rank, allowed in taxa.items()):
                counts["taxon"] += 1
                continue
            country = _clean(row.get(cols["country"])) if cols["country"] else ""
            if countries and country not in countries:
                counts["country"] += 1
                continue
            writer.writerow({"record_id": _clean(row.get(cols["record_id"])), "unit_id": unit,
                             "country": country, "sequence": seq,
                             **{k: values[k] for k in ("genus", "family", "order", "class", "phylum", "species")}})
            counts["kept"] += 1
    return counts


def _parse_taxa(items: list[str]) -> dict[str, set[str]]:
    taxa: dict[str, set[str]] = {}
    for item in items:
        rank, _, names = item.partition("=")
        if rank not in RANKS or not names:
            raise SystemExit(f"--taxon expects rank=Name1,Name2 with rank one of {RANKS}")
        taxa.setdefault(rank, set()).update(n.strip() for n in names.split(",") if n.strip())
    return taxa


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("source", help="BCDM TSV from a BOLD data package")
    parser.add_argument("target", help="record table to write")
    parser.add_argument("--marker", default="COI-5P", help="marker code to keep (empty for all)")
    parser.add_argument("--taxon", action="append", default=[], help="rank=Name1,Name2; repeatable")
    parser.add_argument("--country", default="", help="comma-separated country/ocean values to keep")
    args = parser.parse_args(argv)
    countries = {c.strip() for c in args.country.split(",") if c.strip()} or None
    counts = convert(args.source, args.target, args.marker, _parse_taxa(args.taxon), countries)
    print("\t".join(f"{k}={v}" for k, v in counts.items()))


if __name__ == "__main__":
    main()
