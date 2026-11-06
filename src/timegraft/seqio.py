"""Reading and writing the tabular and sequence formats used by the pipeline."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Iterator

# Record table columns. Extra columns are carried along untouched.
RECORD_COLUMNS = ("record_id", "unit_id", "genus", "sequence")


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: str | Path, rows: Iterable[dict], columns: Iterable[str]) -> None:
    columns = list(columns)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_records(path: str | Path) -> list[dict[str, str]]:
    """Read a record table and check that the required columns are present."""
    rows = read_tsv(path)
    if rows:
        missing = [c for c in RECORD_COLUMNS if c not in rows[0]]
        if missing:
            raise ValueError(f"{path}: missing required columns {missing}")
    for row in rows:
        row["sequence"] = row["sequence"].strip().upper().replace("U", "T")
        row["genus"] = row["genus"].strip()
        row["unit_id"] = row["unit_id"].strip()
    return rows


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def write_fasta(path: str | Path, items: Iterable[tuple[str, str]]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        for name, seq in items:
            handle.write(f">{name}\n{seq}\n")
