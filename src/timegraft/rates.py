"""Substitution rates per taxon, for crown ages taken from the barcodes.

One rate for all taxa is the weakest assumption in rate-based crown ages
(issue #18): mitochondrial rates differ between classes (fishes are slower
than most invertebrates, sharks slower still). A rate is therefore resolved
per genus from its taxonomy, from the genus upwards, and falls back on the
marker's default rate.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping

from .seqio import read_tsv

RANKS = ("genus", "family", "order", "class", "phylum")


def genus_taxonomy(records: Iterable[Mapping[str, str]]) -> dict[str, dict[str, str]]:
    """Family, order, class and phylum per genus: the most common value among its records."""
    counts: dict[str, dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
    for record in records:
        genus = (record.get("genus") or "").strip()
        if not genus:
            continue
        for rank in RANKS[1:]:
            name = (record.get(rank) or "").strip()
            if name:
                counts[genus][rank][name] += 1
    return {
        genus: {"genus": genus, **{rank: (c[rank].most_common(1)[0][0] if c[rank] else "") for rank in RANKS[1:]}}
        for genus, c in counts.items()
    }


def read_rate_table(path: str | Path) -> dict[str, dict[str, float]]:
    """Rates from a table with columns ``rank``, ``name`` and ``rate`` (as written by timegraft.calibrate).

    Rows of rank ``all`` are a pooled value for information and are skipped.
    """
    table: dict[str, dict[str, float]] = defaultdict(dict)
    for row in read_tsv(path):
        rank = row["rank"].strip()
        if rank == "all":
            continue
        if rank not in RANKS:
            raise ValueError(f"{path}: unknown rank {rank!r}; use {list(RANKS)}")
        table[rank][row["name"].strip()] = _positive(float(row["rate"]))
    return dict(table)


def _positive(rate: float) -> float:
    if not rate > 0:
        raise ValueError(f"a substitution rate must be positive, not {rate}")
    return rate


class Rates:
    """Resolve the substitution rate of a genus from its taxonomy.

    ``overrides`` maps rank to {name: rate}, e.g. ``{"class": {"Teleostei": 0.0045}}``,
    and wins over ``table`` (a path, or a mapping of the same shape).
    """

    def __init__(
        self,
        default: float | None,
        overrides: Mapping[str, Mapping[str, float]] | None = None,
        table: str | Path | Mapping[str, Mapping[str, float]] | None = None,
    ):
        self.default = None if default is None else _positive(float(default))
        self.overrides = {
            rank: {name: _positive(float(rate)) for name, rate in (names or {}).items()}
            for rank, names in (overrides or {}).items()
        }
        unknown = set(self.overrides) - set(RANKS)
        if unknown:
            raise ValueError(f"rates for unknown ranks {sorted(unknown)}; use {list(RANKS)}")
        if table is None or isinstance(table, Mapping):
            self.table = {rank: dict(names) for rank, names in (table or {}).items()}
        else:
            self.table = read_rate_table(table)

    def resolve(self, taxa: Mapping[str, str] | None) -> tuple[float | None, str]:
        """Return (rate, source); source is ``rank:name`` or ``default``."""
        for rank in RANKS:
            name = ((taxa or {}).get(rank) or "").strip()
            if not name:
                continue
            for source in (self.overrides, self.table):
                if name in source.get(rank, {}):
                    return source[rank][name], f"{rank}:{name}"
        return self.default, "default"
