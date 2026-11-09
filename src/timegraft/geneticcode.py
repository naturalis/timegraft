"""Genetic codes: NCBI translation tables, translation, and a code per record.

Taxonomically broad data sets span several mitochondrial genetic codes, so the
code is resolved per record from its taxonomy rather than set once per run.
For each record, from the most specific rank upwards (genus, subfamily, family,
order, class, phylum, kingdom), the first name found in the user's overrides
or in the mapping shipped with the package decides the code; at a given rank,
a user override wins over the shipped mapping. A ``genetic_code`` column in the
record table overrides all of this, and records that match nothing get the
default code.

The shipped mapping (``data/genetic_codes.tsv``) can be regenerated from the
NCBI taxonomy dump:

    python -m timegraft.geneticcode taxdump.tar.gz src/timegraft/data/genetic_codes.tsv

Names are matched regardless of rank, because BOLD and NCBI place some taxa at
different ranks (e.g. Elasmobranchii is a class in BOLD and a subclass in NCBI).
A name whose descendants use more than one code (e.g. Chordata) is marked
``mixed`` and is skipped, so that a more specific name, or the default, decides.
"""

from __future__ import annotations

import argparse
import csv
import io
import tarfile
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Iterable, Mapping

# Amino acids per codon, in NCBI order (first, second, third base each TCAG);
# '*' is a stop. Generated from Biopython's Bio.Data.CodonTable, which follows
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi; tests check them.
TABLES: dict[int, str] = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Standard
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",  # Vertebrate Mitochondrial
    3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Yeast Mitochondrial
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Mold, Protozoan, Coelenterate Mito.
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",  # Invertebrate Mitochondrial
    9: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",  # Echinoderm and Flatworm Mito.
    13: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",  # Ascidian Mitochondrial
    14: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",  # Alternative Flatworm Mito.
    21: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",  # Trematode Mitochondrial
    24: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",  # Pterobranchia Mitochondrial
    33: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",  # Cephalodiscidae Mitochondrial
}

_BASES = "TCAG"
_CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
_MAPS: dict[int, dict[str, str]] = {code: dict(zip(_CODONS, aas)) for code, aas in TABLES.items()}

DEFAULT_TABLE = Path(__file__).parent / "data" / "genetic_codes.tsv"

# Record columns searched for a taxon name, most specific first.
RANKS = ("genus", "subfamily", "family", "order", "class", "phylum", "kingdom")

MIXED = "mixed"


def _check(code: int) -> int:
    if code not in TABLES:
        raise ValueError(f"genetic code {code} not supported; choose from {sorted(TABLES)}")
    return code


def _degap(seq: str) -> str:
    return seq.upper().replace("-", "").replace(".", "").replace("U", "T")


def stop_codons(code: int) -> frozenset[str]:
    return frozenset(codon for codon, aa in _MAPS[_check(code)].items() if aa == "*")


def translate(seq: str, code: int, frame: int = 0) -> str:
    """Translate one forward frame. Codons with ambiguous bases become 'X'."""
    table = _MAPS[_check(code)]
    seq = _degap(seq)
    return "".join(table.get(seq[i : i + 3], "X") for i in range(frame, len(seq) - 2, 3))


def frame_stops(seq: str, code: int) -> list[int]:
    """Number of stop codons in each of the three forward frames."""
    stops = stop_codons(code)
    seq = _degap(seq)
    return [sum(seq[i : i + 3] in stops for i in range(frame, len(seq) - 2, 3)) for frame in range(3)]


def best_frame(seq: str, code: int) -> tuple[int, int]:
    """The forward frame with the fewest stop codons, and that number."""
    counts = frame_stops(seq, code)
    frame = min(range(3), key=counts.__getitem__)
    return frame, counts[frame]


def has_open_frame(seq: str, code: int) -> bool:
    """True if any of the three forward frames has no stop codon."""
    return best_frame(seq, code)[1] == 0


# Code per record -------------------------------------------------------------


def read_table(path: str | Path) -> dict[str, int | None]:
    """Read a name-to-code mapping; ``mixed`` names map to None."""
    out: dict[str, int | None] = {}
    with open(path, newline="") as handle:
        rows = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        for row in rows:
            value = row["code"].strip()
            out[row["name"].strip()] = None if value == MIXED else _check(int(value))
    return out


class GeneticCodes:
    """Resolve the genetic code of a record from its taxonomy.

    ``overrides`` maps rank to {name: code}, e.g. ``{"class": {"Ascidiacea": 13}}``.
    ``table`` is the path of a name-to-code mapping; None for the shipped one,
    False for none at all.
    """

    def __init__(
        self,
        default: int = 5,
        overrides: Mapping[str, Mapping[str, int]] | None = None,
        table: str | Path | None | bool = None,
    ):
        self.default = _check(int(default))
        self.overrides = {
            rank: {name: _check(int(code)) for name, code in (names or {}).items()}
            for rank, names in (overrides or {}).items()
        }
        unknown = set(self.overrides) - set(RANKS)
        if unknown:
            raise ValueError(f"genetic code overrides for unknown ranks {sorted(unknown)}; use {list(RANKS)}")
        if table is False:
            self.table: dict[str, int | None] = {}
        else:
            self.table = read_table(table or DEFAULT_TABLE)

    def resolve(self, record: Mapping[str, str]) -> tuple[int, str, str]:
        """Return (code, source, taxon); source is column, config, table or default."""
        given = (record.get("genetic_code") or "").strip()
        if given:
            return _check(int(given)), "column", ""
        for rank in RANKS:
            name = (record.get(rank) or "").strip()
            if not name:
                continue
            if name in self.overrides.get(rank, {}):
                return self.overrides[rank][name], "config", name
            if self.table.get(name) is not None:
                return self.table[name], "table", name
        return self.default, "default", ""


# Building the mapping from the NCBI taxonomy dump ----------------------------

# Ranks written to the mapping. Ranks between these are included because
# BOLD and NCBI do not always agree on rank (e.g. Elasmobranchii, Actinopterygii).
EMIT_RANKS = frozenset({
    "kingdom", "subkingdom", "phylum", "subphylum", "superclass", "class", "subclass",
    "infraclass", "superorder", "order",
})
# Ranks whose codes count when deciding if a name is mixed.
CODE_RANKS = EMIT_RANKS | {"suborder", "infraorder", "superfamily", "family", "subfamily", "genus"}
NAME_CLASSES = frozenset({"scientific name", "synonym", "equivalent name"})


def _dmp(handle: Iterable[str]) -> Iterable[list[str]]:
    for line in handle:
        yield line.rstrip("\n").rstrip("\t|").split("\t|\t")


def _open_taxdump(path: Path, member: str):
    if path.is_dir():
        return open(path / member, encoding="utf-8")
    archive = tarfile.open(path)
    return io.TextIOWrapper(archive.extractfile(member), encoding="utf-8")


def build_table(taxdump: str | Path, root: int = 2759) -> list[dict[str, str]]:
    """Name-to-mitochondrial-code rows for taxa under ``root`` (default Eukaryota).

    ``taxdump`` is NCBI's taxdump.tar.gz or a directory with nodes.dmp and names.dmp.
    """
    taxdump = Path(taxdump)
    parent, rank, mgc = {}, {}, {}
    with _open_taxdump(taxdump, "nodes.dmp") as handle:
        for f in _dmp(handle):
            taxid = int(f[0])
            parent[taxid], rank[taxid], mgc[taxid] = int(f[1]), f[2], int(f[8])

    children = defaultdict(list)
    for taxid, up in parent.items():
        if taxid != up:
            children[up].append(taxid)
    order, stack = [], [root]
    while stack:
        node = stack.pop()
        order.append(node)
        stack.extend(children.get(node, ()))

    # codes used at or below each node, at the ranks that count
    # (children come before their parent in reversed preorder)
    below: dict[int, set[int]] = {}
    emitted: dict[int, frozenset[int]] = {}
    for node in reversed(order):
        codes = below.pop(node, set())
        if rank[node] in CODE_RANKS and mgc[node]:
            codes.add(mgc[node])
        if rank[node] in EMIT_RANKS:
            emitted[node] = frozenset(codes)
        if node != root and codes:
            below.setdefault(parent[node], set()).update(codes)

    names = defaultdict(list)
    with _open_taxdump(taxdump, "names.dmp") as handle:
        for f in _dmp(handle):
            taxid = int(f[0])
            if taxid in emitted and f[3] in NAME_CLASSES:
                names[f[1]].append(taxid)

    rows = []
    for name in sorted(names):
        taxids = sorted(set(names[name]))
        codes = sorted(set().union(*(emitted[t] for t in taxids)))
        if not codes:
            continue
        rows.append({
            "name": name,
            "rank": ";".join(sorted({rank[t] for t in taxids})),
            "code": str(codes[0]) if len(codes) == 1 else MIXED,
            "codes": ";".join(str(c) for c in codes),
            "ncbi_taxid": ";".join(str(t) for t in taxids),
        })
    return rows


def write_table(rows: list[dict[str, str]], path: str | Path, source: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        handle.write(f"# Mitochondrial genetic code per taxon name, from {source} ({date.today().isoformat()}).\n")
        handle.write("# Built with: python -m timegraft.geneticcode. 'mixed': descendants use more than one code.\n")
        writer = csv.DictWriter(handle, fieldnames=["name", "rank", "code", "codes", "ncbi_taxid"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Build the taxon-to-genetic-code mapping from the NCBI taxdump.")
    parser.add_argument("taxdump", help="taxdump.tar.gz, or a directory with nodes.dmp and names.dmp")
    parser.add_argument("output", help="output TSV")
    parser.add_argument("--root", type=int, default=2759, help="NCBI taxid to restrict to (default 2759, Eukaryota)")
    args = parser.parse_args(argv)
    rows = build_table(args.taxdump, args.root)
    write_table(rows, args.output, f"NCBI taxdump {Path(args.taxdump).name}, root {args.root}")
    print(f"{len(rows)} names, {sum(r['code'] == MIXED for r in rows)} mixed")


if __name__ == "__main__":
    main()
