"""Record filters, the 'filter and extract' marker plug-in.

COI: length, ambiguity, and a reading frame free of stop codons under the
record's genetic code, which the caller resolves per record from its taxonomy
(see ``timegraft.geneticcode``). No profile HMM is used.

ITS: length and ambiguity. Extraction of ITS1/ITS2 with ITSx and chimera
checks are expected upstream, or as a separate rule.
"""

from __future__ import annotations

from timegraft import geneticcode

NUCLEOTIDES = frozenset("ACGT")


def _degap(seq: str) -> str:
    return seq.replace("-", "").replace(".", "")


def ambiguous_fraction(seq: str) -> float:
    seq = _degap(seq)
    if not seq:
        return 1.0
    return sum(1 for c in seq if c not in NUCLEOTIDES) / len(seq)


def has_open_frame(seq: str, genetic_code: int = 5) -> bool:
    """True if any of the three forward frames has no stop codon."""
    return geneticcode.has_open_frame(seq, genetic_code)


def check_coi(seq: str, min_length: int = 200, max_ambiguous: float = 0.01, genetic_code: int = 5, **_) -> str | None:
    """Return None if the record passes, otherwise the reason for rejection."""
    if len(_degap(seq)) < min_length:
        return "too_short"
    if ambiguous_fraction(seq) > max_ambiguous:
        return "too_ambiguous"
    if not has_open_frame(seq, genetic_code):
        return "stop_codon"
    return None


def check_its(seq: str, min_length: int = 100, max_ambiguous: float = 0.01, **_) -> str | None:
    if len(_degap(seq)) < min_length:
        return "too_short"
    if ambiguous_fraction(seq) > max_ambiguous:
        return "too_ambiguous"
    return None


CHECKS = {"coi": check_coi, "its": check_its}

# Markers whose check takes a genetic code, resolved per record.
USES_GENETIC_CODE = frozenset({"coi"})


def get_check(marker: str):
    try:
        return CHECKS[marker]
    except KeyError:
        raise ValueError(f"no filter for marker {marker!r}; choose from {sorted(CHECKS)}")
