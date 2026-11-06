"""Pairwise distance functions, the 'distances' marker plug-in.

Every method takes two sequences and returns a non-negative float. Methods are
looked up by name from the configuration, so a new method (for instance one
backed by precomputed embeddings) only needs to be registered here.
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

import numpy as np

PURINES = frozenset("AG")
PYRIMIDINES = frozenset("CT")
NUCLEOTIDES = frozenset("ACGT")

# Distance assigned when a pair is saturated or the formula is undefined.
DEFAULT_MAX_DISTANCE = 2.0


def k2p(a: str, b: str, max_distance: float = DEFAULT_MAX_DISTANCE, **_) -> float:
    """Kimura two-parameter distance on aligned sequences, pairwise deletion.

    Sites where either sequence has a gap or an ambiguity code are skipped.
    """
    if len(a) != len(b):
        raise ValueError("k2p needs aligned sequences of equal length")
    sites = transitions = transversions = 0
    for x, y in zip(a, b):
        if x not in NUCLEOTIDES or y not in NUCLEOTIDES:
            continue
        sites += 1
        if x == y:
            continue
        if (x in PURINES and y in PURINES) or (x in PYRIMIDINES and y in PYRIMIDINES):
            transitions += 1
        else:
            transversions += 1
    if sites == 0:
        raise ValueError("no comparable sites between sequences")
    p = transitions / sites
    q = transversions / sites
    w1 = 1.0 - 2.0 * p - q
    w2 = 1.0 - 2.0 * q
    if w1 <= 0 or w2 <= 0:
        return max_distance
    return min(-0.5 * math.log(w1) - 0.25 * math.log(w2), max_distance)


def _kmers(seq: str, k: int) -> set[str]:
    return {
        seq[i : i + k]
        for i in range(len(seq) - k + 1)
        if set(seq[i : i + k]) <= NUCLEOTIDES
    }


def kmer(a: str, b: str, k: int = 8, max_distance: float = DEFAULT_MAX_DISTANCE, **_) -> float:
    """Alignment-free distance from shared k-mers (Mash formula).

    Gaps are removed first, so aligned and unaligned input give the same result.
    """
    sa = _kmers(a.replace("-", ""), k)
    sb = _kmers(b.replace("-", ""), k)
    if not sa or not sb:
        raise ValueError("sequence shorter than k")
    jaccard = len(sa & sb) / len(sa | sb)
    if jaccard == 0:
        return max_distance
    return min(-1.0 / k * math.log(2.0 * jaccard / (1.0 + jaccard)), max_distance)


METHODS: dict[str, Callable[..., float]] = {"k2p": k2p, "kmer": kmer}
NEEDS_ALIGNMENT = {"k2p": True, "kmer": False}


def get_method(name: str) -> Callable[..., float]:
    try:
        return METHODS[name]
    except KeyError:
        raise ValueError(f"unknown distance method {name!r}; choose from {sorted(METHODS)}")


def matrix(seqs: Sequence[str], method: str = "k2p", **options) -> np.ndarray:
    """Symmetric distance matrix for a list of sequences."""
    func = get_method(method)
    n = len(seqs)
    out = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            out[i, j] = out[j, i] = func(seqs[i], seqs[j], **options)
    return out
