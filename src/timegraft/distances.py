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


_CODE = np.full(256, 4, dtype=np.uint8)  # 0..3 = A, C, G, T; 4 = gap or ambiguity
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i


def k2p_matrix(seqs: Sequence[str], max_distance: float = DEFAULT_MAX_DISTANCE) -> np.ndarray:
    """All pairwise K2P distances at once, vectorised; same result as ``k2p`` per pair."""
    if len({len(s) for s in seqs}) > 1:
        raise ValueError("k2p needs aligned sequences of equal length")
    codes = np.stack([_CODE[np.frombuffer(s.encode(), dtype=np.uint8)] for s in seqs])
    valid = codes < 4
    purine = (codes == 0) | (codes == 2)
    n = len(seqs)
    out = np.zeros((n, n))
    for i in range(n - 1):
        a, rest = codes[i], codes[i + 1 :]
        both = valid[i] & valid[i + 1 :]
        sites = both.sum(axis=1)
        diff = both & (rest != a)
        transitions = (diff & (purine[i + 1 :] == purine[i])).sum(axis=1)
        transversions = diff.sum(axis=1) - transitions
        if (sites == 0).any():
            raise ValueError("no comparable sites between sequences")
        p = transitions / sites
        q = transversions / sites
        w1, w2 = 1.0 - 2.0 * p - q, 1.0 - 2.0 * q
        with np.errstate(divide="ignore", invalid="ignore"):
            d = -0.5 * np.log(w1) - 0.25 * np.log(w2)
        d = np.where((w1 <= 0) | (w2 <= 0), max_distance, np.minimum(d, max_distance))
        out[i, i + 1 :] = d
        out[i + 1 :, i] = d
    return out


def _kmer_codes(seq: str, k: int) -> np.ndarray:
    """The distinct k-mers of ``seq`` without gaps, as integers (base 4), skipping those with an ambiguity code;
    the same k-mers as :func:`_kmers`."""
    codes = _CODE[np.frombuffer(seq.replace("-", "").encode(), dtype=np.uint8)].astype(np.int64)
    if len(codes) < k:
        return np.zeros(0, dtype=np.int64)
    windows = np.lib.stride_tricks.sliding_window_view(codes, k)
    windows = windows[(windows < 4).all(axis=1)]
    return np.unique(windows @ (4 ** np.arange(k - 1, -1, -1, dtype=np.int64)))


# above this many cells (sequences times distinct k-mers) the records-by-k-mers matrix of kmer_matrix() is not
# built, and intersections come from the k-mer sets instead: the same counts in less memory
DENSE_CELLS = 50_000_000  # 200 MB of float32


def kmer_matrix(seqs: Sequence[str], k: int = 8, max_distance: float = DEFAULT_MAX_DISTANCE,
                dense_cells: int = DENSE_CELLS) -> np.ndarray:
    """All pairwise k-mer distances at once; the same result as ``kmer`` per pair, bit for bit (#112).

    Each sequence's k-mer set is built once. The intersections are exact counts: from one product of a
    records-by-k-mers matrix, or from the sets where that matrix would exceed ``dense_cells`` cells (large
    genera with many distinct k-mers). The Mash formula is applied as in ``kmer``.
    """
    n = len(seqs)
    if n < 2:
        return np.zeros((n, n))
    sets = [_kmer_codes(s, k) for s in seqs]
    if any(len(x) == 0 for x in sets):
        raise ValueError("sequence shorter than k")
    columns, index = np.unique(np.concatenate(sets), return_inverse=True)
    if n * len(columns) <= dense_cells:
        present = np.zeros((n, len(columns)), dtype=np.float32)  # exact for counts below 2**24
        rows = np.repeat(np.arange(n), [len(x) for x in sets])
        present[rows, index] = 1.0
        inter = (present @ present.T).astype(np.int64)
    else:
        as_sets = [set(x.tolist()) for x in sets]
        inter = np.zeros((n, n), dtype=np.int64)
        for i in range(n):
            inter[i, i] = len(as_sets[i])
            for j in range(i + 1, n):
                inter[i, j] = inter[j, i] = len(as_sets[i] & as_sets[j])
    size = np.array([len(x) for x in sets], dtype=np.int64)
    union = size[:, None] + size[None, :] - inter
    jaccard = inter / union  # true division of exact counts: correctly rounded, as in kmer()
    ratio = 2.0 * jaccard / (1.0 + jaccard)
    # math.log, not numpy.log, which can differ in the last bit and so change which record is the medoid; taken
    # once per distinct value
    values, where = np.unique(ratio, return_inverse=True)
    logs = np.array([math.log(v) if v > 0 else -math.inf for v in values])[where].reshape(ratio.shape)
    d = -1.0 / k * logs
    out = np.where(jaccard == 0, max_distance, np.minimum(d, max_distance))
    np.fill_diagonal(out, 0.0)
    return out


def matrix(seqs: Sequence[str], method: str = "k2p", **options) -> np.ndarray:
    """Symmetric distance matrix for a list of sequences."""
    if method == "k2p":
        return k2p_matrix(seqs, **{k: v for k, v in options.items() if k == "max_distance"})
    if method == "kmer":
        return kmer_matrix(seqs, **{k: v for k, v in options.items() if k in ("k", "max_distance")})
    func = get_method(method)
    n = len(seqs)
    out = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            out[i, j] = out[j, i] = func(seqs[i], seqs[j], **options)
    return out
