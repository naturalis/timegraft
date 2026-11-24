"""The vectorised k-mer distance matrix (#112): the same numbers as the pairwise ``distances.kmer``."""

import numpy as np
import pytest

from timegraft import distances


def _pairwise(seqs, **kw):
    n = len(seqs)
    out = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            out[i, j] = out[j, i] = distances.kmer(seqs[i], seqs[j], **kw)
    return out


def test_kmer_matrix_equals_pairwise_kmer():
    # #112: gaps, ambiguity codes, identical sequences, and a pair with no shared k-mers
    seqs = [
        "ACGTACGTTAGCATGCAAGT",
        "ACGT-ACGTTAGC--ATGCAAGT",     # the first with gaps: identical k-mers
        "ACGTNCGTTAGCATGCRAGT",        # ambiguity codes break k-mers
        "TTTTTTTTTTGGGGGGGGGG",        # no k-mer shared with the others
        "ACGTACGTTAGCATGGAAGTCC",
    ]
    for k in (4, 8):
        got = distances.kmer_matrix(seqs, k=k)
        assert np.array_equal(got, _pairwise(seqs, k=k))  # bit for bit, not approximately
        assert np.array_equal(distances.matrix(seqs, method="kmer", k=k), got)
    assert distances.kmer_matrix(seqs, k=8)[0, 3] == distances.DEFAULT_MAX_DISTANCE
    assert distances.kmer_matrix(seqs, k=8)[0, 1] == 0.0
    capped = distances.kmer_matrix(seqs, k=8, max_distance=0.1)
    assert np.array_equal(capped, _pairwise(seqs, k=8, max_distance=0.1))


def test_kmer_matrix_random_sequences_bit_for_bit():
    rng = np.random.default_rng(112)
    base = "".join(rng.choice(list("ACGT"), 300))
    seqs = []
    for _ in range(30):
        s = list(base)
        for i in rng.choice(300, 25, replace=False):
            s[i] = rng.choice(list("ACGTN-"))
        seqs.append("".join(s))
    assert np.array_equal(distances.kmer_matrix(seqs, k=8), _pairwise(seqs, k=8))
    # without the dense matrix (large inputs), from the k-mer sets: the same numbers
    assert np.array_equal(distances.kmer_matrix(seqs, k=8, dense_cells=0), _pairwise(seqs, k=8))


def test_kmer_matrix_edge_cases():
    assert distances.kmer_matrix([], k=8).shape == (0, 0)
    assert np.array_equal(distances.kmer_matrix(["ACGTACGTAC"], k=8), np.zeros((1, 1)))
    with pytest.raises(ValueError):
        distances.kmer_matrix(["ACGTACGTAC", "ACG"], k=8)  # shorter than k, as kmer() raises
    with pytest.raises(ValueError):
        distances.kmer_matrix(["ACGTACGTAC", "NNNNNNNNNNNN"], k=8)  # no valid k-mer
