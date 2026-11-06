"""Validation of the informative window."""

import dendropy
import numpy as np
import pytest

from timegraft.validation import binned_medians, compare, plateau_onset, sample_pairs, spearman, stratum


def test_stratum_is_lowest_shared_rank():
    ranks = ["genus", "family", "order"]
    a = {"genus": "G1", "family": "F1", "order": "O1"}
    assert stratum(a, {"genus": "G1", "family": "F1", "order": "O1"}, ranks) == "genus"
    assert stratum(a, {"genus": "G2", "family": "F1", "order": "O1"}, ranks) == "family"
    assert stratum(a, {"genus": "G3", "family": "F9", "order": "O2"}, ranks) == "above_order"


def test_sample_pairs_caps_each_stratum():
    units = {f"u{i}": {"genus": "G" if i < 5 else f"H{i}", "family": "F"} for i in range(10)}
    pairs = sample_pairs(units, ["genus", "family"], per_stratum=4)
    by = {}
    for _, _, s in pairs:
        by[s] = by.get(s, 0) + 1
    assert by == {"genus": 4, "family": 4}
    assert pairs == sample_pairs(units, ["genus", "family"], per_stratum=4)  # reproducible


def test_spearman_and_plateau():
    x = np.linspace(0, 10, 50)
    assert spearman(x, x ** 2) == pytest.approx(1.0)
    saturating = 1 - np.exp(-x)
    bins = binned_medians(x, saturating, n_bins=10)
    onset = plateau_onset(bins)
    assert 1.5 < onset < 3.5  # 1 - exp(-x) reaches 90% of its maximum near x = 2.3


def test_compare_on_a_small_tree():
    tree = dendropy.Tree.get(data="((a:1,b:1):4,(c:2,d:2):3);", schema="newick", rooting="force-rooted")
    seq = "ACGTTGCAAGGCTTACGATCGATCGTAGCTAGCTAGGCATCGATCGGATCGAT"
    units = {
        "a": {"sequence": seq, "genus": "G1", "family": "F1"},
        "b": {"sequence": seq[:-2] + "AA", "genus": "G1", "family": "F1"},
        "c": {"sequence": seq[::-1], "genus": "G2", "family": "F1"},
        "d": {"sequence": seq[::-1][:-3] + "TTT", "genus": "G2", "family": "F1"},
    }
    rows, summary, bins = compare(tree, units, {"kmer": {"method": "kmer", "k": 5}}, ["genus", "family"])
    assert len(rows) == 6
    assert {r["stratum"] for r in rows} == {"genus", "family"}
    assert next(r for r in rows if (r["unit_a"], r["unit_b"]) == ("a", "b"))["reference"] == pytest.approx(2.0)
    assert {s["stratum"] for s in summary} == {"all", "genus", "family"}
