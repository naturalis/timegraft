"""One dated genus subtree, with the spread of crown age and PD over roots (#34)."""

import numpy as np
import pytest

from timegraft.build import build_genus, tree_length
from timegraft.dating import is_ultrametric

BASES = np.array(list("ACGT"))


def mutate(rng, seq, n):
    seq = seq.copy()
    sites = rng.choice(len(seq), size=n, replace=False)
    seq[sites] = [rng.choice([b for b in "ACGT" if b != seq[i]]) for i in sites]
    return seq


@pytest.fixture
def genus():
    rng = np.random.default_rng(1)
    root = BASES[rng.integers(0, 4, 600)]
    # two clades in the genus, an outgroup further away, and one genus unit among the outgroup
    out_anc = mutate(rng, root, 90)
    clade1, clade2 = mutate(rng, root, 20), mutate(rng, root, 20)
    ingroup = [(f"G{i}", "".join(mutate(rng, c, 6))) for i, c in enumerate([clade1] * 3 + [clade2] * 3)]
    out1, out2 = mutate(rng, out_anc, 15), mutate(rng, out_anc, 15)
    outgroup = [(f"O{i}", "".join(mutate(rng, c, 6))) for i, c in enumerate([out1, out1, out2, out2])]
    ingroup.append(("G_odd", "".join(mutate(rng, out1, 6))))  # nested in the outgroup, next to O0 and O1
    return ingroup, outgroup


MARKER = {
    "align": "none",
    "distance": {"method": "k2p"},
    "crown_age": {"method": "rate", "rate": 0.0115, "cap_fraction": 0.9},
    "dating": {"method": "mpl"},
}


def test_provenance_records_root_uncertainty(genus):
    ingroup, outgroup = genus
    tree, prov = build_genus("Genus", ingroup, outgroup, 40.0, MARKER, {})
    assert is_ultrametric(tree)
    assert prov["rooting"] == "outgroup_tolerant"
    assert prov["misplaced"] == "G_odd"
    assert prov["root_alternatives"] >= 1
    assert 0.0 <= float(prov["root_ambiguity"]) <= 1.0
    lo, crown, hi = float(prov["crown_age_min"]), float(prov["crown_age"]), float(prov["crown_age_max"])
    assert lo <= crown <= hi
    # leaving the odd unit out gives a much younger genus
    assert lo < 0.8 * crown
    assert float(prov["pd_rel_range"]) >= 0.0


def test_roots_that_cannot_be_dated_are_skipped(genus, monkeypatch):
    import timegraft.build as build
    from timegraft.dating import CrownAge

    ingroup, outgroup = genus
    real = build.crown_age
    calls = []

    def crown_age(genus, stem, tree, *args, **kwargs):
        calls.append(tree)
        age = real(genus, stem, tree, *args, **kwargs)
        # every root but the chosen one comes out negative, as with negative BioNJ branch lengths
        return age if len(calls) == 1 else CrownAge(-1.0, age.source)

    monkeypatch.setattr(build, "crown_age", crown_age)
    _, prov = build.build_genus("Genus", ingroup, outgroup, 40.0, MARKER, {})
    assert len(calls) > 1
    assert prov["crown_age_min"] == prov["crown_age_max"] == prov["crown_age"]
    assert prov["pd_rel_range"] == "0.0000"


def test_crown_range_only_for_ages_from_the_barcodes(genus):
    ingroup, outgroup = genus
    _, prov = build_genus("Genus", ingroup, outgroup, 40.0, MARKER, {"Genus": 5.0})
    assert prov["crown_source"] == "table"
    assert prov["crown_age_min"] == prov["crown_age_max"] == ""
    assert prov["pd_rel_range"] != ""  # the topology still depends on the root


def test_uncertainty_can_be_switched_off(genus):
    ingroup, outgroup = genus
    tree_on, on = build_genus("Genus", ingroup, outgroup, 40.0, MARKER, {})
    tree_off, off = build_genus("Genus", ingroup, outgroup, 40.0, MARKER, {}, rooting={"uncertainty": False})
    for key in ("root_ambiguity", "root_alternatives", "crown_age_min", "crown_age_max", "pd_rel_range"):
        assert off[key] == ""
    # the subtree itself does not change
    assert on["crown_age"] == off["crown_age"]
    assert tree_length(tree_on) == pytest.approx(tree_length(tree_off))


@pytest.fixture
def long_branch():
    """Eight close units and one with many more substitutions of its own (#55)."""
    rng = np.random.default_rng(3)
    root = BASES[rng.integers(0, 4, 600)]
    out = [(f"O{i}", "".join(mutate(rng, mutate(rng, root, 80), 5))) for i in range(3)]
    ingroup = [(f"G{i}", "".join(mutate(rng, root, 8))) for i in range(8)]
    ingroup.append(("G_long", "".join(mutate(rng, root, 70))))
    return ingroup, out


def test_rate_outliers_are_listed_or_set_aside(long_branch):
    ingroup, outgroup = long_branch
    marker = {**MARKER, "screening": {"rate_outlier_z": 3.5}}
    _, flagged = build_genus("Genus", ingroup, outgroup, 40.0, marker, {})
    assert flagged["rate_outliers"] == "G_long" and flagged["rate_outlier_action"] == "flag"
    assert flagged["n_units"] == 9

    marker["screening"]["rate_outlier_action"] = "set_aside"
    tree, aside = build_genus("Genus", ingroup, outgroup, 40.0, marker, {})
    assert aside["rate_outliers"] == "G_long" and aside["rate_outlier_action"] == "set_aside"
    assert aside["n_units"] == 8 and "G_long" not in {l.taxon.label for l in tree.leaf_node_iter()}
    assert float(aside["crown_age"]) < float(flagged["crown_age"])  # the long branch stretched the genus

    with pytest.raises(ValueError):
        build_genus("Genus", ingroup, outgroup, 40.0, {**MARKER, "screening": {"rate_outlier_action": "drop"}}, {})
