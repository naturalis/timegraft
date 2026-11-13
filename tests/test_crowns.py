"""Crown ages from the species sampled, and a plausibility test against the barcodes (#54)."""

import math
from pathlib import Path

import dendropy
import pytest

from timegraft.build import crown_decision
from timegraft.dating import crown_age, saturation_factor, timetree_plausibility
from timegraft.otol_static import CROWN_COLUMNS, Taxonomy, crown_ages
from timegraft.timetree import TimeTree

DATA = Path(__file__).parent / "data" / "otol"
ON = {"enabled": True, "rate_factor": 10, "saturation": 0.15}


def test_saturation_factor():
    assert saturation_factor(0.0, 0.15) == 1.0
    assert saturation_factor(0.0133, 0.15) == pytest.approx(1.047, abs=1e-3)  # saturation hides almost nothing
    assert saturation_factor(0.137, 0.15) == pytest.approx(2.678, abs=1e-3)
    assert saturation_factor(0.15, 0.15) == math.inf
    # grows with distance
    assert saturation_factor(0.05, 0.15) < saturation_factor(0.1, 0.15) < saturation_factor(0.14, 0.15)


@pytest.mark.parametrize("tt, barcode, d, verdict", [
    (69.05, 1.16, 0.0133, "implausible"),   # Ruditapes in the analysis: 59 times the barcode age at 1.3%
    (443.7, 11.93, 0.137, "implausible"),   # Haliclona: 37 times, saturation allows about 27
    (193.1, 8.12, 0.0934, "implausible"),   # Hydra
    (20.0, 5.0, 0.05, "plausible"),         # a typical genus: four times the barcode minimum
    (60.0, 5.0, 0.14, "plausible"),         # near saturation, much older crowns are possible
    (3.0, 5.0, 0.05, "plausible"),          # younger than the barcode minimum, within a faster rate
    (1.15, 13.8, 0.16, "implausible_young"),  # Gammarus: 4 of 184 species known to TimeTree, all close
])
def test_timetree_plausibility(tt, barcode, d, verdict):
    got, ratio, allowed = timetree_plausibility(tt, barcode, d)
    assert got == verdict
    assert ratio == pytest.approx(tt / barcode)


def _tree(d):
    return dendropy.Tree.get(data=f"(A:{d},B:{d});", schema="newick", rooting="force-rooted")


def test_implausible_timetree_crown_is_demoted_to_an_upper_bound():
    # barcode age 1.3 Ma (d 0.013 at rate 0.01); TimeTree says 69 Ma, stem 232 Ma
    kw = dict(method="rate", rate=0.01, suggested={"G": 69.0}, estimate="timetree_midpoint")
    used = crown_age("G", 232, _tree(0.013), {}, **kw)
    assert (used.source, used.age, used.plausibility) == ("timetree", 69.0, "")  # no test, as before
    demoted = crown_age("G", 232, _tree(0.013), {}, plausibility=ON, **kw)
    assert demoted.source == "midpoint_timetree_demoted" and demoted.plausibility == "implausible"
    assert demoted.age == pytest.approx(math.sqrt(1.3 * 69.0))
    assert demoted.upper == 69.0 and demoted.upper_source == "timetree"  # still the upper bound
    assert demoted.timetree_ratio == pytest.approx(69.0 / 1.3) and demoted.ratio_allowed < 11
    # a plausible TimeTree crown stays the answer
    kept = crown_age("G", 232, _tree(0.05), {}, plausibility=ON, **{**kw, "suggested": {"G": 20.0}})
    assert (kept.source, kept.age, kept.plausibility) == ("timetree", 20.0, "plausible")
    off = crown_age("G", 232, _tree(0.013), {}, plausibility={"enabled": False}, **kw)
    assert off.source == "timetree" and off.plausibility == ""
    # the fraction method has no age from the barcodes to test against
    frac = crown_age("G", 232, _tree(0.013), {}, method="fraction", fraction=0.5, suggested={"G": 69.0},
                     estimate="timetree_midpoint", plausibility=ON)
    assert frac.source == "timetree" and frac.plausibility == ""
    # with other estimates TimeTree is only ever an upper bound; the test is still recorded
    mid = crown_age("G", 232, _tree(0.013), {}, plausibility=ON, **{**kw, "estimate": "midpoint"})
    assert mid.source == "midpoint" and mid.plausibility == "implausible"


def test_too_young_timetree_crown_is_not_used():
    # barcode age 13.8 Ma; TimeTree 1.15 Ma, from species that do not span the crown; stem 53 Ma
    kw = dict(method="rate", rate=0.01, suggested={"G": 1.15}, estimate="timetree_midpoint", plausibility=ON)
    young = crown_age("G", 53, _tree(0.138), {}, **kw)
    assert young.plausibility == "implausible_young" and young.upper_source == "stem"
    assert young.source == "midpoint" and young.age == pytest.approx(math.sqrt(13.8 * 0.9 * 53))
    assert "too young" in crown_decision(young, "sampled") and "stem cap" in crown_decision(young, "sampled")
    before = crown_age("G", 53, _tree(0.138), {}, **{**kw, "plausibility": None})
    assert (before.source, before.age) == ("timetree", 1.15)


def test_crown_decision():
    kw = dict(method="rate", rate=0.01, suggested={"G": 69.0}, estimate="timetree_midpoint", plausibility=ON)
    demoted = crown_age("G", 232, _tree(0.013), {}, **kw)
    text = crown_decision(demoted, "sampled")
    assert text.startswith("timetree[sampled] implausible (ratio 53.1 > ") and "midpoint of barcodes" in text
    kept = crown_age("G", 232, _tree(0.05), {}, **{**kw, "suggested": {"G": 20.0}})
    assert crown_decision(kept, "global_fallback") == "timetree[global_fallback] (plausible)"
    none = crown_age("G", 232, _tree(0.05), {}, method="rate", rate=0.01, estimate="timetree_midpoint")
    assert crown_decision(none) == "midpoint of barcodes and stem cap"
    assert crown_decision(crown_age("G", 50, _tree(0.05), {"G": 20.0}, method="rate", rate=0.01)) == "table"


def _fake(ages):
    def fetch(url, payload=None):
        ids = [int(i) for i in url.split("/mrca/id/")[1].split("/")[0].split("+")]
        return {"found_ids": ids, "precomputed_age": ages.get(frozenset(ids))}

    return fetch


@pytest.fixture(scope="module")
def taxonomy():
    return Taxonomy(DATA / "taxonomy.tsv", {"Mytilus", "Ostrea"}, context="Metazoa")


def test_sampled_species_crown(taxonomy):
    # TimeTree knows the two Mytilus species of OTT (6550, 6551), and a third (6549, from COL) that is
    # the one sampled with M. edulis
    tt = TimeTree(fetch=_fake({frozenset({6550, 6551}): 5.0, frozenset({6549, 6550}): 1.5}))
    rows = crown_ages(taxonomy, tt, ["Mytilus"], sampled={"Mytilus": {"Mytilus edulis", "Mytilus galloprovincialis"}},
                      species_ncbi={"Mytilus galloprovincialis": 6549})
    assert len(rows) == 1 and set(rows[0]) == set(CROWN_COLUMNS)
    row = rows[0]
    assert (row["crown_age"], row["crown_basis"], row["n_sampled_species"]) == ("1.5000", "sampled", 2)
    assert row["sampled_crown_age"] == "1.5000" and row["global_crown_age"] == "5.0000"
    assert row["query_ids"] == "6549+6550"
    # the global basis, as before #54
    row = crown_ages(taxonomy, tt, ["Mytilus"], sampled={"Mytilus": {"Mytilus edulis"}}, basis="global")[0]
    assert (row["crown_age"], row["crown_basis"], row["sampled_crown_age"]) == ("5.0000", "global", "")


def test_sampled_species_fallbacks(taxonomy):
    tt = TimeTree(fetch=_fake({frozenset({6550, 6551}): 5.0}))
    # one species sampled: its subtree is within a species, so no TimeTree crown
    row = crown_ages(taxonomy, tt, ["Mytilus"], sampled={"Mytilus": {"Mytilus edulis"}})[0]
    assert (row["crown_age"], row["crown_basis"], row["n_sampled_species"]) == ("", "single_sampled_species", 1)
    row = crown_ages(taxonomy, tt, ["Mytilus"], sampled={"Mytilus": {"Mytilus edulis"}}, single_species="global")[0]
    assert (row["crown_age"], row["crown_basis"]) == ("5.0000", "global_fallback")
    # sampled species TimeTree cannot date: the global crown
    row = crown_ages(taxonomy, tt, ["Mytilus"], sampled={"Mytilus": {"Mytilus edulis", "Mytilus x"}},
                     species_ncbi={"Mytilus x": 999})[0]
    assert (row["crown_age"], row["crown_basis"]) == ("5.0000", "global_fallback")
    # no species names in the records: the global crown
    row = crown_ages(taxonomy, tt, ["Mytilus"], sampled={})[0]
    assert (row["crown_age"], row["crown_basis"], row["n_sampled_species"]) == ("5.0000", "global", "")
    with pytest.raises(ValueError):
        crown_ages(taxonomy, tt, ["Mytilus"], basis="regional")
