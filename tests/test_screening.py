"""Screening units for misidentifications and outliers (issue #38)."""

import shutil

import dendropy
import pytest

from timegraft.dating import is_ultrametric, root_age
from timegraft.placement import attach_to_genus_stem
from timegraft.screening import screen, search_units, set_aside


def decisions(rows):
    return {r["unit_id"]: r["decision"] for r in rows}


def test_unit_closer_to_another_genus_is_set_aside():
    genus = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B"}
    hits = {
        "A1": [("A2", 0.93), ("A3", 0.91)],
        "A2": [("A1", 0.93), ("A3", 0.92)],
        "A3": [("B1", 0.95), ("A2", 0.92)],  # 3% closer to B: a misidentified B?
        "B1": [("A3", 0.95), ("B2", 0.94)],  # 1% closer to A: within the margin, kept
        "B2": [("B1", 0.94)],
    }
    rows = screen(genus, hits, margin=0.02, min_identity=0.80)
    assert decisions(rows) == {"A1": "keep", "A2": "keep", "A3": "other_genus", "B1": "keep", "B2": "keep"}
    row = next(r for r in rows if r["unit_id"] == "A3")
    assert (row["best_congener"], row["best_congener_identity"], row["best_other_genus"], row["best_other_identity"]) == \
        ("A2", "0.9200", "B", "0.9500")
    assert set_aside(rows) == {"A3"}


def test_isolated_unit_is_set_aside_only_when_the_genus_has_a_core():
    genus = {"A1": "A", "A2": "A", "A3": "A", "C1": "C", "C2": "C", "S1": "S"}
    hits = {
        "A1": [("A2", 0.95)], "A2": [("A1", 0.95)],
        "A3": [],  # nothing at 80%: isolated from the A1-A2 core
        "C1": [], "C2": [],  # two units far apart: no core, nothing to be isolated from
        "S1": [("A1", 0.85)],  # single-unit genus: not screened
    }
    assert decisions(screen(genus, hits)) == {
        "A1": "keep", "A2": "keep", "A3": "isolated", "C1": "keep", "C2": "keep", "S1": "keep"}


def test_a_genus_never_loses_all_its_units():
    genus = {"R1": "Radix", "R2": "Radix", "G1": "Gyraulus", "G2": "Gyraulus", "G3": "Gyraulus"}
    hits = {
        "R1": [("G1", 0.99), ("R2", 0.90)], "R2": [("G2", 0.98), ("R1", 0.90)],
        "G1": [("R1", 0.99), ("G2", 0.95)], "G2": [("G1", 0.95), ("R2", 0.98)], "G3": [("G1", 0.94)],
    }
    rows = decisions(screen(genus, hits))
    assert rows["R1"] == rows["R2"] == "kept_whole_genus"
    assert rows["G1"] == "other_genus" and rows["G3"] == "keep"
    # self hits and hits to unknown units are ignored, and a hit found one way counts both ways
    rows = decisions(screen({"A1": "A", "A2": "A", "A3": "A"}, {"A1": [("A1", 1.0), ("X9", 0.99), ("A2", 0.9)]}))
    assert rows == {"A1": "keep", "A2": "keep", "A3": "isolated"}


@pytest.mark.skipif(shutil.which("vsearch") is None, reason="needs vsearch")
def test_search_units_drops_self_hits():
    import random

    rng = random.Random(5)
    base = [rng.choice("ACGT") for _ in range(400)]

    def mutate(seq, n):
        seq = list(seq)
        for i in rng.sample(range(len(seq)), n):
            seq[i] = rng.choice([b for b in "ACGT" if b != seq[i]])
        return "".join(seq)

    reps = [("u1", mutate(base, 4)), ("u2", mutate(base, 8)), ("u3", mutate(base, 200))]
    hits = search_units(reps, min_identity=0.80)
    assert all(t != q for q, h in hits.items() for t, _ in h)
    assert [t for t, _ in hits["u1"]] == ["u2"] and not hits.get("u3")


def test_isolated_unit_goes_on_the_stem_of_its_genus():
    # genus G: crown at 4, stem node at 10; genus H: one unit, stem node at 10
    tree = dendropy.Tree.get(data="(((G1:1,G2:1):3,G3:4):6,H1:10);", schema="newick", rooting="force-rooted")
    age = attach_to_genus_stem(tree, "G4", ["G1", "G2", "G3", "G4"], fraction=0.5)
    assert age == pytest.approx(7.0)  # halfway between crown (4) and stem node (10)
    assert is_ultrametric(tree) and root_age(tree) == pytest.approx(10.0)
    assert {l.taxon.label for l in tree.mrca(taxon_labels=["G1", "G4"]).leaf_iter()} == {"G1", "G2", "G3", "G4"}
    assert attach_to_genus_stem(tree, "H2", ["H1", "H2"]) == pytest.approx(5.0)  # single unit: its tip is the crown
    assert attach_to_genus_stem(tree, "Z2", ["Z1", "Z2"]) is None
