"""Screening units for misidentifications and outliers (issue #38)."""

import shutil

import dendropy
import pytest

from timegraft.dating import is_ultrametric, root_age
from timegraft.placement import attach_to_genus_stem
from timegraft.build import rate_outliers
from timegraft.geneticcode import translate
from timegraft.screening import (
    apply_translation, assess_translation, calibrate, compare_codons, in_scope, judge_codons, read_threshold_table,
    screen, search_units, set_aside,
)


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


# Issue #55: reference scope, thresholds per class, codon diagnostics, reasons


def test_every_decision_has_its_reason_and_thresholds():
    genus = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B", "S1": "S"}
    hits = {"A1": [("A2", 0.93)], "A3": [("B1", 0.95), ("A2", 0.92)], "B1": [("B2", 0.94)]}
    rows = {r["unit_id"]: r for r in screen(genus, hits)}
    assert "by 0.0300 >= margin 0.02" in rows["A3"]["reason"] and "B at 0.9500" in rows["A3"]["reason"]
    assert rows["S1"]["reason"].startswith("only unit of its genus")
    assert all(r["reason"] and r["margin"] == "0.02" and r["min_identity"] == "0.8" for r in rows.values())
    assert rows["A3"]["best_other_scope"] == "run" and rows["A3"]["reference_scope"] == "run"


def test_thresholds_per_class():
    genus = {"A1": "A", "A2": "A", "A3": "A", "C1": "C", "C2": "C", "C3": "C"}
    hits = {"A1": [("A2", 0.95)], "A3": [("A1", 0.84)], "C1": [("C2", 0.95)], "C3": [("C1", 0.84)]}
    classes = {u: ("Copepoda" if u.startswith("C") else "Insecta") for u in genus}
    rows = {r["unit_id"]: r for r in screen(genus, hits, class_of=classes,
                                            thresholds={"Copepoda": {"min_identity": 0.86}})}
    assert rows["A3"]["decision"] == "keep" and rows["A3"]["threshold_class"] == "default"
    assert rows["C3"]["decision"] == "isolated" and rows["C3"]["threshold_class"] == "Copepoda"
    assert rows["C3"]["min_identity"] == "0.86" and rows["C3"]["margin"] == "0.02"


def test_threshold_table(tmp_path):
    path = tmp_path / "t.tsv"
    path.write_text("class\tmargin\tmin_identity\nCopepoda\t0.03\t0.82\nInsecta\t\t0.85\nGastropoda\t\t\n")
    assert read_threshold_table(path) == {"Copepoda": {"margin": 0.03, "min_identity": 0.82},
                                          "Insecta": {"min_identity": 0.85}}


def test_external_reference_adds_evidence_but_is_not_decided_on():
    genus = {"A1": "A", "A2": "A", "A3": "A"}
    # A3 is closer to an external unit of genus B than to its congeners; A2's nearest congener is external
    hits = {"A1": [("A2", 0.93)], "A3": [("x:B9", 0.97), ("A1", 0.92)], "A2": [("x:A7", 0.99)]}
    rows = {r["unit_id"]: r for r in screen(genus, hits, external={"x:B9": "B", "x:A7": "A"},
                                            reference_scope="run+external:family")}
    assert set(rows) == {"A1", "A2", "A3"}
    assert rows["A3"]["decision"] == "other_genus" and rows["A3"]["best_other_scope"] == "external"
    assert rows["A2"]["best_congener"] == "x:A7" and rows["A2"]["best_congener_scope"] == "external"
    assert rows["A1"]["reference_scope"] == "run+external:family"
    # without the reference, A3 has no other genus to be closer to
    assert {r["unit_id"]: r["decision"] for r in screen(genus, hits)}["A3"] == "keep"


def test_reference_rows_in_scope():
    ref = [{"unit_id": "U1", "genus": "A", "family": "F1"}, {"unit_id": "U2", "genus": "B", "family": "F2"},
           {"unit_id": "R1", "genus": "A", "family": "F1"}]
    run = {"family": {"F1"}}
    assert [r["unit_id"] for r in in_scope(ref, run, "family", run_units=["R1"])] == ["U1"]
    assert [r["unit_id"] for r in in_scope(ref, run, "all", run_units=["R1"])] == ["U1", "U2"]


# no stop in frame 0 under code 5; ATAGTACCTAAC puts stops in frames 1 (TAG) and 2 (TAA)
CODING = "ATGTTTCTGAAAGGTACCCGTATTGAAGCTCCGTGGAGCGATATGCATTATAACATAGTACCTAAC" * 4


def mutate_codon_positions(seq, positions, every=1):
    """Change the base at the given codon positions (0, 1, 2), in every ``every``-th codon."""
    swap = {"A": "G", "G": "A", "C": "T", "T": "C"}
    seq = list(seq)
    for i in range(0, len(seq) - 2, 3 * every):
        for p in positions:
            seq[i + p] = swap[seq[i + p]]
    return "".join(seq)


def test_codon_comparison_tells_synonymous_from_amino_acid_changes():
    third = mutate_codon_positions(CODING, [2], every=2)  # transitions at third positions: mostly silent
    first = mutate_codon_positions(CODING, [0], every=2)  # first positions: mostly change the amino acid
    s3 = compare_codons(1, third, CODING, third, 5)
    s1 = compare_codons(1, first, CODING, first, 5)
    assert s3["codons_compared"] == s1["codons_compared"] == len(CODING) // 3
    assert s3["codons_changed"] == s1["codons_changed"] == len(CODING) // 6
    assert s3["aa_changes"] / s3["codons_changed"] < 0.2 < 0.8 < s1["aa_changes"] / s1["codons_changed"]
    assert judge_codons(s3)[0] == "ok" and judge_codons(s1)[0] == "numt_like"
    # the count follows the genetic code: check against translation
    aa = sum(a != b for a, b in zip(translate(first, 5), translate(CODING, 5)))
    assert s1["aa_changes"] == aa


def test_frameshift_needs_an_indel_that_stays_out_of_frame():
    q = CODING[:60] + "-" + CODING[60:]
    t = CODING[:60] + "A" + CODING[60:]
    assert compare_codons(1, q, t, CODING, 5)["frameshift_stretches"] == 1
    # two indels that cancel out within a few bases are an artefact of the alignment
    q2 = CODING[:60] + "-" + CODING[60:66] + CODING[66:]
    t2 = CODING[:60] + CODING[60:66] + "-" + CODING[66:]
    assert compare_codons(1, q2, t2, CODING, 5)["frameshift_stretches"] == 0
    # a codon-sized indel keeps the frame
    q3 = CODING[:60] + "---" + CODING[60:]
    t3 = CODING[:60] + "AAA" + CODING[60:]
    assert compare_codons(1, q3, t3, CODING, 5)["frameshift_stretches"] == 0
    verdict, why = judge_codons({"codons_changed": 1, "aa_changes": 0, "frameshift_stretches": 1})
    assert verdict == "frameshift" and "out of frame" in why


def test_assess_translation_and_set_aside_action():
    genus = {"A1": "A", "A2": "A", "A3": "A"}
    hits = {"A1": [("A2", 0.97)], "A3": [("A1", 0.95)]}
    numt = mutate_codon_positions(CODING, [0], every=2)
    seqs = {"A1": CODING, "A2": CODING, "A3": numt}
    alignments = {("A1", "A2"): (1, CODING, CODING), ("A1", "A3"): (1, CODING, numt)}  # A3's pair seen from A1
    rows = screen(genus, hits)
    assess_translation(rows, alignments, seqs, {u: 5 for u in genus})
    by = {r["unit_id"]: r for r in rows}
    assert by["A3"]["translation"] == "numt_like" and by["A3"]["translation_ref"] == "A1"
    assert by["A1"]["translation"] == by["A2"]["translation"] == "ok"  # A2 against A1: the pair seen from A1
    apply_translation(rows, "flag")
    assert by["A3"]["decision"] == "keep"
    apply_translation(rows, "set_aside")
    assert by["A3"]["decision"] == "numt_like" and set_aside(rows) >= {"A3"}
    assert by["A3"]["reason"].startswith("codons unlike those of A1")
    with pytest.raises(ValueError):
        apply_translation(rows, "drop")


def test_numt_like_never_empties_a_genus():
    rows = [{"unit_id": u, "genus": "A", "decision": "keep", "translation": "numt_like", "translation_ref": "x",
             "translation_reason": "r", "reason": ""} for u in ("A1", "A2")]
    apply_translation(rows, "set_aside")
    assert {r["decision"] for r in rows} == {"kept_whole_genus"}


def test_calibration_summarises_identities_per_class():
    rows = [{"unit_id": f"u{i}", "genus": f"G{i % 10}", "genus_units": "3", "decision": "keep",
             "best_congener_identity": f"{0.80 + i / 1000:.4f}", "best_other_identity": "0.85",
             "nonsynonymous_share": "0.1", "translation": "ok"} for i in range(100)]
    rows.append({"unit_id": "x", "genus": "G1", "genus_units": "3", "decision": "other_genus",
                 "best_congener_identity": "0.80", "best_other_identity": "0.90"})
    out = calibrate(rows, {r["unit_id"]: "Insecta" for r in rows}, min_units=50)
    (row,) = out
    assert row["class"] == "Insecta" and row["n_units"] == 101 and row["set_aside"] == 1
    assert row["suggested_min_identity"] == "0.80"  # the 1% quantile (0.801) rounds down, and is held at the floor
    assert row["other_closer"] == 51  # 50 kept units below 0.85, and the set-aside one
    assert calibrate(rows, {}, min_units=500)[0]["suggested_min_identity"] == ""


def test_rate_outliers_flags_a_long_branch():
    tree = dendropy.Tree.get(data="((a:1,b:1.1):1,(c:0.9,d:1):1.05,e:6.0);", schema="newick")
    out = rate_outliers(tree, 3.5)
    assert out["units"] == ["e"] and float(out["max_z"]) > 3.5
    assert rate_outliers(tree, 3.5, min_units=10) == {"units": [], "max_z": ""}
    even = dendropy.Tree.get(data="((a:1,b:1):1,(c:1,d:1):1);", schema="newick")
    assert rate_outliers(even)["units"] == []


@pytest.mark.skipif(shutil.which("vsearch") is None, reason="needs vsearch")
def test_search_keeps_the_alignment_of_the_best_hit_per_side():
    import random

    rng = random.Random(7)
    base = CODING * 2
    def mutate(seq, n):
        seq = list(seq)
        for i in rng.sample(range(len(seq)), n):
            seq[i] = rng.choice([b for b in "ACGT" if b != seq[i]])
        return "".join(seq)
    reps = [("a1", mutate(base, 5)), ("a2", mutate(base, 10)), ("b1", mutate(base, 30))]
    external = [("x:a9", mutate(base, 3))]
    genus = {"a1": "A", "a2": "A", "b1": "B", "x:a9": "A"}
    hits, alignments = search_units(reps, 0.8, external=external, genus_of=genus)
    assert "x:a9" in [t for t, _ in hits["a1"]]
    for (q, t), (qilo, qrow, trow) in alignments.items():
        assert q in genus and len(qrow) == len(trow) and qilo >= 1
    sides = {(q, genus[t] == genus[q]) for q, t in alignments}
    assert len(sides) == len(alignments)  # one per query and side


def test_ambiguous_frame_takes_the_one_with_fewest_amino_acid_changes():
    seq = "ATGTTTCTGAAAGGTACCCGTATTGAAGCTCCGTGGAGCGATATGCATTATAAC" * 4  # several frames without stops
    from timegraft.geneticcode import frame_stops

    assert sorted(frame_stops(seq, 9))[:2] == [0, 0]
    changed = mutate_codon_positions(seq, [2], every=2)  # silent-ish in frame 0, not in the others
    stats = compare_codons(1, changed, seq, changed, 9)
    assert stats["frame"] == 0
    frame0 = sum(a != b for a, b in zip(translate(changed, 9), translate(seq, 9)))
    assert stats["aa_changes"] == frame0


def test_mutual_pairs_further_keeps_the_unit_closer_to_its_genus():
    genus = {"R1": "Radix", "R2": "Radix", "R3": "Radix", "G1": "Gyraulus", "G2": "Gyraulus", "G3": "Gyraulus"}
    hits = {"R1": [("G1", 0.99), ("R2", 0.90)], "R2": [("R3", 0.95)],
            "G1": [("G2", 0.96)], "G2": [("G3", 0.95)]}
    both = decisions(screen(genus, hits))
    assert both["R1"] == both["G1"] == "other_genus"  # each is the other's nearest: the pair cannot tell
    rows = {r["unit_id"]: r for r in screen(genus, hits, mutual_pairs="further")}
    assert rows["R1"]["decision"] == "other_genus"  # 0.09 above its best congener
    assert rows["G1"]["decision"] == "keep" and "closer of the two" in rows["G1"]["reason"]  # 0.03 above
    with pytest.raises(ValueError):
        screen(genus, hits, mutual_pairs="one")
