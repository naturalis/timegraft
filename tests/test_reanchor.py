"""Re-anchoring flagged genera before dating (#50)."""

from pathlib import Path

import dendropy
import pytest

from timegraft import concordance as cc
from timegraft import reanchor as ra
from timegraft.otol_static import SynthTree, Taxonomy, genus_anchors, genus_topology_newick

DATA = Path(__file__).parent / "data" / "otol"
GENERA = ["Mytilus", "Ostrea", "Magallana", "Asterias"]


@pytest.fixture(scope="module")
def taxonomy():
    return Taxonomy(DATA / "taxonomy.tsv", set(GENERA), context="Metazoa")


@pytest.fixture(scope="module")
def synth():
    return SynthTree(DATA / "labelled_supertree.tre")


def test_read_overrides(tmp_path):
    path = tmp_path / "overrides.tsv"
    path.write_text("genus\taction\ttarget\tsource\tnote\n"
                    "Onchidella\treanchor\tfamily:Onchidiidae\tRV 2026-11-12\twrong NCBI merge in OTT\n"
                    "Magelona\tkeep\t\tRV\tCOI artefact, placement fine\n")
    table = ra.read_overrides(path)
    assert table["Onchidella"]["target"] == "family:Onchidiidae" and table["Magelona"]["action"] == "keep"
    assert ra.read_overrides(None) == {}
    path.write_text("genus\taction\ttarget\nX\tmove\t\n")
    with pytest.raises(ValueError):
        ra.read_overrides(path)
    path.write_text("genus\taction\ttarget\nX\treanchor\t\n")
    with pytest.raises(ValueError):
        ra.read_overrides(path)


def test_find_target(taxonomy):
    names = cc.name_index(taxonomy)
    assert ra.find_target(taxonomy, names, {"family": "Ostreidae", "class": "Bivalvia"}) == (
        8, "family:Ostreidae", "reference_lineage")
    # a family OTT does not hold: the class instead, if its lineage agrees with the reference phylum
    assert ra.find_target(taxonomy, names, {"family": "Nowhereidae", "class": "Asteroidea",
                                            "phylum": "Metazoa"})[:2] == (13, "class:Asteroidea")
    assert ra.find_target(taxonomy, names, {}) == (None, "", "no_target")
    # a unique name whose lineage disagrees with the reference is not taken (it may be a homonym elsewhere)
    # (the last name tried is reported)
    assert ra.find_target(taxonomy, names, {"family": "Ostreidae", "class": "Asteroidea"}) == (
        None, "class:Asteroidea", "no_agreeing_target")
    assert ra.find_target(taxonomy, names, {"class": "Asteroidea"}) == (None, "class:Asteroidea",
                                                                        "no_agreeing_target")
    assert ra.find_target(taxonomy, names, {}, "class:Asteroidea") == (13, "class:Asteroidea", "override")
    assert ra.find_target(taxonomy, names, {}, "ott:8") == (8, "family:Ostreidae", "override")
    assert ra.find_target(taxonomy, names, {}, "family:Nope")[2] == "override_not_found"
    # override ids that the loaded taxonomy does not hold are refused, not looked up
    for bad in ("ott:-1", "ott:999999", "ott:abc", f"ott:{len(taxonomy.parent) - 1 + 1000}"):
        assert ra.find_target(taxonomy, names, {}, bad) == (None, bad, "override_not_found")


def test_find_target_passes_over_an_unresolvable_rank(taxonomy):
    # pretend OTT also holds an Ostreidae among the starfish (a homonym)
    names = dict(cc.name_index(taxonomy))
    names[("family", "Ostreidae")] = [8, 13]
    # the reference class tells the two apart
    assert ra.find_target(taxonomy, names, {"family": "Ostreidae", "class": "Bivalvia"}) == (
        8, "family:Ostreidae", "reference_lineage")
    # only the shared kingdom to go on: neither Ostreidae agrees uniquely, so the class is used instead
    assert ra.find_target(taxonomy, names, {"family": "Ostreidae", "class": "Asteroidea", "phylum": "Metazoa"}) == (
        13, "class:Asteroidea", "reference_lineage")


@pytest.mark.parametrize("value, expected", [(False, "off"), (True, "apply"), ("report", "report"),
                                             ("off", "off")])
def test_yaml_off_is_off(value, expected):
    assert ra.normalise_mode(value) == expected


def test_synthesis_node(taxonomy, synth):
    node, how = ra.synthesis_node(synth, taxonomy, 13)
    assert synth.label[node] == "ott13" and how == "taxon_in_synthesis"
    # Ostreidae is not a node: the MRCA of its members (Ostrea and Magallana's species)
    node, how = ra.synthesis_node(synth, taxonomy, 8)
    assert synth.label[node] == "mrcaott6ott12" and how == "mrca_of_members"
    # Mytilidae: its only other genus is not in the synthesis and Mytilus itself is left out
    node, how = ra.synthesis_node(synth, taxonomy, 4, exclude=5)
    assert synth.label[node] == "ott2" and how == "ancestor_in_synthesis"


def test_decide():
    concordance = [
        {"genus": "A", "concordance_flag": "flag", "severity": "5", "n_flags": "3", "n_evaluable": "3"},
        {"genus": "B", "concordance_flag": "flag", "severity": "3", "n_flags": "2", "n_evaluable": "3"},
        {"genus": "C", "concordance_flag": "review", "severity": "5", "n_flags": "1", "n_evaluable": "3"},
        {"genus": "D", "concordance_flag": "flag", "severity": "5", "n_flags": "2", "n_evaluable": "2"},
        {"genus": "E", "concordance_flag": "", "severity": "0", "n_flags": "0", "n_evaluable": "3"},
    ]
    overrides = {"D": {"action": "keep", "target": ""}, "E": {"action": "reanchor", "target": "family:X"}}
    out = ra.decide(concordance, overrides, cc.LEVELS["class"], {g: 1 for g in "ABCDE"})
    assert out["A"]["decision"] == "reanchor_auto" and "severity=5" in out["A"]["trigger"]
    assert "B" not in out and "C" not in out  # below the severity, or only one failed check
    assert out["D"]["decision"] == "kept_override"
    assert out["E"]["decision"] == "reanchor_override"


def test_reanchored_genus_moves_in_the_topology(taxonomy, synth):
    # Mytilus moved to the starfish (as if OTT had misplaced it), through the curated table
    anchors, status = genus_anchors(synth, taxonomy, GENERA)
    names = cc.name_index(taxonomy)
    uid, _, _ = ra.find_target(taxonomy, names, {}, "class:Asteroidea")
    node, _ = ra.synthesis_node(synth, taxonomy, uid, exclude=taxonomy.genus["Mytilus"]["ott_id"])
    anchors["Mytilus"], status["Mytilus"] = node, f"reanchored:{synth.label[node]}"
    newick, placed = genus_topology_newick(synth, taxonomy, GENERA, (anchors, status))
    assert placed["Mytilus"] == "reanchored:ott13"
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True)
    sisters = cc.backbone_sisters(tree)
    assert sisters["Mytilus"] == ["Asterias"] and sisters["Asterias"] == ["Mytilus"]
    assert set(sisters["Ostrea"]) == {"Magallana"}
    # the enclosing taxa for TimeTree are now those of the new place
    assert [u for u, _ in ra.lineage_from(taxonomy, uid)] == [13, 2, 1]


def _assess(taxonomy, synth, anchors_status, hits):
    """Concordance on the topology from the given anchors, with Mytilus's reference lineage that of a mussel."""
    genera = GENERA
    newick, placed = genus_topology_newick(synth, taxonomy, genera, anchors_status)
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True)
    matches = {g: {**taxonomy.genus[g], "status": placed[g]} for g in genera}
    reference = {"Mytilus": ("col", {"family": "Mytilidae", "class": "Bivalvia", "phylum": "Mollusca"}),
                 "Ostrea": ("col", {"family": "Ostreidae", "class": "Bivalvia", "phylum": "Mollusca"}),
                 "Magallana": ("col", {"family": "Ostreidae", "class": "Bivalvia", "phylum": "Mollusca"}),
                 "Asterias": ("col", {"family": "Asteriidae", "class": "Asteroidea", "phylum": "Echinodermata"})}
    phylum = {g: ref["phylum"] for g, (_, ref) in reference.items()}
    gastropod = [{"name": "Mollusca", "rank": "phylum"}, {"name": "Bivalvia", "rank": "class"}]
    starfish = [{"name": "Echinodermata", "rank": "phylum"}, {"name": "Asteroidea", "rank": "class"}]
    merged = {"Mytilus": ["6548", "999"]}
    names = cc.name_index(taxonomy)
    from timegraft.otol_static import lineage_names

    rows = cc.assess(taxonomy, names, genera, matches, reference, merged, {"6548": gastropod, "999": starfish},
                     hits, cc.backbone_sisters(tree), phylum, lineage_names, flag_level=cc.LEVELS["class"],
                     max_share=0.0, min_hits=5, min_flags=2)
    return rows, matches, names, reference, tree


@pytest.mark.parametrize("mode", ["apply", "report"])
def test_misplaced_genus_is_flagged_and_moved(tmp_path, mode):
    # OTT has put Mytilus among the starfish (as it put Onchidella among the sea urchins)
    taxonomy = Taxonomy(DATA / "taxonomy.tsv", set(GENERA), context="Metazoa")
    path = tmp_path / "synth.tre"
    path.write_text("(((ott10)ott9,ott12)ott3,((ott15)ott14,(ott6,ott7)ott5)ott13)ott2;")
    synth = SynthTree(path)
    anchors, status = genus_anchors(synth, taxonomy, GENERA)
    hits = {"Mytilus": ["Ostrea"] * 4 + ["Magallana"] * 2}
    rows, matches, names, reference, tree = _assess(taxonomy, synth, (anchors, status), hits)
    by = {r["genus"]: r for r in rows}
    # OTT's taxonomy is right here, only the synthesis is wrong: merged ids and barcode neighbours flag it
    assert by["Mytilus"]["concordance_flag"] == "flag" and by["Mytilus"]["n_flags"] == "2"
    assert by["Mytilus"]["lineage_flag"] == "" and by["Mytilus"]["neighbour_flag"] == "yes"
    assert [r["genus"] for r in rows if r["concordance_flag"] == "flag"] == ["Mytilus"]
    provenance, moved = ra.plan(synth, taxonomy, names, matches, anchors, rows,
                                {g: ref for g, (_, ref) in reference.items()}, {}, mode, cc.LEVELS["class"])
    row = provenance["Mytilus"]
    assert row["original_anchor"] == "ott5" and row["reanchor_source"] == "auto"
    assert row["reanchor_target"] == "family:Mytilidae (ott4)"
    assert row["reanchor_how"] == "reference_lineage;ancestor_in_synthesis"
    assert "severity=5" in row["reanchor_trigger"]
    assert set(provenance) == {"Mytilus"}  # nothing else is touched
    if mode == "report":
        assert row["placement_action"] == "would_reanchor" and moved == {}
        return
    assert row["placement_action"] == "reanchored" and synth.label[moved["Mytilus"][0]] == "ott3"
    anchors["Mytilus"], status["Mytilus"] = moved["Mytilus"][0], "reanchored:ott3"
    after, *_ = _assess(taxonomy, synth, (anchors, status), hits)
    after = {r["genus"]: r for r in after}
    # no longer flagged: the neighbours agree; the merged ids still differ, as they are a property of the OTT record
    assert after["Mytilus"]["concordance_flag"] == "review" and after["Mytilus"]["neighbour_concordance"] == "1.000"
    assert after["Asterias"]["concordance_flag"] == ""


def test_override_keeps_a_flagged_genus(tmp_path):
    taxonomy = Taxonomy(DATA / "taxonomy.tsv", set(GENERA), context="Metazoa")
    path = tmp_path / "synth.tre"
    path.write_text("(((ott10)ott9,ott12)ott3,((ott15)ott14,(ott6,ott7)ott5)ott13)ott2;")
    synth = SynthTree(path)
    anchors, status = genus_anchors(synth, taxonomy, GENERA)
    rows, matches, names, reference, _ = _assess(taxonomy, synth, (anchors, status),
                                                 {"Mytilus": ["Ostrea"] * 6})
    keep = {"Mytilus": {"action": "keep", "target": "", "source": "RV", "note": "checked"}}
    provenance, moved = ra.plan(synth, taxonomy, names, matches, anchors, rows,
                                {g: ref for g, (_, ref) in reference.items()}, keep, "apply", cc.LEVELS["class"])
    assert moved == {} and provenance["Mytilus"]["placement_action"] == "kept_override"
    assert provenance["Mytilus"]["reanchor_source"] == "override:RV"
    assert ra.plan(synth, taxonomy, names, matches, anchors, rows, {}, {}, "off", 4) == ({}, {})
    assert ra.plan(synth, taxonomy, names, matches, anchors, rows, {}, {}, False, 4) == ({}, {})
