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


def _assess(taxonomy, synth, anchors_status, hits, reference_of=None, merged=None, lineage_alone=None):
    """Concordance on the topology from the given anchors, with Mytilus's reference lineage that of a mussel."""
    genera = GENERA
    newick, placed = genus_topology_newick(synth, taxonomy, genera, anchors_status)
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True)
    matches = {g: {**taxonomy.genus[g], "status": placed[g]} for g in genera}
    reference = {"Mytilus": ("col", {"family": "Mytilidae", "class": "Bivalvia", "phylum": "Mollusca"}),
                 "Ostrea": ("col", {"family": "Ostreidae", "class": "Bivalvia", "phylum": "Mollusca"}),
                 "Magallana": ("col", {"family": "Ostreidae", "class": "Bivalvia", "phylum": "Mollusca"}),
                 "Asterias": ("col", {"family": "Asteriidae", "class": "Asteroidea", "phylum": "Echinodermata"}),
                 **(reference_of or {})}
    phylum = {g: ref["phylum"] for g, (_, ref) in reference.items()}
    gastropod = [{"name": "Mollusca", "rank": "phylum"}, {"name": "Bivalvia", "rank": "class"}]
    starfish = [{"name": "Echinodermata", "rank": "phylum"}, {"name": "Asteroidea", "rank": "class"}]
    merged = {"Mytilus": ["6548", "999"]} if merged is None else merged
    names = cc.name_index(taxonomy)
    from timegraft.otol_static import lineage_names

    rows = cc.assess(taxonomy, names, genera, matches, reference, merged, {"6548": gastropod, "999": starfish},
                     hits, cc.backbone_sisters(tree), phylum, lineage_names, flag_level=cc.LEVELS["class"],
                     max_share=0.0, min_hits=5, min_flags=2, lineage_alone=lineage_alone)
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
    # Mytilidae is not in this synthesis: the genus goes to the enclosing node and is not placed within it (#70)
    assert row["reanchor_how"] == "reference_lineage;ancestor_in_synthesis;enclosing_node"
    assert row["reanchor_species_inside"] == ""
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


@pytest.mark.parametrize("alone", [None, "class"])
def test_wrong_name_only_is_moved_on_the_lineage_check(tmp_path, alone):
    # OTT gave the name of the run's Asterias to a starfish, without merging any NCBI ids, and there are no
    # barcode neighbours: only the lineage check can fail. With lineage_alone it flags and moves the genus (#71)
    taxonomy = Taxonomy(DATA / "taxonomy.tsv", set(GENERA), context="Metazoa")
    path = tmp_path / "synth.tre"
    path.write_text("(((ott10)ott9,ott12)ott3,((ott15)ott14,(ott6,ott7)ott5)ott13)ott2;")
    synth = SynthTree(path)
    anchors, status = genus_anchors(synth, taxonomy, GENERA)
    mussel = {"family": "Mytilidae", "class": "Bivalvia", "phylum": "Mollusca"}
    rows, matches, names, reference, _ = _assess(
        taxonomy, synth, (anchors, status), {}, reference_of={"Asterias": ("records", mussel)}, merged={},
        lineage_alone=None if alone is None else cc.LEVELS[alone])
    by = {r["genus"]: r for r in rows}
    assert by["Asterias"]["lineage_flag"] == "yes" and by["Asterias"]["n_flags"] == "1"
    provenance, moved = ra.plan(synth, taxonomy, names, matches, anchors, rows,
                                {g: ref for g, (_, ref) in reference.items()}, {}, "apply", cc.LEVELS["class"])
    if alone is None:  # the earlier rule: listed for review, left in place
        assert by["Asterias"]["concordance_flag"] == "review" and by["Asterias"]["flag_rule"] == ""
        assert moved == {}
        return
    assert by["Asterias"]["concordance_flag"] == "flag" and by["Asterias"]["flag_rule"] == "lineage_alone"
    assert [r["genus"] for r in rows if r["concordance_flag"] == "flag"] == ["Asterias"]
    assert provenance["Asterias"]["placement_action"] == "reanchored"
    assert provenance["Asterias"]["reanchor_target"] == "family:Mytilidae (ott4)"


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


# A resolved family (#70). OTT merged the mussel Mytilus into a starfish genus of the same name: the merged
# taxon sits in Asteroidea and holds the mussels' species (6, 7, 8) and more starfish species (30, 31, 32 and
# 33). The synthesis, built from trees, puts the mussels inside Mytilidae, sister to Perna.
RESOLVED_TAXONOMY = """uid	|	parent_uid	|	name	|	rank	|	sourceinfo	|	uniqname	|	flags	|
1	|		|	life	|	no rank	|	ncbi:1	|		|		|
2	|	1	|	Metazoa	|	kingdom	|	ncbi:33208	|		|		|
3	|	2	|	Bivalvia	|	class	|	ncbi:6544	|		|		|
4	|	3	|	Mytilidae	|	family	|	ncbi:6545	|		|		|
20	|	4	|	Perna	|	genus	|	ncbi:6549	|		|		|
21	|	20	|	Perna viridis	|	species	|	ncbi:73031	|		|		|
22	|	20	|	Perna perna	|	species	|	ncbi:73032	|		|		|
23	|	4	|	Modiolus	|	genus	|	ncbi:6560	|		|		|
24	|	23	|	Modiolus modiolus	|	species	|	ncbi:6561	|		|		|
13	|	2	|	Asteroidea	|	class	|	ncbi:7588	|		|		|
5	|	13	|	Mytilus	|	genus	|	ncbi:6548	|		|		|
6	|	5	|	Mytilus edulis	|	species	|	ncbi:6550	|		|		|
7	|	5	|	Mytilus trossulus	|	species	|	ncbi:6551	|		|		|
8	|	5	|	Mytilus galloprovincialis	|	species	|	ncbi:29158	|		|		|
30	|	5	|	Mytilus stella	|	species	|		|		|		|
31	|	5	|	Mytilus astra	|	species	|		|		|		|
32	|	5	|	Mytilus radiata	|	species	|		|		|		|
33	|	5	|	Mytilus spina	|	species	|		|		|		|
14	|	13	|	Asterias	|	genus	|	ncbi:7601	|		|		|
15	|	14	|	Asterias rubens	|	species	|	ncbi:7604	|		|		|
"""
RESOLVED_SYNTH = ("((((((ott6,ott7)mrcaott6ott7,ott8)mrcaott6ott8,(ott21,ott22)ott20)mrcaott6ott21,(ott24)ott23)ott4)ott3,"
                  "(((ott30,ott31,ott32,ott33)starfish,(ott15)ott14)ott13)ott12)ott2;")
RESOLVED_GENERA = ["Mytilus", "Perna", "Modiolus", "Asterias"]


def _resolved(tmp_path, species=(6, 7, 8)):
    text = "\n".join(line for line in RESOLVED_TAXONOMY.splitlines()
                     if not (line.split("\t|\t")[0].isdigit() and int(line.split("\t|\t")[0]) in {6, 7, 8} - set(species)))
    (tmp_path / "taxonomy.tsv").write_text(text + "\n")
    (tmp_path / "synth.tre").write_text(RESOLVED_SYNTH)
    taxonomy = Taxonomy(tmp_path / "taxonomy.tsv", set(RESOLVED_GENERA), context="Metazoa")
    synth = SynthTree(tmp_path / "synth.tre")
    anchors, status = genus_anchors(synth, taxonomy, RESOLVED_GENERA, rule="majority")
    return taxonomy, synth, anchors, status


def _move(taxonomy, synth, anchors, status, within):
    matches = {g: {**taxonomy.genus[g], "status": status[g]} for g in RESOLVED_GENERA}
    flag = [{"genus": "Mytilus", "concordance_flag": "flag", "severity": str(cc.LEVELS["phylum"])}]
    mussel = {"Mytilus": {"family": "Mytilidae", "class": "Bivalvia", "phylum": "Mollusca"}}
    provenance, moved = ra.plan(synth, taxonomy, cc.name_index(taxonomy), matches, anchors, flag, mussel, {},
                                "apply", cc.LEVELS["class"], within=within)
    anchors, status = dict(anchors), dict(status)
    anchors["Mytilus"], status["Mytilus"] = moved["Mytilus"][0], f"reanchored:{synth.label[moved['Mytilus'][0]]}"
    newick, _ = genus_topology_newick(synth, taxonomy, RESOLVED_GENERA, (anchors, status))
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True)
    return provenance["Mytilus"], synth.label[moved["Mytilus"][0]], cc.backbone_sisters(tree)


def test_reanchored_within_a_resolved_family(tmp_path):
    taxonomy, synth, anchors, status = _resolved(tmp_path)
    # the starfish species are the majority of the merged taxon: anchored among the starfish
    assert synth.label[anchors["Mytilus"]] == "starfish"
    row, node, sisters = _move(taxonomy, synth, anchors, status, "majority")
    # attached within Mytilidae at the majority clade of its three species there (two of them, as #51 has it),
    # so Perna keeps Mytilus as its sister
    assert node == "mrcaott6ott7" and row["reanchor_species_inside"] == "3"
    assert row["reanchor_how"] == "reference_lineage;taxon_in_synthesis;majority_clade_in_target"
    assert sisters["Perna"] == ["Mytilus"] and sisters["Mytilus"] == ["Perna"]
    assert set(sisters["Modiolus"]) == {"Mytilus", "Perna"}


def test_reanchored_at_the_family_node_as_before(tmp_path):
    taxonomy, synth, anchors, status = _resolved(tmp_path)
    row, node, sisters = _move(taxonomy, synth, anchors, status, "off")
    # within_target off: at the family's node, the genus's stem becomes the family's crown, and Perna loses
    # Mytilus as its only sister: the family becomes a polytomy of the three genera
    assert node == "ott4" and row["reanchor_species_inside"] == ""
    assert row["reanchor_how"] == "reference_lineage;taxon_in_synthesis"
    assert set(sisters["Mytilus"]) == {"Modiolus", "Perna"} and set(sisters["Perna"]) == {"Modiolus", "Mytilus"}


@pytest.mark.parametrize("species,node,how", [((6, 7), "mrcaott6ott7", "mrca_in_target"),
                                              ((8,), "mrcaott6ott8", "species_in_target")])
def test_few_species_inside_the_target(tmp_path, species, node, how):
    # fewer than min_species inside the target: their MRCA, or the parent of a single one, as for anchoring (#51)
    taxonomy, synth, anchors, status = _resolved(tmp_path, species)
    row, moved_to, _ = _move(taxonomy, synth, anchors, status, "majority")
    assert moved_to == node and row["reanchor_how"].endswith(how)
    assert row["reanchor_species_inside"] == str(len(species))


def test_within_target_prefers_the_sampled_species(tmp_path):
    taxonomy, synth, *_ = _resolved(tmp_path)
    target = synth.index["ott4"]
    node, how, n = ra.within_target(synth, taxonomy, "Mytilus", target, sampled={"Mytilus edulis", "Mytilus trossulus"})
    assert synth.label[node] == "mrcaott6ott7" and how == "mrca_in_target" and n == 2
    # one sampled species is too few (min_sampled 2): all three species inside the target
    node, how, n = ra.within_target(synth, taxonomy, "Mytilus", target, sampled={"Mytilus edulis"})
    assert how == "majority_clade_in_target" and n == 3
    # the starfish species of the merged taxon are never counted: they are outside the target
    node, how, n = ra.within_target(synth, taxonomy, "Mytilus", synth.index["ott13"])
    assert how == "majority_clade_in_target" and synth.label[node] == "starfish" and n == 4


def test_within_target_must_be_known():
    with pytest.raises(ValueError, match="within_target"):
        ra.plan(None, None, {}, {}, {}, [], {}, {}, "apply", 3, within="family")


def test_target_not_in_the_synthesis_is_not_searched(tmp_path):
    # Mytilidae has no node and no members in the synthesis: the mussels hang directly below Bivalvia. The genus
    # goes to the enclosing node, as before #70, and is not placed among the mussel species inside it
    taxonomy, synth, anchors, status = _resolved(tmp_path)
    (tmp_path / "synth.tre").write_text("(((ott6,ott7)mrcaott6ott7,ott8,(ott50,ott51)other)ott3,"
                                        "(((ott30,ott31,ott32,ott33)starfish,(ott15)ott14)ott13)ott12)ott2;")
    synth = SynthTree(tmp_path / "synth.tre")
    for g in ("Perna", "Modiolus"):
        anchors.pop(g, None)
    anchors.update(genus_anchors(synth, taxonomy, ["Mytilus", "Asterias"], rule="majority")[0])
    matches = {g: {**taxonomy.genus[g], "status": status[g]} for g in ("Mytilus", "Asterias")}
    flag = [{"genus": "Mytilus", "concordance_flag": "flag", "severity": str(cc.LEVELS["phylum"])}]
    mussel = {"Mytilus": {"family": "Mytilidae", "class": "Bivalvia", "phylum": "Mollusca"}}
    provenance, moved = ra.plan(synth, taxonomy, cc.name_index(taxonomy), matches, anchors, flag, mussel, {},
                                "apply", cc.LEVELS["class"])
    row = provenance["Mytilus"]
    assert synth.label[moved["Mytilus"][0]] == "ott3"
    assert row["reanchor_how"] == "reference_lineage;ancestor_in_synthesis;enclosing_node"
    assert row["reanchor_species_inside"] == ""


# Only the name moved (#101). OTT gave the name Mytilus to a starfish genus (5) without the mussels' species:
# those stay in Mytilidae under the mussels' own taxon, renamed (40), as the `name` variant of inject.py does.
NAME_ONLY_TAXONOMY = "\n".join(
    line.replace("|\t5\t|\tMytilus edulis", "|\t40\t|\tMytilus edulis")
        .replace("|\t5\t|\tMytilus trossulus", "|\t40\t|\tMytilus trossulus")
        .replace("|\t5\t|\tMytilus galloprovincialis", "|\t40\t|\tMytilus galloprovincialis")
    for line in RESOLVED_TAXONOMY.splitlines()) + "\n40\t|\t4\t|\tMytilus zzmerged\t|\tgenus\t|\tncbi:6548\t|\t\t|\t\t|\n"


def _name_only(tmp_path, by_name=True):
    (tmp_path / "taxonomy.tsv").write_text(NAME_ONLY_TAXONOMY)
    (tmp_path / "synth.tre").write_text(RESOLVED_SYNTH)
    taxonomy = Taxonomy(tmp_path / "taxonomy.tsv", set(RESOLVED_GENERA), context="Metazoa", species_by_name=by_name)
    synth = SynthTree(tmp_path / "synth.tre")
    anchors, status = genus_anchors(synth, taxonomy, RESOLVED_GENERA, rule="majority")
    return taxonomy, synth, anchors, status


def test_name_only_species_are_indexed(tmp_path):
    taxonomy, *_ = _name_only(tmp_path)
    # the matched taxon holds only the starfish; the mussels are found by name, and nothing else is
    assert {uid for uid, _ in taxonomy.species["Mytilus"]} == {30, 31, 32, 33}
    assert sorted(taxonomy.species_named["Mytilus"]) == [(6, "Mytilus edulis", 6550), (7, "Mytilus trossulus", 6551),
                                                         (8, "Mytilus galloprovincialis", 29158)]
    assert set(taxonomy.species_named) == {"Mytilus"}
    assert not _name_only(tmp_path, by_name=False)[0].species_named


def test_name_only_genus_is_attached_within_the_target(tmp_path):
    # the genus's species are found by name inside Mytilidae and the genus goes to its own clade: Perna keeps
    # Mytilus as its sister, as in the merge case of #70
    row, node, sisters = _move(*_name_only(tmp_path), "majority")
    assert node == "mrcaott6ott7" and row["reanchor_species_inside"] == "3"
    assert row["reanchor_how"] == "reference_lineage;taxon_in_synthesis;majority_clade_in_target;by_name"
    assert sisters["Perna"] == ["Mytilus"] and sisters["Mytilus"] == ["Perna"]


def test_name_only_without_the_name_search(tmp_path):
    # without the search the genus goes to the family's node, as before #101, and Perna loses its sister
    row, node, sisters = _move(*_name_only(tmp_path, by_name=False), "majority")
    assert node == "ott4" and row["reanchor_how"].endswith(";target_node")
    assert set(sisters["Perna"]) == {"Modiolus", "Mytilus"}


def test_name_search_stays_within_the_target(tmp_path):
    taxonomy, synth, *_ = _name_only(tmp_path)
    # species named after the genus but outside the target taxon in OTT are not taken, even if the synthesis
    # places them inside the target node
    node, how, n = ra.within_target(synth, taxonomy, "Mytilus", synth.index["ott4"], target_uid=13)
    assert how == "target_node" and n == 0
    # the sampled species are preferred among those found by name
    node, how, n = ra.within_target(synth, taxonomy, "Mytilus", synth.index["ott4"], target_uid=4,
                                    sampled={"Mytilus edulis", "Mytilus trossulus"})
    assert synth.label[node] == "mrcaott6ott7" and how == "mrca_in_target;by_name" and n == 2
    # where the matched taxon has species inside the target, the name search is not made
    node, how, n = ra.within_target(synth, taxonomy, "Mytilus", synth.index["ott13"], target_uid=13)
    assert how == "majority_clade_in_target" and n == 4


def test_name_only_genus_is_dated_with_its_own_species(tmp_path):
    from timegraft.otol_static import lookup_ids

    taxonomy, *_ = _name_only(tmp_path)
    # the matched taxon has no NCBI ids of the mussels' species; those found by name inside Mytilidae do
    assert lookup_ids(taxonomy, "Mytilus") == [6548]
    own = [(u, n) for u, _, n in ra.named_species(taxonomy, "Mytilus", 4)]
    assert lookup_ids(taxonomy, "Mytilus", species=own) == [6548, 6550, 6551, 29158]
    assert lookup_ids(taxonomy, "Mytilus", species=own, keep=lambda n: n != 6551) == [6548, 6550, 29158]
    assert ra.named_species(taxonomy, "Mytilus", 13) == []


def test_name_only_ids_without_the_identifier_check(tmp_path):
    from timegraft.otol_static import lookup_ids

    taxonomy, *_ = _name_only(tmp_path)
    # with check_ids: false the species ids are still those found by name; the genus id is kept as it is
    own = [(u, n) for u, _, n in ra.named_species(taxonomy, "Mytilus", 4)]
    species = lookup_ids(taxonomy, "Mytilus", species=own)[1:]
    assert ra.with_species(lookup_ids(taxonomy, "Mytilus"), species) == [6548, 6550, 6551, 29158]
    assert ra.with_species([6548, 999], species) == [6548, 6550, 6551, 29158]
    assert ra.with_species(6548, species) == [6548, 6550, 6551, 29158]
    assert ra.with_species([None], []) == [None]
