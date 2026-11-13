"""Majority-clade anchoring of genera that are not a node of the synthesis (#51)."""

import pytest

from timegraft.otol_static import ANCHOR_COLUMNS, SynthTree, Taxonomy, genus_anchors, majority_clade

# Genus Testa (ott100) is not a node of the synthesis. Four of its five species form a clade within
# family A; the fifth (ott105) has been placed in family B, far away, so the MRCA of all five is the root.
TAXONOMY = """uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t
1\t|\t\t|\tlife\t|\tno rank\t|\t\t|\t\t|\t\t|\t
2\t|\t1\t|\tMetazoa\t|\tkingdom\t|\tncbi:33208\t|\t\t|\t\t|\t
3\t|\t2\t|\tAidae\t|\tfamily\t|\tncbi:3\t|\t\t|\t\t|\t
100\t|\t3\t|\tTesta\t|\tgenus\t|\tncbi:100\t|\t\t|\t\t|\t
101\t|\t100\t|\tTesta alpha\t|\tspecies\t|\tncbi:101\t|\t\t|\t\t|\t
102\t|\t100\t|\tTesta beta\t|\tspecies\t|\tncbi:102\t|\t\t|\t\t|\t
103\t|\t100\t|\tTesta gamma\t|\tspecies\t|\tncbi:103\t|\t\t|\t\t|\t
104\t|\t100\t|\tTesta delta\t|\tspecies\t|\tncbi:104\t|\t\t|\t\t|\t
105\t|\t100\t|\tTesta stray\t|\tspecies\t|\tncbi:105\t|\t\t|\t\t|\t
110\t|\t3\t|\tOther\t|\tgenus\t|\tncbi:110\t|\t\t|\t\t|\t
111\t|\t110\t|\tOther one\t|\tspecies\t|\tncbi:111\t|\t\t|\t\t|\t
4\t|\t2\t|\tBidae\t|\tfamily\t|\tncbi:4\t|\t\t|\t\t|\t
120\t|\t4\t|\tFar\t|\tgenus\t|\tncbi:120\t|\t\t|\t\t|\t
121\t|\t120\t|\tFar one\t|\tspecies\t|\tncbi:121\t|\t\t|\t\t|\t
130\t|\t3\t|\tDuo\t|\tgenus\t|\tncbi:130\t|\t\t|\t\t|\t
131\t|\t130\t|\tDuo one\t|\tspecies\t|\tncbi:131\t|\t\t|\t\t|\t
132\t|\t130\t|\tDuo two\t|\tspecies\t|\tncbi:132\t|\t\t|\t\t|\t
"""
SYNTH = ("((((ott101,ott102)mrcaott101ott102,(ott103,ott104)mrcaott103ott104)mrcaott101ott103,ott131,ott111)ott3,"
         "((ott121)ott120,ott105,ott132)ott4)ott2;")
GENERA = ["Testa", "Other", "Far", "Duo"]


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    d = tmp_path_factory.mktemp("anchoring")
    (d / "taxonomy.tsv").write_text(TAXONOMY)
    (d / "synth.tre").write_text(SYNTH)
    return Taxonomy(d / "taxonomy.tsv", set(GENERA), context="Metazoa"), SynthTree(d / "synth.tre")


def test_majority_clade(data):
    _, synth = data
    nodes = [synth.index[f"ott{u}"] for u in (101, 102, 103, 104, 105)]
    node, inside = majority_clade(synth, nodes)
    assert synth.label[node] == "mrcaott101ott103" and inside == 4
    # with a stricter share no smaller clade qualifies than the root
    node, inside = majority_clade(synth, nodes, share=0.8)
    assert synth.label[node] == "ott2" and inside == 5


def test_one_stray_species_no_longer_drags_the_anchor(data):
    taxonomy, synth = data
    old, old_status = genus_anchors(synth, taxonomy, GENERA, rule="mrca")
    assert synth.label[old["Testa"]] == "ott2" and old_status["Testa"] == "attached_at_mrca"
    diagnostics = {}
    new, status = genus_anchors(synth, taxonomy, GENERA, rule="majority", diagnostics=diagnostics)
    assert synth.label[new["Testa"]] == "mrcaott101ott103" and status["Testa"] == "attached_at_majority"
    d = diagnostics["Testa"]
    assert set(d) == set(ANCHOR_COLUMNS)
    assert d["anchor_rule"] == "majority" and d["anchor_basis"] == "all" and d["anchor_fallback"] == ""
    assert (d["anchor_species"], d["anchor_species_inside"], d["anchor_outside_share"]) == ("5", "4", "0.200")
    assert d["anchor_outside_species"] == "Testa stray"
    assert d["anchor_mrca_node"] == "ott2" and d["anchor_mrca_clade_genera"] == "4"
    # genera that are nodes of the synthesis are untouched
    assert status["Far"] == "in_synthesis" and new["Far"] == old["Far"]


def test_sampled_species_are_preferred(data):
    taxonomy, synth = data
    diagnostics = {}
    sampled = {"Testa": {"Testa gamma", "Testa delta", "Testa nova"}}
    anchors, status = genus_anchors(synth, taxonomy, GENERA, rule="majority", sampled=sampled,
                                    diagnostics=diagnostics)
    # only two sampled species in the synthesis: their MRCA, too few for a majority
    assert synth.label[anchors["Testa"]] == "mrcaott103ott104" and status["Testa"] == "attached_at_mrca"
    assert diagnostics["Testa"]["anchor_basis"] == "sampled"
    assert diagnostics["Testa"]["anchor_fallback"] == "too_few_species"
    # a sampled set that includes the stray: still the majority
    sampled = {"Testa": {"Testa alpha", "Testa beta", "Testa stray"}}
    anchors, _ = genus_anchors(synth, taxonomy, GENERA, rule="majority", sampled=sampled)
    assert synth.label[anchors["Testa"]] == "mrcaott101ott102"
    # fewer than min_sampled sampled species in the synthesis: all species instead
    anchors, _ = genus_anchors(synth, taxonomy, GENERA, rule="majority", sampled={"Testa": {"Testa alpha"}},
                               diagnostics=(diag := {}))
    assert synth.label[anchors["Testa"]] == "mrcaott101ott103" and diag["Testa"]["anchor_basis"] == "all"


def test_controlled_fallbacks(data):
    taxonomy, synth = data
    diagnostics = {}
    # Duo's two species are split between families: two species are too few for a majority
    anchors, status = genus_anchors(synth, taxonomy, GENERA, rule="majority", diagnostics=diagnostics)
    assert synth.label[anchors["Duo"]] == "ott2" and status["Duo"] == "attached_at_mrca"
    assert diagnostics["Duo"]["anchor_fallback"] == "too_few_species"
    assert diagnostics["Duo"]["anchor_outside_share"] == "0.000"
    # no clade holds more than 80 % of Testa's species
    anchors, status = genus_anchors(synth, taxonomy, GENERA, rule="majority", share=0.8, diagnostics=diagnostics)
    assert status["Testa"] == "attached_at_mrca" and diagnostics["Testa"]["anchor_fallback"] == "no_smaller_majority"
    with pytest.raises(ValueError):
        genus_anchors(synth, taxonomy, GENERA, rule="median")
