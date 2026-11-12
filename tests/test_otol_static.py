"""Genus topology from static OToL files, on a tiny hand-made taxonomy and tree."""

from pathlib import Path

import dendropy
import pytest

from timegraft.backbone_build import _higher_sides, date_topology
from timegraft.otol_static import (SynthTree, Taxonomy, clade_richness, genus_anchors, genus_topology_newick, lineage,
                                   lookup_ids)
from timegraft.timetree import TimeTree

DATA = Path(__file__).parent / "data" / "otol"
GENERA = ["Mytilus", "Ostrea", "Magallana", "Asterias", "Lostgenus", "Nowhere"]


@pytest.fixture(scope="module")
def taxonomy():
    return Taxonomy(DATA / "taxonomy.tsv", set(GENERA), context="Metazoa")


def test_context_resolves_homonyms(taxonomy):
    assert taxonomy.genus["Mytilus"] == {"ott_id": 5, "ncbi_id": 6548, "status": "matched"}
    assert taxonomy.genus["Nowhere"]["status"] == "unmatched"
    assert Taxonomy(DATA / "taxonomy.tsv", {"Mytilus"}).genus["Mytilus"]["status"] == "ambiguous"


def test_species_and_lineage(taxonomy):
    assert sorted(n for _, n in taxonomy.species["Mytilus"]) == [6550, 6551]
    assert lookup_ids(taxonomy, "Mytilus") == [6548, 6550, 6551]
    assert lookup_ids(taxonomy, "Magallana") == [0, 29159]  # genus without an NCBI id
    # only ranks from family to phylum carry NCBI ids here; kingdom and above are too broad
    assert [n for _, n in lineage(taxonomy, "Asterias") if n] == [7588]


def test_genus_topology(taxonomy):
    synth = SynthTree(DATA / "labelled_supertree.tre")
    newick, status = genus_topology_newick(synth, taxonomy, GENERA)
    assert status["Mytilus"] == "in_synthesis"
    assert status["Magallana"] == "attached_at_mrca"
    assert status["Lostgenus"].startswith("attached_at_taxon:")  # via Mytilidae, which is not in the tree either
    assert status["Nowhere"] == "unmatched"
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True)
    assert {l.taxon.label for l in tree.leaf_node_iter()} == {"Mytilus", "Ostrea", "Magallana", "Asterias", "Lostgenus"}
    asterias = tree.find_node_with_taxon_label("Asterias")
    assert asterias.parent_node is tree.seed_node


def test_higher_taxa_stand_in_only_when_they_stay_on_one_side(taxonomy):
    tree = dendropy.Tree.get(data="((Mytilus,Ostrea),Asterias);", schema="newick", rooting="force-rooted")
    higher = {g: lineage(taxonomy, g) for g in ["Mytilus", "Ostrea", "Asterias"]}
    sides = _higher_sides(tree.seed_node, higher, limit=10)
    # Bivalvia (6544) for one side, Asteroidea (7588) for the other; Metazoa is shared, so never used
    assert sorted(map(set, sides), key=len) == [{7588}, {6545, 6563, 6544}]

    def fetch(url, payload=None):
        ids = [int(i) for i in url.split("/mrca/id/")[1].split("/")[0].split("+")]
        known = {6544, 7588}  # TimeTree knows the classes, not the genera
        found = [i for i in ids if i in known]
        return {"found_ids": found, "precomputed_age": 600.0 if len(found) == 2 else None}

    ncbi = {g: lookup_ids(taxonomy, g) for g in ["Mytilus", "Ostrea", "Asterias"]}
    tree, report = date_topology(tree, ncbi, TimeTree(fetch=fetch), higher=higher, root_age=None)
    root = next(r for r in report if r["n_genera"] == 3)
    assert root["source"] == "timetree_higher_taxa" and float(root["age"]) == pytest.approx(600.0)


def test_clade_richness_counts_genera_missing_from_the_data():
    genera = ["Mytilus", "Ostrea", "Asterias"]
    taxonomy = Taxonomy(DATA / "taxonomy.tsv", set(genera), context="Metazoa")
    synth = SynthTree(DATA / "labelled_supertree.tre")
    newick, _ = genus_topology_newick(synth, taxonomy, genera)
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True)
    anchor, _ = genus_anchors(synth, taxonomy, genera)
    bivalves = tree.find_node_with_taxon_label("Mytilus").parent_node
    # Magallana is not in the data but is in the bivalve clade of the synthesis
    assert clade_richness(synth, anchor, tree, taxonomy) == {bivalves: 3, tree.seed_node: 4}
    assert clade_richness(synth, anchor, tree) == {bivalves: 4, tree.seed_node: 5}
    # subtree ranges follow from the preorder numbering
    sizes = synth.subtree_sizes()
    assert sizes[0] == len(synth.parent)
