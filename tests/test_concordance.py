"""Placement concordance diagnostics (#49): report only."""

from pathlib import Path

import dendropy
import pytest

from timegraft import concordance as cc
from timegraft.otol_static import Taxonomy

DATA = Path(__file__).parent / "data" / "otol"


def test_split_level():
    genus = [("g", ""), ("sub", "subfamily"), ("fam", "family"), ("ord", "order"), ("cls", "class"), ("root", "")]
    assert cc.split_level(genus, [("fam", "family"), ("ord", "order")]) == (0, "fam")  # family contains the genus
    other_family = [("fam2", "family"), ("ord", "order"), ("cls", "class")]
    assert cc.split_level(genus, other_family) == (1, "ord")
    # parting at an unranked node between two phyla is still a phylum-level split
    a = [("x", ""), ("p1", "phylum"), ("bilateria", ""), ("metazoa", "kingdom")]
    b = [("y", "family"), ("p2", "phylum"), ("bilateria", ""), ("metazoa", "kingdom")]
    assert cc.split_level(a, b) == (5, "bilateria")


@pytest.fixture(scope="module")
def taxonomy():
    return Taxonomy(DATA / "taxonomy.tsv", {"Mytilus", "Ostrea", "Asterias"}, context="Metazoa")


def test_lineage_distance(taxonomy):
    names = cc.name_index(taxonomy)
    mytilus = taxonomy.genus["Mytilus"]["ott_id"]
    assert cc.lineage_distance(taxonomy, names, mytilus, {"family": "Mytilidae"})["lineage_distance"] == "0"
    # another bivalve family: they part below Bivalvia
    row = cc.lineage_distance(taxonomy, names, mytilus, {"family": "Ostreidae"})
    assert row == {"reference_taxon": "family:Ostreidae", "lineage_mrca": "Bivalvia", "lineage_split": "family",
                   "lineage_distance": "1"}
    # a starfish class for a mussel: different classes
    row = cc.lineage_distance(taxonomy, names, mytilus, {"family": "Nowhereidae", "class": "Asteroidea"})
    assert row["reference_taxon"] == "class:Asteroidea" and row["lineage_split"] == "class"
    assert row["lineage_mrca"] == "Metazoa"
    assert cc.lineage_distance(taxonomy, names, mytilus, {})["lineage_distance"] == ""


def _ncbi(*names_ranks):
    return [{"name": n, "rank": r} for n, r in names_ranks]


def test_merged_id_distance():
    gastropod = _ncbi(("Metazoa", "kingdom"), ("Bilateria", "clade"), ("Mollusca", "phylum"),
                      ("Gastropoda", "class"), ("Onchidiidae", "family"))
    echinoid = _ncbi(("Metazoa", "kingdom"), ("Bilateria", "clade"), ("Echinodermata", "phylum"),
                     ("Echinoidea", "class"), ("Laganidae", "family"))
    synonym = _ncbi(("Metazoa", "kingdom"), ("Mollusca", "phylum"), ("Gastropoda", "class"), ("Onchidiidae", "family"))
    row = cc.merged_id_distance({"36932": gastropod, "262330": echinoid, "1": synonym}, ["262330", "36932", "1"])
    assert row == {"merged_mrca": "Bilateria", "merged_split": "phylum", "merged_id_distance": "5"}
    assert cc.merged_id_distance({"36932": gastropod, "1": synonym}, ["36932", "1"])["merged_id_distance"] == "0"
    assert cc.merged_id_distance({"36932": gastropod}, ["36932", "2"])["merged_id_distance"] == ""
    assert cc.merged_id_distance({}, ["36932"])["merged_id_distance"] == ""


def test_neighbours_and_sisters():
    tree = dendropy.Tree.get(data="((Onchidella,Arbacia),(Limax,Oxychilus));", schema="newick")
    sisters = cc.backbone_sisters(tree)
    assert sisters["Onchidella"] == ["Arbacia"]
    phylum = {"Onchidella": "Mollusca", "Arbacia": "Echinodermata", "Limax": "Mollusca", "Oxychilus": "Mollusca"}
    row = cc.neighbour_concordance("Onchidella", ["Limax"] * 3 + ["Oxychilus"] * 2, sisters["Onchidella"], phylum)
    assert row["neighbour_concordance"] == "0.000" and row["neighbour_phyla"] == "Mollusca:5"
    assert row["sister_phyla"] == "Echinodermata"
    row = cc.neighbour_concordance("Limax", ["Oxychilus"] * 4, sisters["Limax"], phylum)
    assert row["neighbour_concordance"] == "" and row["neighbour_hits"] == "4"  # too few hits


def test_aggregate_and_ranking():
    level = cc.LEVELS["class"]
    one = cc.aggregate({"genus": "A", "lineage_distance": "5", "merged_id_distance": "5",
                        "neighbour_concordance": "0.000"}, level, 0.0)
    assert one["n_flags"] == "3" and one["concordance_flag"] == "flag" and one["severity"] == "5"
    two = cc.aggregate({"genus": "B", "lineage_distance": "0", "merged_id_distance": "",
                        "neighbour_concordance": "0.000"}, level, 0.0)
    assert two["concordance_flag"] == "review" and two["neighbour_flag"] == "yes" and two["lineage_flag"] == ""
    three = cc.aggregate({"genus": "C", "lineage_distance": "3", "merged_id_distance": "",
                          "neighbour_concordance": "0.900"}, level, 0.0)
    assert three["concordance_flag"] == "" and three["n_evaluable"] == "2"
    four = cc.aggregate({"genus": "D", "lineage_distance": "", "merged_id_distance": "",
                         "neighbour_concordance": ""}, level, 0.0)
    assert four["concordance_flag"] == "not_evaluable"
    ranked = cc.rank_rows([four, three, two, one])
    assert [r["genus"] for r in ranked] == ["A", "B", "C", "D"] and ranked[0]["rank"] == "1"
