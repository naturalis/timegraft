"""NCBI ids of re-anchored genera: only those of the right animal reach TimeTree (#69).

The case is that of *Onchidella*: OTT merged a gastropod genus and a sand dollar
genus of the same name, so its taxon holds species, and NCBI ids, of both. It sits
among the sea urchins; re-anchoring moves it to the gastropod family Onchidiidae.
"""

import dendropy
import pytest

from timegraft import concordance as cc
from timegraft import reanchor as ra
from timegraft.backbone_build import date_topology
from timegraft.otol_static import Taxonomy, crown_ages, lineage, lookup_ids
from timegraft.timetree import TimeTree

ROWS = [
    (1, "", "life", "no rank", "ncbi:1"),
    (2, 1, "Metazoa", "kingdom", "ncbi:33208"),
    (20, 2, "Mollusca", "phylum", "ncbi:6447"),
    (21, 20, "Gastropoda", "class", "ncbi:6448"),
    (22, 21, "Onchidiidae", "family", "ncbi:7000"),
    (23, 21, "Siphonariidae", "family", "ncbi:7001"),
    (24, 23, "Siphonaria", "genus", "ncbi:3000"),
    (25, 24, "Siphonaria pectinata", "species", "ncbi:3001"),
    (26, 24, "Siphonaria lessonii", "species", "ncbi:3002"),
    (30, 2, "Echinodermata", "phylum", "ncbi:7586"),
    (31, 30, "Echinoidea", "class", "ncbi:7625"),
    (32, 31, "Laganidae", "family", "ncbi:7700"),
    (33, 31, "Echinocardium", "genus", "ncbi:4000"),
    (34, 33, "Echinocardium cordatum", "species", "ncbi:4001"),
    (35, 33, "Echinocardium flavescens", "species", "ncbi:4002"),
    # the merged genus, misplaced among the urchins; the sand dollar's NCBI id first, as for Onchidella in OTT
    (40, 32, "Onchidella", "genus", "ncbi:222,ncbi:111"),
    (41, 40, "Onchidella celtica", "species", "ncbi:1001"),
    (42, 40, "Onchidella floridana", "species", "ncbi:1002"),
    (43, 40, "Onchidella nigricans", "species", "ncbi:1003"),  # NCBI puts it in another gastropod family
    (44, 40, "Onchidella marginata", "species", "ncbi:1004"),  # NCBI does not classify it
    (45, 40, "Onchidella evelinae", "species", ""),
    (46, 40, "Peronella japonica", "species", "ncbi:2001"),
    (47, 40, "Peronella lesueuri", "species", "ncbi:2002"),
]
GENERA = {"Onchidella", "Siphonaria", "Echinocardium"}

TOP = [{"name": "Metazoa", "rank": "kingdom"}]
GASTROPOD = TOP + [{"name": "Mollusca", "rank": "phylum"}, {"name": "Gastropoda", "rank": "class"}]
URCHIN = TOP + [{"name": "Echinodermata", "rank": "phylum"}, {"name": "Echinoidea", "rank": "class"}]
CLASSIFICATIONS = {
    7000: GASTROPOD,                                                    # Onchidiidae
    111: GASTROPOD + [{"name": "Onchidiidae", "rank": "family"}],       # the gastropod genus
    1001: GASTROPOD + [{"name": "Onchidiidae", "rank": "family"}],
    1002: GASTROPOD + [{"name": "Onchidiidae", "rank": "family"}],
    1003: GASTROPOD + [{"name": "Veronicellidae", "rank": "family"}],   # parts at family level: kept
    222: URCHIN + [{"name": "Laganidae", "rank": "family"}],            # the sand dollar genus
    2001: URCHIN + [{"name": "Laganidae", "rank": "family"}],
    2002: URCHIN + [{"name": "Laganidae", "rank": "family"}],
}
SAND_DOLLAR = {222, 2001, 2002}


@pytest.fixture()
def taxonomy(tmp_path):
    path = tmp_path / "taxonomy.tsv"
    lines = ["uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t"]
    lines += [f"{u}\t|\t{p}\t|\t{n}\t|\t{r}\t|\t{s}\t|\t\t|\t\t|\t" for u, p, n, r, s in ROWS]
    path.write_text("\n".join(lines) + "\n")
    return Taxonomy(path, GENERA, context="Metazoa")


class Classify:
    """ChecklistBank stand-in that counts its calls."""

    def __init__(self, table):
        self.table, self.calls = table, []

    def __call__(self, ncbi_id):
        self.calls.append(int(ncbi_id))
        return self.table.get(int(ncbi_id), [])


def recording_timetree(queries):
    """TimeTree stand-in: gastropods alone 50 Ma, with urchins 100 Ma; a query that mixes gastropods and
    sand dollars gets 600 Ma, the age of their split."""

    def fetch(url, payload=None):
        ids = set(int(i) for i in url.split("/mrca/id/")[1].split("/")[0].split("+"))
        queries.append(ids)
        urchins = {4000, 4001, 4002}
        age = 600.0 if SAND_DOLLAR & ids and ids - SAND_DOLLAR - urchins else 100.0 if urchins & ids else 50.0
        return {"found_ids": sorted(ids), "precomputed_age": age}

    return TimeTree(fetch=fetch)


def screen_for(taxonomy, classify):
    return ra.IdScreen(classify, {"Onchidella": ra.id_reference(taxonomy, classify, 22)}, cc.LEVELS["class"])


def test_reference_is_the_ncbi_lineage_of_the_new_place(taxonomy):
    reference, against = ra.id_reference(taxonomy, Classify(CLASSIFICATIONS), 22)
    assert against == "family:Onchidiidae (ncbi:7000)"
    assert reference == [("Onchidiidae", "family"), ("Gastropoda", "class"), ("Mollusca", "phylum"),
                         ("Metazoa", "kingdom")]
    # a target NCBI does not classify: the nearest taxon above it that it does
    _, against = ra.id_reference(taxonomy, Classify({6448: TOP + [{"name": "Mollusca", "rank": "phylum"}]}), 22)
    assert against == "class:Gastropoda (ncbi:6448)"
    assert ra.id_reference(taxonomy, Classify({}), 22) == (None, "")


def test_only_the_right_animal_is_kept(taxonomy):
    classify = Classify(CLASSIFICATIONS)
    screen = screen_for(taxonomy, classify)
    before = lookup_ids(taxonomy, "Onchidella")
    assert before[0] == 222 and SAND_DOLLAR & set(before)  # what OTT gives: the sand dollar's ids too
    ids = screen.screen("Onchidella", before, lambda keep: lookup_ids(taxonomy, "Onchidella", keep=keep)[1:])
    assert ids == [0, 1001, 1002, 1003]  # sand dollar genus id dropped; 1004 cannot be checked
    assert screen.genera_row("Onchidella") == {"reanchor_ids_checked": "7", "reanchor_ids_dropped": "4"}
    audit = screen.audit_row("Onchidella")
    assert audit["reanchor_ids_against"] == "family:Onchidiidae (ncbi:7000)"
    assert audit["reanchor_kept_ids"] == "1001;1002;1003"
    assert audit["reanchor_dropped_ids"] == ("genus:222:parts_at_phylum;1004:unclassified;"
                                             "2001:parts_at_phylum;2002:parts_at_phylum")
    # the gastropod genus id is kept when it is the one OTT lists first
    assert screen.screen("Onchidella", [111], lambda keep: [])[0] == 111
    assert "genus:111" in screen.audit_row("Onchidella")["reanchor_kept_ids"]
    # genera that were not re-anchored are not checked
    calls = len(classify.calls)
    assert screen.keep("Echinocardium", 4001) and screen.keep("Siphonaria", 2001)
    assert len(classify.calls) == calls and "Siphonaria" not in screen


def test_threshold_is_the_reanchoring_level(taxonomy):
    classify = Classify(CLASSIFICATIONS)
    strict = ra.IdScreen(classify, {"Onchidella": ra.id_reference(taxonomy, classify, 22)}, cc.LEVELS["family"])
    assert not strict.keep("Onchidella", 1003) and strict.decisions["Onchidella"][1003] == "parts_at_family"
    assert strict.keep("Onchidella", 1001)


def test_node_queries_never_carry_the_other_animal(taxonomy):
    newick = "((Onchidella,Siphonaria),Echinocardium);"
    others = {"Siphonaria": lookup_ids(taxonomy, "Siphonaria"), "Echinocardium": lookup_ids(taxonomy, "Echinocardium")}
    higher = {g: lineage(taxonomy, g) for g in GENERA}
    higher["Onchidella"] = ra.lineage_from(taxonomy, 22)

    def run(onchidella):
        queries = []
        tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted")
        _, report = date_topology(tree, {**others, "Onchidella": onchidella}, recording_timetree(queries),
                                  per_genus=12, higher=higher, method="precedence")
        return queries, report

    # without the check, OTT's ids put the sand dollar into the gastropod node: the age of their split
    queries, _ = run(lookup_ids(taxonomy, "Onchidella"))
    assert any(SAND_DOLLAR & q and {3000, 3001, 3002} & q for q in queries)  # gastropod node, sand dollar ids
    screen = screen_for(taxonomy, Classify(CLASSIFICATIONS))
    ids = screen.screen("Onchidella", lookup_ids(taxonomy, "Onchidella"),
                        lambda keep: lookup_ids(taxonomy, "Onchidella", keep=keep)[1:])
    queries, report = run(ids)
    assert queries and not any(SAND_DOLLAR & q for q in queries)
    assert {1001, 1002, 1003} <= set().union(*queries)
    by = {r["example_genera"]: r for r in report}
    assert float(by["Onchidella;Siphonaria"]["age"]) == 50.0 and by["Onchidella;Siphonaria"]["source"] == "timetree"


def test_fallback_to_the_new_place_when_no_id_survives(taxonomy):
    # NCBI classifies only the target: every id of the genus is dropped as unclassified
    screen = screen_for(taxonomy, Classify({7000: GASTROPOD}))
    ids = screen.screen("Onchidella", lookup_ids(taxonomy, "Onchidella"),
                        lambda keep: lookup_ids(taxonomy, "Onchidella", keep=keep)[1:])
    assert ids == [0]
    assert screen.genera_row("Onchidella")["reanchor_ids_dropped"] == "7"
    queries = []
    higher = {g: lineage(taxonomy, g) for g in GENERA}
    higher["Onchidella"] = ra.lineage_from(taxonomy, 22)
    tree = dendropy.Tree.get(data="((Onchidella,Siphonaria),Echinocardium);", schema="newick",
                             rooting="force-rooted")
    _, report = date_topology(tree, {"Onchidella": ids, "Siphonaria": lookup_ids(taxonomy, "Siphonaria"),
                                     "Echinocardium": lookup_ids(taxonomy, "Echinocardium")},
                              recording_timetree(queries), higher=higher, method="precedence")
    by = {r["example_genera"]: r for r in report}
    # Onchidella is stood in for by its new family, which holds no genus of the other side
    assert by["Onchidella;Siphonaria"]["source"] == "timetree_higher_taxa"
    assert any(7000 in q for q in queries) and not any(SAND_DOLLAR & q for q in queries)
    # no reference at all: nothing can be checked, nothing is trusted
    blind = ra.IdScreen(Classify(CLASSIFICATIONS), {"Onchidella": (None, "")}, cc.LEVELS["class"])
    assert not blind.keep("Onchidella", 1001) and blind.decisions["Onchidella"][1001] == "no_reference"
    assert blind.audit_row("Onchidella")["reanchor_ids_against"] == "none"


@pytest.mark.parametrize("basis", ["global", "sampled"])
def test_crown_queries_never_carry_the_other_animal(taxonomy, basis):
    screen = screen_for(taxonomy, Classify(CLASSIFICATIONS))
    # sampled: names from the records, one of them a sand dollar filed under the genus
    sampled = {"Onchidella": {"Onchidella celtica", "Onchidella floridana", "Peronella japonica"}}
    for keep, bad in ((None, True), (screen.keep, False)):
        queries = []
        rows = crown_ages(taxonomy, recording_timetree(queries), ["Onchidella"], sampled=sampled, basis=basis,
                          keep=keep)
        assert any(SAND_DOLLAR & q for q in queries) is bad
    row = rows[0]
    assert set(map(int, row["query_ids"].split("+"))) == ({1001, 1002} if basis == "sampled" else {1001, 1002, 1003})
    assert row["crown_age"] == "50.0000"
    if basis == "sampled":
        assert row["n_sampled_species"] == 2 and row["crown_basis"] == "sampled"
