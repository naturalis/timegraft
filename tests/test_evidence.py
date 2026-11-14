"""Evidence and provenance kept for the consistency layer (#47): none of it may change an age."""

import json
from pathlib import Path

import dendropy
import pytest

from timegraft.backbone_build import (GENERA_COLUMNS, REPORT_COLUMNS, date_topology, genera_rows, records_lineage,
                                      timetree_summary)
from timegraft.otol_static import (CROWN_COLUMNS, SUPPORT_KINDS, SynthTree, Taxonomy, crown_ages, genus_anchor_info,
                                   genus_anchors, genus_topology_newick, lineage_names, synth_nodes, synthesis_support)
from timegraft.timetree import EVIDENCE_COLUMNS, TimeTree, is_legacy, parse_response
from timegraft.webapi import JsonCache

DATA = Path(__file__).parent / "data" / "otol"

# shaped like an answer of api.timetree.org/mrca/id/9606+10090/json (shortened)
RESPONSE = {
    "found_ids": [9606, 10090],
    "missing_ids": [],
    "mrca_ttid": 266386,
    "taxon_id": 314146,
    "time_estimates": "{86,91,94.5,70.14}",
    "precomputed_age": 83.51395605057974,
    "precomputed_ci_low": 80.13669528675152,
    "precomputed_ci_high": 86.89121681440795,
    "adjusted_age": 83.51396,
    "species_data": {"9606": {"scientific_name": "Homo sapiens"}},
    "study_data": [
        {"first_author": "Douady", "title": "A long title", "pub_year": 2003, "tree_count": 1, "age_estimate": 86},
        {"first_author": "Poux", "title": "Another", "pub_year": 2006, "tree_count": 2, "age_estimate": 91},
    ],
    "version": 20260820,
}


def test_parse_response_keeps_the_evidence():
    entry = parse_response(RESPONSE)
    assert entry["v"] == 2
    assert entry["found"] == [9606, 10090] and entry["age"] == pytest.approx(83.514, abs=1e-3)
    assert entry["mrca_ncbi"] == 314146
    assert entry["estimates"] == [86.0, 91.0, 94.5, 70.14]
    assert entry["ci_low"] == pytest.approx(80.137, abs=1e-3)
    assert entry["studies"][1] == {"author": "Poux", "year": 2006, "trees": 2, "age": 91.0}
    assert "title" not in entry["studies"][0]  # titles are not needed and would bloat the cache
    assert entry["api_version"] == "20260820"


def test_parse_response_same_age_rules_as_before():
    # a single estimate: no interval; a zero age means no age, as in cache version 1
    single = parse_response({"found_ids": [1, 2], "precomputed_age": 1.8, "precomputed_ci_low": 0,
                             "precomputed_ci_high": 0, "time_estimates": "{1.8}"})
    assert single["age"] == 1.8 and single["ci_low"] is None and single["ci_high"] is None
    assert parse_response({"found_ids": [1], "precomputed_age": 0})["age"] is None
    nothing = parse_response({"found_ids": [1], "missing_ids": [2], "error": "No TimeTree study info"})
    assert nothing["age"] is None and nothing["missing"] == [2] and nothing["error"].startswith("No TimeTree")
    assert parse_response(None) == {"v": 2, "found": [], "age": None, "error": "no answer"}


def test_legacy_cache_entries_are_read_unchanged():
    cache = JsonCache(None)
    cache.set("mrca:1+2", {"found": [1, 2], "age": 50.0})

    def fetch(url, payload=None):
        raise AssertionError("a cached query must not be asked again without upgrade_legacy")

    tt = TimeTree(cache, fetch=fetch)
    assert tt.split_age([[1], [2]]) == 50.0
    assert is_legacy(cache.get("mrca:1+2"))
    evidence = tt.evidence([1, 2])
    assert evidence["evidence"] == "legacy_cache" and evidence["found_ids"] == "1+2"
    assert evidence["n_estimates"] == ""
    assert tt.evidence([3, 4])["evidence"] == "none"


@pytest.mark.parametrize("prefetch", [False, True])
def test_upgrade_keeps_the_cached_age(prefetch):
    cache = JsonCache(None)
    cache.set("mrca:1+2", {"found": [1, 2], "age": 50.0})
    cache.set("mrca:3+4", {"found": [3, 4], "age": 20.0})

    def fetch(url, payload=None):
        ids = [int(i) for i in url.split("/mrca/id/")[1].split("/")[0].split("+")]
        age = 55.0 if ids == [1, 2] else 20.0  # TimeTree has moved on for 1+2
        return {"found_ids": ids, "precomputed_age": age, "time_estimates": "{19,21}", "taxon_id": 99,
                "study_data": [{"first_author": "X", "age_estimate": 19}, {"first_author": "Y", "age_estimate": 21}]}

    tt = TimeTree(cache, fetch=fetch, upgrade_legacy=True)
    if prefetch:
        tt.prefetch([[1, 2], [3, 4]])
    assert tt.split_age([[1], [2]]) == 50.0
    assert tt.split_age([[3], [4]]) == 20.0
    moved, same = cache.get("mrca:1+2"), cache.get("mrca:3+4")
    assert moved["v"] == 2 and moved["upgraded_from"] == 1
    assert moved["age"] == 50.0 and moved["fresh_age"] == 55.0
    assert "fresh_age" not in same and "fresh_found" not in same
    evidence = tt.evidence([3, 4])
    assert evidence["evidence"] == "full" and evidence["n_studies"] == "2" and evidence["mrca_ncbi"] == "99"
    assert evidence["estimate_min"] == "19.0000" and evidence["estimate_max"] == "21.0000"


def test_upgrade_survives_network_trouble():
    cache = JsonCache(None)
    cache.set("mrca:1+2", {"found": [1, 2], "age": 50.0})

    def fetch(url, payload=None):
        raise OSError("tunnel dropped")

    tt = TimeTree(cache, fetch=fetch, upgrade_legacy=True)
    assert tt.split_age([[1], [2]]) == 50.0
    assert is_legacy(cache.get("mrca:1+2"))


def test_upgrade_keeps_the_legacy_entry_on_an_http_error():
    from timegraft.webapi import ApiError

    cache = JsonCache(None)
    cache.set("mrca:1+2", {"found": [1, 2], "age": 50.0})

    def fetch(url, payload=None):
        raise ApiError(503, "service unavailable")

    tt = TimeTree(cache, fetch=fetch, upgrade_legacy=True)
    assert tt.split_age([[1], [2]]) == 50.0
    assert is_legacy(cache.get("mrca:1+2")) and "error" not in cache.get("mrca:1+2")
    tt.prefetch([[1, 2]])  # the concurrent path too
    assert is_legacy(cache.get("mrca:1+2")) and "error" not in cache.get("mrca:1+2")


def _fake(ages, missing=()):
    def fetch(url, payload=None):
        ids = [int(i) for i in url.split("/mrca/id/")[1].split("/")[0].split("+")]
        found = [i for i in ids if i not in missing]
        age = ages.get(frozenset(found))
        return {"found_ids": found, "precomputed_age": age, "time_estimates": f"{{{age}}}" if age else None,
                "taxon_id": 1000 + len(found), "version": 1}

    return fetch


def test_report_carries_evidence_without_changing_ages():
    newick = "(((A,B),C),(D,E));"
    ncbi = {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5}
    ages = {frozenset({1, 2, 4, 5}): 100.0, frozenset({1, 2, 3}): 120.0, frozenset({1, 2}): 30.0}
    plain = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted")
    _, before = date_topology(plain, ncbi, TimeTree(fetch=_fake(ages)), method="precedence")
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted")
    labels = {n: f"mrcaott{i}" for i, n in enumerate(tree.preorder_node_iter()) if not n.is_leaf()}
    support = {labels[tree.seed_node]: {"supported_by": 3, "conflicts_with": 1}}
    _, report = date_topology(tree, ncbi, TimeTree(fetch=_fake(ages)), synth_nodes=labels, support=support,
                                method="precedence")
    assert [(r["age"], r["source"]) for r in report] == [(r["age"], r["source"]) for r in before]
    assert all(set(r) == set(REPORT_COLUMNS) for r in report)
    by = {r["example_genera"]: r for r in report}
    abc, ab = by["A;B;C"], by["A;B"]
    root = next(r for r in report if r["parent"] == "")
    assert root["timetree_evidence"] == "full" and root["timetree_n_estimates"] == "1"
    assert root["synth_supported_by"] == "3" and root["synth_conflicts_with"] == "1"
    assert ab["synth_supported_by"] == "0"  # no annotations: taxonomy only
    # (A,B,C) at 120 is older than the root at 100: dropped, and the report says against which node
    assert abc["source"] == "interpolated_after_drop" and abc["conflict_with"] == root["node"]
    assert abc["timetree_attempt"] == "0" and abc["timetree_evidence"] == "full"
    assert by["D;E"]["timetree_evidence"] == "" and by["D;E"]["timetree_attempt"] == ""
    assert all(r["fitted_age"] == "" and r["conflict_class"] == "" for r in report)


def test_report_without_synthesis_support_leaves_it_empty():
    tree = dendropy.Tree.get(data="(A,B);", schema="newick", rooting="force-rooted")
    _, report = date_topology(tree, {"A": 1, "B": 2}, TimeTree(fetch=_fake({frozenset({1, 2}): 5.0})))
    assert report[0]["synth_node"] == "" and report[0]["synth_supported_by"] == ""


@pytest.fixture()
def merged_taxonomy(tmp_path):
    # Mytilus with two NCBI ids, as OTT does for a wrong merge (cf. Onchidella in #20)
    text = (DATA / "taxonomy.tsv").read_text().replace(
        "5\t|\t4\t|\tMytilus\t|\tgenus\t|\tncbi:6548\t|",
        "5\t|\t4\t|\tMytilus\t|\tgenus\t|\tncbi:6548,worms:138228,ncbi:262330\t|",
    )
    path = tmp_path / "taxonomy.tsv"
    path.write_text(text)
    return path


def test_match_info_and_lineage(merged_taxonomy):
    genera = ["Mytilus", "Ostrea", "Magallana", "Asterias", "Nowhere"]
    taxonomy = Taxonomy(merged_taxonomy, set(genera), context="Metazoa")
    assert taxonomy.genus["Mytilus"] == {"ott_id": 5, "ncbi_id": 6548, "status": "matched"}  # match unchanged
    info = taxonomy.match_info["Mytilus"]
    assert info["ott_candidates"] == 2 and info["ott_candidates_in_context"] == 1  # the fungal homonym
    assert info["ott_ncbi_ids"] == "6548;262330" and info["ott_worms_ids"] == "138228"
    assert taxonomy.match_info["Nowhere"]["ott_candidates"] == 0
    assert lineage_names(taxonomy, "Mytilus") == {"ott_family": "Mytilidae", "ott_order": "", "ott_class": "Bivalvia",
                                                  "ott_phylum": ""}
    ambiguous = Taxonomy(merged_taxonomy, {"Mytilus"})
    assert ambiguous.match_info["Mytilus"]["ott_candidate_ids"] == "5;17"


def test_anchor_info_and_synthesis_nodes(tmp_path):
    genera = ["Mytilus", "Ostrea", "Magallana", "Asterias"]
    taxonomy = Taxonomy(DATA / "taxonomy.tsv", set(genera), context="Metazoa")
    synth = SynthTree(DATA / "labelled_supertree.tre")
    anchors, _ = genus_anchors(synth, taxonomy, genera)
    info = genus_anchor_info(synth, taxonomy, anchors)
    assert info["Mytilus"] == {"anchor_node": "ott5", "anchor_clade_genera": "1", "n_species_ott": "2",
                               "n_species_in_synthesis": "2"}
    assert info["Magallana"]["anchor_node"] == "mrcaott6ott12"  # attached above its single species
    newick, _ = genus_topology_newick(synth, taxonomy, genera)
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True)
    labels = {n: synth.label[m] for n, m in synth_nodes(synth, anchors, tree).items() if not n.is_leaf()}
    assert labels[tree.seed_node] == "ott2"

    path = tmp_path / "annotations.json"
    path.write_text(json.dumps({"nodes": {
        "ott2": {"supported_by": {"a@t1": "n1", "b@t2": "n2"}, "conflicts_with": {"c@t3": ["n9"]}},
        "mrcaott6ott12": {"resolves": {"a@t1": "n3"}},
        "ott99": None,
    }}))
    support = synthesis_support(path, {"ott2", "mrcaott6ott12", "ott99"})
    assert support == {"ott2": {"supported_by": 2, "conflicts_with": 1}, "mrcaott6ott12": {"resolves": 1}}


@pytest.mark.parametrize("chunk", [1, 2, 3, 7, 1 << 20])
def test_synthesis_support_streams_like_json_load(tmp_path, monkeypatch, chunk):
    # #72: read incrementally; same result as json.load at any chunk boundary
    from timegraft import otol_static
    monkeypatch.setattr(otol_static._JSONStream, "CHUNK", chunk)
    document = {
        "date_completed": "2025-12-20", "num_tips": 12345678, "sources": ["ot_1@t1", "ot_2@t2"],
        "source_id_map": {"ot_1@t1": {"study_id": "ot_1", "nested": [1, 2.5e-3, True, None, "\\u00e9 \\" ]}},
        "nodes": {
            "ott2": {"supported_by": {"a@t1": "n1", "b@t2": "n2"}, "conflicts_with": {"c@t3": ["n9"]}},
            "mrcaott6ott12": {"resolves": {"a@t1": "n3"}, "was_uncontested": True},
            "ott99": None,
            "ott7": {},
            "ott1234567": {"terminal": {"x@t9": "ott1234567"}, "supported_by": []},
        },
        "tree_id": "opentree16.1", "root_ott_id": 93302,
    }
    path = tmp_path / "annotations.json"
    for text in (json.dumps(document), json.dumps(document, indent=1), json.dumps(document, separators=(",", ":"))):
        path.write_text(text, encoding="utf-8")
        expected = {k: {kind: len(v[kind]) for kind in SUPPORT_KINDS if isinstance(v.get(kind), dict)}
                    for k, v in json.loads(text)["nodes"].items() if isinstance(v, dict)}
        assert synthesis_support(path) == expected
        assert list(synthesis_support(path)) == list(expected)
        assert synthesis_support(path, {"ott2", "ott1234567", "ott404"}) == {
            "ott2": {"supported_by": 2, "conflicts_with": 1}, "ott1234567": {"terminal": 1}}


def test_synthesis_support_edge_cases(tmp_path):
    path = tmp_path / "annotations.json"
    path.write_text('{"tree_id": "x"}')
    assert synthesis_support(path) == {}  # no nodes
    path.write_text('{"nodes": {}, "tree_id": "x"}')
    assert synthesis_support(path) == {}
    path.write_text('{"nodes": null}')
    assert synthesis_support(path) == {}
    # the last of duplicate keys wins, as with json.load
    path.write_text('{"nodes": {"ott1": {"resolves": {"a": 1}}, "ott1": {"supported_by": {"b": 1, "c": 2}}}}')
    assert synthesis_support(path) == {"ott1": {"supported_by": 2}}
    path.write_text('{"nodes": {"ott1": {"resolves": {"a": 1}}, "ott1": null}}')
    assert synthesis_support(path) == {}
    for broken in ('{"nodes": {"ott1": {"resolves": {"a": 1}}', '{"nodes": {"ott1" {}}}', '["nodes"]', ''):
        path.write_text(broken)
        with pytest.raises(ValueError):
            synthesis_support(path)


def test_crown_ages_carry_evidence():
    taxonomy = Taxonomy(DATA / "taxonomy.tsv", {"Mytilus"}, context="Metazoa")
    tt = TimeTree(fetch=_fake({frozenset({6550, 6551}): 4.5}))
    rows = crown_ages(taxonomy, tt, ["Mytilus"])
    assert len(rows) == 1 and set(rows[0]) == set(CROWN_COLUMNS)
    assert rows[0]["crown_age"] == "4.5000" and rows[0]["query_ids"] == "6550+6551"
    assert rows[0]["timetree_evidence"] == "full" and rows[0]["timetree_mrca_ncbi"] == "1002"
    assert set(EVIDENCE_COLUMNS) <= {c[len("timetree_"):] for c in CROWN_COLUMNS if c.startswith("timetree_")}


def test_genera_rows_and_records_lineage():
    records = [
        {"genus": "Onchidella", "family": "Onchidiidae", "class": "Gastropoda", "phylum": "Mollusca"},
        {"genus": "Onchidella", "family": "Onchidiidae", "class": "Gastropoda", "phylum": "Mollusca"},
        {"genus": "Onchidella", "family": "Laganidae", "class": "Echinoidea", "phylum": "Echinodermata"},
        {"genus": "", "family": "X"},
    ]
    lineage = records_lineage(records)
    assert lineage["Onchidella"]["records_n"] == "3"
    assert lineage["Onchidella"]["records_family"] == "Onchidiidae"
    assert lineage["Onchidella"]["records_family_share"] == "0.667"
    assert lineage["Onchidella"]["records_order"] == "" and lineage["Onchidella"]["records_order_share"] == ""
    matches = {"Onchidella": {"status": "matched", "ott_id": 1, "ncbi_id": None}, "Other": {"status": "unmatched"}}
    rows = genera_rows(matches, lineage, {"Onchidella": {"ott_ncbi_ids": "36932;262330", "not_a_column": "x"}})
    assert [r["genus"] for r in rows] == ["Onchidella", "Other"]
    assert all(list(r) == list(GENERA_COLUMNS) for r in rows)
    assert rows[0]["ncbi_id"] == "" and rows[0]["records_phylum"] == "Mollusca"
    assert rows[0]["ott_ncbi_ids"] == "36932;262330" and rows[1]["records_n"] == ""


def test_timetree_summary():
    report = [{"query_ids": "1+2", "timetree_evidence": "full", "timetree_api_version": "20260820"},
              {"query_ids": "3+4", "timetree_evidence": "legacy_cache", "timetree_api_version": ""},
              {"query_ids": "", "timetree_evidence": ""}]
    assert timetree_summary(report) == {"queries_by_evidence": {"full": 1, "legacy_cache": 1},
                                        "api_versions": {"20260820": 1}}
