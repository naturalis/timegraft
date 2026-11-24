"""Which node did TimeTree date? (#108, step 2)"""

import dendropy
import pytest

from timegraft import ttcheck
from timegraft.backbone_build import date_topology
from timegraft.timetree import TimeTree

# NCBI: 1 > 10 > 30 > A (101), B (102), C (103); 10 > D (104). The backbone is (((A,B),C),D), so NCBI's taxon
# 30 is the node of A, B and C, and the node of A and B has no NCBI taxon of its own.
LINEAGE = {"101": ["1", "10", "30", "101"], "102": ["1", "10", "30", "102"], "103": ["1", "10", "30", "103"],
           "104": ["1", "10", "104"], "30": ["1", "10", "30"], "10": ["1", "10"], "1": ["1"]}
IDS = {"A": 101, "B": 102, "C": 103, "D": 104}


def lineage(taxon):
    return LINEAGE.get(str(taxon))


def test_common_ancestor():
    assert ttcheck.common_ancestor([LINEAGE["101"], LINEAGE["102"]]) == "30"
    assert ttcheck.common_ancestor([LINEAGE["101"], LINEAGE["104"]]) == "10"
    assert ttcheck.common_ancestor([LINEAGE["101"]]) == "101"
    assert ttcheck.common_ancestor([["1", "2"], ["3"]]) == ""


def test_check_outcomes():
    known = ttcheck.genus_lineages(IDS, lineage)
    assert ttcheck.check([101, 103], {"A", "B", "C"}, known, lineage)["outcome"] == "match"
    wider = ttcheck.check([101, 102], {"A", "B"}, known, lineage)
    assert wider["outcome"] == "wider" and wider["ancestor"] == "30" and wider["outside"] == {"C"}
    narrower = ttcheck.check([101, 102], {"A", "B", "C", "D"}, known, lineage)
    assert narrower["outcome"] == "narrower" and narrower["missing"] == {"D"}
    assert ttcheck.check([101, 999], {"A", "B"}, known, lineage)["outcome"] == "unchecked"
    assert ttcheck.check([], {"A"}, known, lineage)["outcome"] == "unchecked"
    # genera without a known lineage count on neither side
    assert ttcheck.check([101, 103], {"A", "B", "C", "X"}, known, lineage)["outcome"] == "match"


def _query(answers):
    def query(ids):
        return answers.get(tuple(sorted(ids)), {"found": [], "mrca_ttid": None})
    return query


def test_probe():
    ids = [101, 102]
    inside = _query({(101, 102, 103): {"found": [101, 102, 103], "mrca_ttid": 7}})
    assert ttcheck.probe(ids, 7, ["C"], IDS, inside) == ("wider", ["C"])
    outside = _query({(101, 102, 103): {"found": [101, 102, 103], "mrca_ttid": 8}})
    assert ttcheck.probe(ids, 7, ["C"], IDS, outside) == ("probed_outside", ["C"])
    # a probe TimeTree does not know is skipped; with none left, the answer stays unchecked
    unknown = _query({(101, 102, 103): {"found": [101, 102], "mrca_ttid": 7}})
    assert ttcheck.probe(ids, 7, ["C"], IDS, unknown) == ("unchecked", [])
    both = _query({(101, 102, 103): {"found": [101, 102], "mrca_ttid": 7},
                   (101, 102, 104): {"found": [101, 102, 104], "mrca_ttid": 7}})
    assert ttcheck.probe(ids, 7, ["C", "D"], IDS, both) == ("wider", ["D"])
    # without TimeTree's node id nothing can be probed
    assert ttcheck.probe(ids, None, ["C"], IDS, inside) == ("unchecked", [])
    # at most `probes` known genera are tried
    many = _query({(101, 102, 103): {"found": [101, 102, 103], "mrca_ttid": 8},
                   (101, 102, 104): {"found": [101, 102, 104], "mrca_ttid": 7}})
    assert ttcheck.probe(ids, 7, ["C", "D"], IDS, many, probes=1) == ("probed_outside", ["C"])
    assert ttcheck.probe(ids, 7, ["C", "D"], IDS, many, probes=2) == ("wider", ["C", "D"])


def test_nearest_outside():
    tree = dendropy.Tree.get(data="((((A,B),C),D),E);", schema="newick", rooting="force-rooted")
    ab = tree.mrca(taxon_labels=["A", "B"])
    assert ttcheck.nearest_outside(ab, {"C", "D", "E"}) == ["C", "D", "E"]
    assert ttcheck.nearest_outside(ab, {"E", "C"}) == ["C", "E"]
    assert ttcheck.nearest_outside(ab, {"A"}) == []


def _fetch(named_ab: bool, ab_ttid: int, probe_ttid: int):
    """TimeTree's answers for (((A,B),C),D): root 100, ABC 50, AB 40. The probe of (A,B) with C added is the
    query A+B+C, whose node id is ``probe_ttid``."""

    def fetch(url, payload=None):
        ids = {int(i) for i in url.split("/mrca/id/")[1].split("/")[0].split("+")}
        if 104 in ids:
            age, taxon, ttid = 100.0, 10, 1
        elif 103 in ids:
            age, taxon, ttid = 50.0, 30, probe_ttid if ids == {101, 102, 103} else 3
        else:
            age, taxon, ttid = 40.0, 30 if named_ab else 0, ab_ttid
        return {"found_ids": sorted(ids), "precomputed_age": age, "taxon_id": taxon, "mrca_ttid": ttid,
                "time_estimates": "{" + str(age) + "}"}

    return fetch


def _date(mode, named_ab=True, ab_ttid=5, probe_ttid=3):
    tree = dendropy.Tree.get(data="(((A,B),C),D);", schema="newick", rooting="force-rooted")
    _, report = date_topology(tree, dict(IDS), TimeTree(fetch=_fetch(named_ab, ab_ttid, probe_ttid)),
                              per_side=2, node_check={"mode": mode, "probes": 3}, lineage=lineage)
    by = {r["example_genera"]: r for r in report if int(r["n_genera"]) < 4}
    by["root"] = report[0]
    return by


def test_named_wider_answer_is_rejected():
    # TimeTree's answer for (A,B) is NCBI taxon 30, which also holds C: it dated the node of A, B and C
    by = _date("reject")
    ab = by["A;B"]
    assert ab["timetree_node_check"] == "wider" and ab["timetree_ncbi_ancestor"] == "30"
    assert ab["timetree_genera_outside"] == "1" and ab["timetree_genera_missing"] == "0"
    assert ab["source"] == "interpolated_timetree_other_node" and ab["calibration_sigma"] == ""
    assert ab["timetree_age"] == "40.0000" and 0 < float(ab["age"]) < 50
    assert by["A;B;C"]["timetree_node_check"] == "match" and float(by["A;B;C"]["age"]) == pytest.approx(50.0)
    assert by["root"]["timetree_node_check"] == "match"


def test_report_mode_keeps_the_answer():
    ab = _date("report")["A;B"]
    assert ab["timetree_node_check"] == "wider" and float(ab["age"]) == pytest.approx(40.0)
    assert ab["source"].startswith("timetree")


def test_off_records_nothing():
    ab = _date("off")["A;B"]
    assert ab["timetree_node_check"] == "" and float(ab["age"]) == pytest.approx(40.0)


def test_unnamed_answer_is_probed():
    # unnamed (TimeTree made the node when resolving NCBI's taxon 30): adding C gives TimeTree's node 3, not 5,
    # so C is outside the node TimeTree dated, and the answer is kept
    outside = _date("reject", named_ab=False, ab_ttid=5, probe_ttid=3)["A;B"]
    assert outside["timetree_node_check"] == "probed_outside" and outside["timetree_probes"] == "C"
    assert float(outside["age"]) == pytest.approx(40.0)
    # adding C gives the same node: TimeTree's node holds C, and the answer is not used
    inside = _date("reject", named_ab=False, ab_ttid=3, probe_ttid=3)["A;B"]
    assert inside["timetree_node_check"] == "wider" and inside["source"] == "interpolated_timetree_other_node"


def test_root_keeps_its_answer():
    # (((A,B),C),(D,G)), with G (107) outside NCBI's taxon 10 and unknown to TimeTree: the root's answer
    # spans A to D only, so it is narrower than the root. It is recorded but kept: nothing else dates the root
    tree = dendropy.Tree.get(data="(((A,B),C),(D,G));", schema="newick", rooting="force-rooted")
    lineages = {**LINEAGE, "107": ["1", "107"]}

    def fetch(url, payload=None):
        ids = [int(i) for i in url.split("/mrca/id/")[1].split("/")[0].split("+")]
        return {"found_ids": [i for i in ids if i != 107], "precomputed_age": 100.0, "taxon_id": 10,
                "mrca_ttid": 1, "time_estimates": "{100.0}"}

    _, report = date_topology(tree, {**IDS, "G": 107}, TimeTree(fetch=fetch), node_check={"mode": "reject"},
                              lineage=lambda t: lineages.get(str(t)))
    root = report[0]
    assert root["timetree_node_check"] == "narrower" and root["timetree_genera_missing"] == "1"
    assert float(root["age"]) == pytest.approx(100.0)


def test_unknown_mode():
    tree = dendropy.Tree.get(data="(A,B);", schema="newick", rooting="force-rooted")
    with pytest.raises(ValueError):
        date_topology(tree, {"A": 101, "B": 102}, TimeTree(fetch=_fetch(True, 5, 3)), node_check={"mode": "drop"},
                      lineage=lineage)


def test_ncbi_lineage_from_checklistbank():
    from timegraft.col import ChecklistBank

    calls = []

    def fetch(url, payload=None):
        calls.append(url)
        if url.endswith("/taxon/7608/classification"):
            return [{"id": "1", "name": "root"}, {"id": "33208", "name": "Metazoa"}, {"id": "7600", "name": "Asteriidae"}]
        return []

    client = ChecklistBank(316691, fetch=fetch)
    assert client.ncbi_lineage(7608) == ["1", "33208", "7600", "7608"]
    assert client.ncbi_lineage(7608) == ["1", "33208", "7600", "7608"] and len(calls) == 1  # cached
    assert client.ncbi_lineage(999) is None
    assert "/dataset/2169/taxon/7608/classification" in calls[0]
