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
    # both genera outside and genera missing: a focal genus is left out, so it is narrower whatever is outside
    both = ttcheck.check([101, 102], {"A", "B", "D"}, known, lineage)
    assert both["outcome"] == "narrower" and both["outside"] == {"C"} and both["missing"] == {"D"}


def test_named_answers_are_checked_against_the_reported_taxon():
    # today's NCBI has a clade 31 of A and B inside 30, which TimeTree's NCBI release lacked: TimeTree reports 30
    # for A+B. Checked against 30 (what TimeTree dated) the answer is wider; against 31 it would pass
    newer = {**LINEAGE, "101": ["1", "10", "30", "31", "101"], "102": ["1", "10", "30", "31", "102"],
             "31": ["1", "10", "30", "31"]}
    look = newer.get
    known = ttcheck.genus_lineages(IDS, lambda t: newer.get(str(t)))
    assert ttcheck.check([101, 102], {"A", "B"}, known, lambda t: look(str(t)))["outcome"] == "match"
    named = ttcheck.check([101, 102], {"A", "B"}, known, lambda t: look(str(t)), reported=30)
    assert named["outcome"] == "wider" and named["ancestor"] == "30" and named["outside"] == {"C"}
    # a reported taxon that does not hold the recognised ids in today's NCBI cannot be checked safely
    assert ttcheck.check([101, 104], {"A", "D"}, known, lambda t: look(str(t)), reported=30)["outcome"] == "unchecked"
    assert ttcheck.check([101, 102], {"A", "B"}, known, lambda t: look(str(t)), reported=77)["outcome"] == "unchecked"


def _query(answers):
    def query(ids):
        return answers.get(tuple(sorted(ids)), {"found": [], "mrca_ttid": None})
    return query


def test_probe():
    ids = [101, 102]
    inside = _query({(101, 102, 103): {"found": [101, 102, 103], "mrca_ttid": 7}})
    assert ttcheck.probe(ids, 7, ["C"], IDS, inside) == ("wider", [("C", True)])
    outside = _query({(101, 102, 103): {"found": [101, 102, 103], "mrca_ttid": 8}})
    assert ttcheck.probe(ids, 7, ["C"], IDS, outside) == ("passed", [("C", False)])
    # a probe TimeTree does not know is skipped; with none left, the answer stays unchecked
    unknown = _query({(101, 102, 103): {"found": [101, 102], "mrca_ttid": 7}})
    assert ttcheck.probe(ids, 7, ["C"], IDS, unknown) == ("unchecked", [])
    # so is an answer without a node id: it neither passes nor uses up the limit
    no_id = _query({(101, 102, 103): {"found": [101, 102, 103], "mrca_ttid": None},
                    (101, 102, 104): {"found": [101, 102, 104], "mrca_ttid": 7}})
    assert ttcheck.probe(ids, 7, ["C", "D"], IDS, no_id, probes=1) == ("wider", [("D", True)])
    # without TimeTree's node id for the answer itself nothing can be probed
    assert ttcheck.probe(ids, None, ["C"], IDS, inside) == ("unchecked", [])
    # at most `probes` known genera are tried
    many = _query({(101, 102, 103): {"found": [101, 102, 103], "mrca_ttid": 8},
                   (101, 102, 104): {"found": [101, 102, 104], "mrca_ttid": 7}})
    assert ttcheck.probe(ids, 7, ["C", "D"], IDS, many, probes=1) == ("passed", [("C", False)])
    assert ttcheck.probe(ids, 7, ["C", "D"], IDS, many, probes=2) == ("wider", [("C", False), ("D", True)])
    # probing for coverage: a genus of the node outside TimeTree's node makes it narrower
    assert ttcheck.probe(ids, 7, ["C"], IDS, outside, stop_when_inside=False) == ("narrower", [("C", False)])
    assert ttcheck.probe(ids, 7, ["C"], IDS, inside, stop_when_inside=False) == ("passed", [("C", True)])


def test_represented_and_interleave():
    assert ttcheck.represented([101, 102], {"A", "B", "C"}, IDS, lineage) == {"A", "B"}
    # a species id stands for its genus
    species = {**LINEAGE, "1011": ["1", "10", "30", "101", "1011"]}
    assert ttcheck.represented([1011], {"A", "B"}, IDS, lambda t: species.get(str(t))) == {"A"}
    tree = dendropy.Tree.get(data="((A,B,C),(D,E));", schema="newick", rooting="force-rooted")
    assert ttcheck.interleave(tree.seed_node) == ["A", "D", "B", "E", "C"]


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
    assert outside["timetree_node_check"] == "probed_outside" and outside["timetree_probes"] == "C:out"
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


def _unnamed_fetch(clade: set[int], clade_ttid: int, other_ttid: int, unknown=()):
    """An unnamed TimeTree clade: any query within ``clade`` gets node ``clade_ttid``, any other ``other_ttid``."""

    def fetch(url, payload=None):
        q = {int(i) for i in url.split("/mrca/id/")[1].split("/")[0].split("+")}
        found = sorted(q - set(unknown))
        inside = set(found) <= clade
        return {"found_ids": found, "precomputed_age": 30.0 if inside else 100.0, "taxon_id": 0,
                "mrca_ttid": clade_ttid if inside else other_ttid, "time_estimates": "{30.0}" if inside else "{100.0}"}

    return fetch


def test_unnamed_match_that_leaves_out_a_focal_genus_is_rejected():
    # NCBI's taxon 30 holds A, B, D and Z; TimeTree resolves it with a clade of A, B and D that leaves Z out.
    # Asked about OToL's node ((A,B),(D,Z)) with A and D, TimeTree dates that clade: the NCBI comparison matches,
    # B probes inside, but Z probes outside, so the answer dates a narrower node and is rejected
    lineages = {"101": ["1", "30", "101"], "102": ["1", "30", "102"], "104": ["1", "30", "104"],
                "126": ["1", "30", "126"], "105": ["1", "105"], "30": ["1", "30"], "1": ["1"]}
    ids = {"A": 101, "B": 102, "D": 104, "Z": 126, "E": 105}
    tree = dendropy.Tree.get(data="(((A,B),(D,Z)),E);", schema="newick", rooting="force-rooted")
    _, report = date_topology(tree, ids, TimeTree(fetch=_unnamed_fetch({101, 102, 104}, 9, 3)), per_side=1,
                              per_genus=0, node_check={"mode": "reject"}, lineage=lambda t: lineages.get(str(t)))
    node = next(r for r in report if int(r["n_genera"]) == 4)
    assert set(node["query_genera"].split(";")) <= {"A", "B", "D"}  # Z not among the picks
    assert node["timetree_node_check"] == "narrower" and node["timetree_probes"].endswith("Z:out")
    assert node["source"] == "interpolated_timetree_other_node"
    # with a TimeTree clade that holds Z as well, the same answer is kept
    _, report = date_topology(tree, ids, TimeTree(fetch=_unnamed_fetch({101, 102, 104, 126}, 9, 3)), per_side=1,
                              per_genus=0, node_check={"mode": "reject"}, lineage=lambda t: lineages.get(str(t)))
    node = next(r for r in report if int(r["n_genera"]) == 4)
    assert node["timetree_node_check"] == "match" and "Z:in" in node["timetree_probes"]
    assert node["source"].startswith("timetree")


def test_outside_and_missing_is_rejected_without_probing():
    # NCBI's taxon 30 holds A, B, C and D, not F. The backbone node ((A,B),(D,F)) is asked about with A and D:
    # C is outside the node but inside the ancestor, and F is in the node but outside it. F is left out whatever
    # the probe of C would say, so the answer is rejected and nothing is probed
    lineages = {**LINEAGE, "104": ["1", "10", "30", "104"], "106": ["1", "106"], "105": ["1", "105"]}
    ids = {**IDS, "F": 106, "E": 105}
    tree = dendropy.Tree.get(data="((((A,B),(D,F)),C),E);", schema="newick", rooting="force-rooted")
    calls = []
    fetch = _unnamed_fetch({101, 102, 104}, 9, 3)

    def counting(url, payload=None):
        calls.append(url)
        return fetch(url, payload)

    _, report = date_topology(tree, ids, TimeTree(fetch=counting), per_side=1, per_genus=0,
                              node_check={"mode": "reject"}, lineage=lambda t: lineages.get(str(t)))
    node = next(r for r in report if int(r["n_genera"]) == 4)
    assert set(node["query_genera"].split(";")) == {"A", "D"}
    assert node["timetree_node_check"] == "narrower" and node["timetree_genera_outside"] == "1"
    assert node["timetree_genera_missing"] == "1" and node["timetree_probes"] == ""
    assert node["source"] == "interpolated_timetree_other_node"


def test_yaml_off_and_the_adjusted_age_of_a_rejected_answer():
    # an unquoted `off` in YAML arrives as False: the check is off, not an error
    tree = dendropy.Tree.get(data="(((A,B),C),D);", schema="newick", rooting="force-rooted")
    _, report = date_topology(tree, dict(IDS), TimeTree(fetch=_fetch(True, 5, 3)), per_side=2,
                              node_check={"mode": False}, lineage=lineage)
    assert all(r["timetree_node_check"] == "" for r in report)

    # a rejected answer keeps its adjusted time in the report, and stays out of the calibrations
    def fetch(url, payload=None):
        answer = _fetch(True, 5, 3)(url, payload)
        if answer["precomputed_age"] == 40.0:
            answer["adjusted_age"] = 45.0
        return answer

    tree = dendropy.Tree.get(data="(((A,B),C),D);", schema="newick", rooting="force-rooted")
    _, report = date_topology(tree, dict(IDS), TimeTree(fetch=fetch), per_side=2, node_check={"mode": "reject"},
                              lineage=lineage, timetree_time={"use": "adjusted"})
    ab = next(r for r in report if r["example_genera"] == "A;B")
    assert ab["timetree_node_check"] == "wider" and ab["timetree_adjusted_age"] == "45.0000"
    assert ab["calibration_age"] == "" and ab["source"] == "interpolated_timetree_other_node"
    assert float(ab["age"]) != pytest.approx(45.0)
