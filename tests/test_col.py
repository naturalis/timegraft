"""COL XR as identifier hub (#48), with ChecklistBank mocked."""

from urllib.parse import parse_qs, urlparse

import pytest

from timegraft.backbone_build import join_ncbi
from timegraft.col import (GENUS_COLUMNS, NCBI_DATASET, OTT_DATASET, ChecklistBank, compare_ids, resolve_genera,
                           resolve_species)
from timegraft.webapi import ApiError, JsonCache

RELEASE = 316691


def _cls(*pairs):
    return [{"rank": r, "name": n} for r, n in pairs]


GASTROPOD = _cls(("kingdom", "Animalia"), ("phylum", "Mollusca"), ("class", "Gastropoda"),
                 ("order", "Systellommatophora"), ("family", "Onchidiidae"))
ECHINOID = _cls(("kingdom", "Metazoa"), ("phylum", "Echinodermata"), ("class", "Echinoidea"),
                ("family", "Laganidae"))
CILIATE = _cls(("kingdom", "Chromista"), ("phylum", "Ciliophora"), ("family", "Oxytrichidae"))
MITRIDAE = _cls(("kingdom", "Animalia"), ("phylum", "Mollusca"), ("class", "Gastropoda"), ("family", "Mitridae"))


class FakeCLB:
    """Routes ChecklistBank URLs to canned answers and counts the calls."""

    def __init__(self):
        self.calls = []

    def __call__(self, url, payload=None):
        self.calls.append(url)
        parsed = urlparse(url)
        path, q = parsed.path, {k: v[0] for k, v in parse_qs(parsed.query).items()}
        parts = path.strip("/").split("/")
        if path == f"/dataset/{RELEASE}" or path == "/dataset/3LXR":
            return {"key": RELEASE, "alias": "COL26.10 XR", "version": "2026-10-16 XR", "origin": "xrelease"}
        if path.endswith("/match/nameusage"):
            return self.match(q)
        if parts[-1] == "related":
            return self.related(parts[3], int(q["datasetKey"]))
        if parts[-1] == "classification":
            return {("2169", "36932"): GASTROPOD, ("2169", "262330"): ECHINOID, ("2169", "69650"): MITRIDAE,
                    ("201890", "162893"): GASTROPOD, ("201890", "101519"): MITRIDAE,
                    ("201890", "6383518"): CILIATE, ("2169", "555"): CILIATE, ("2169", "777"): MITRIDAE}.get((parts[1], parts[3]), [])
        if "sector" in parts:
            return {"mode": "merge" if parts[-1] == "99" else "attach"}
        raise ApiError(404, "not found")

    def match(self, q):
        name = q["q"]
        if name == "Onchidella":
            if q.get("family") == "Laganidae":  # wrong context: COL only matches a higher rank
                return {"type": "higherrank", "usage": {"id": "LAG", "rank": "family"}}
            return {"type": "variant", "usage": {"id": "7PB89", "name": name, "rank": "genus", "status": "accepted",
                                                 "sectorKey": 1302, "classification": GASTROPOD}}
        if name == "Mitra":
            if q.get("phylum") == "Mollusca":
                return {"type": "variant", "usage": {"id": "7P94F", "name": name, "rank": "genus",
                                                     "status": "accepted", "sectorKey": 1, "classification": MITRIDAE}}
            return {"type": "ambiguous", "alternatives": [
                {"id": "63C53", "rank": "genus", "status": "accepted", "classification": CILIATE},
                {"id": "7P94F", "rank": "genus", "status": "accepted", "classification": MITRIDAE}]}
        if name == "Boops":  # ChecklistBank says ambiguous, but the others are synonyms
            return {"type": "ambiguous", "usage": {"id": "B1", "name": name, "rank": "genus", "status": "accepted",
                                                   "sectorKey": 1, "classification": MITRIDAE},
                    "alternatives": [{"id": "B2", "rank": "genus", "status": "synonym"}]}
        if name == "Lysidice":  # two accepted animal genera; the record lineage tells them apart
            eunicid = _cls(("phylum", "Annelida"), ("class", "Polychaeta"), ("family", "Eunicidae"))
            eulepethid = _cls(("phylum", "Annelida"), ("class", "Polychaeta"), ("family", "Eulepethidae"))
            return {"type": "ambiguous", "usage": {"id": "L1", "rank": "genus", "status": "accepted",
                                                   "classification": eunicid},
                    "alternatives": [{"id": "L2", "rank": "genus", "status": "accepted",
                                      "classification": eulepethid}]}
        if name == "Xronly":
            return {"type": "exact", "usage": {"id": "XR1", "name": name, "rank": "genus", "status": "accepted",
                                               "sectorKey": 99, "classification": MITRIDAE}}
        if name == "Oldname":
            return {"type": "variant", "usage": {"id": "SYN", "name": name, "rank": "genus", "status": "synonym",
                                                 "accepted": {"id": "7P94F", "name": "Mitra"}, "sectorKey": 1,
                                                 "classification": MITRIDAE}}
        if name == "Mitra mitra":
            return {"type": "exact", "usage": {"id": "SP1", "name": name, "rank": "species", "status": "accepted",
                                               "sectorKey": 1, "classification": MITRIDAE}}
        return {"type": "none"}

    def related(self, usage, dataset):
        table = {
            ("7PB89", NCBI_DATASET): [{"id": "36932", "rank": "genus", "status": "accepted"}],
            ("7PB89", OTT_DATASET): [{"id": "162893", "rank": "genus", "status": "accepted"}],
            ("7P94F", NCBI_DATASET): [{"id": "69650", "rank": "genus", "status": "accepted"},
                                      {"id": "555", "rank": "genus", "status": "accepted"}],
            ("7P94F", OTT_DATASET): [{"id": "6383518", "rank": "genus", "status": "accepted"},
                                     {"id": "101519", "rank": "genus", "status": "accepted"}],
            ("XR1", NCBI_DATASET): [{"id": "555", "rank": "genus", "status": "accepted"}],
            ("SP1", NCBI_DATASET): [{"id": "777", "rank": "species", "status": "accepted"}],
        }
        return table.get((usage, dataset), [])


@pytest.fixture()
def clb():
    return ChecklistBank(RELEASE, JsonCache(None), fetch=FakeCLB())


def lineage(phylum="", family=""):
    return {"records_phylum": phylum, "records_family": family, "records_class": "", "records_order": ""}


def test_context_picks_the_right_homonym_and_its_ids(clb):
    rows = resolve_genera(clb, ["Onchidella"], {"Onchidella": lineage("Mollusca", "Onchidiidae")}, "Animalia")
    row = rows["Onchidella"]
    assert set(row) == set(GENUS_COLUMNS)
    assert row["col_id"] == "7PB89" and row["col_release"] == str(RELEASE) and row["col_context"] == "records"
    assert row["col_family"] == "Onchidiidae" and row["col_phylum"] == "Mollusca"
    assert row["col_ncbi_id"] == "36932" and row["col_ncbi_path"] == "names_index"
    assert row["col_ott_id"] == "162893" and row["col_confidence"] == "high" and row["col_source"] == "base"


def test_wrong_record_lineage_falls_back_to_the_kingdom(clb):
    row = resolve_genera(clb, ["Onchidella"], {"Onchidella": lineage("Echinodermata", "Laganidae")}, "Animalia")
    row = row["Onchidella"]
    assert row["col_context"] == "kingdom" and row["col_id"] == "7PB89"
    assert row["col_confidence"] == "medium"  # the records disagreed with COL


def test_ambiguity_is_kept(clb):
    row = resolve_genera(clb, ["Mitra"], {}, "Animalia")["Mitra"]
    assert row["col_match_type"] == "ambiguous" and row["col_id"] == "" and row["col_confidence"] == ""
    assert row["col_alternatives"] == "63C53;7P94F" and row["col_ncbi_id"] == ""


def test_synonym_only_alternatives_are_not_ambiguity(clb):
    row = resolve_genera(clb, ["Boops"], {"Boops": lineage("Mollusca")}, "Animalia")["Boops"]
    assert row["col_match_type"] == "ambiguous" and row["col_id"] == "B1" and row["col_alternatives"] == ""


def test_record_lineage_separates_accepted_homonyms(clb):
    lin = {"records_phylum": "Annelida", "records_class": "Polychaeta", "records_family": "Eunicidae",
           "records_order": ""}
    row = resolve_genera(clb, ["Lysidice"], {"Lysidice": lin}, "Animalia")["Lysidice"]
    assert row["col_id"] == "L1" and row["col_alternatives"] == "L2" and row["col_confidence"] == "medium"
    # without a lineage that tells them apart, it stays unresolved
    row = resolve_genera(clb, ["Lysidice"], {"Lysidice": lineage("Annelida")}, "Animalia")["Lysidice"]
    assert row["col_id"] == "" and row["col_alternatives"] == "L1;L2"


def test_classification_separates_homonyms_in_ncbi_and_ott(clb):
    row = resolve_genera(clb, ["Mitra"], {"Mitra": lineage("Mollusca")}, "Animalia")["Mitra"]
    assert row["col_ncbi_id"] == "69650" and row["col_ncbi_path"] == "names_index+classification"
    assert row["col_ncbi_candidates"] == "555;69650"
    assert row["col_ott_id"] == "101519"
    assert row["col_confidence"] == "medium"


def test_xr_only_and_disagreeing_crosswalk(clb):
    row = resolve_genera(clb, ["Xronly"], {}, "Animalia")["Xronly"]
    assert row["col_source"] == "xr" and row["col_confidence"] == "low"
    # the only NCBI Xronly is a ciliate: not joined
    assert row["col_ncbi_path"] == "disagrees" and row["col_ncbi_id"] == ""


def test_synonym_resolves_to_accepted(clb):
    row = resolve_genera(clb, ["Oldname"], {"Oldname": lineage("Mollusca")}, "Animalia")["Oldname"]
    assert row["col_status"] == "synonym" and row["col_id"] == "7P94F" and row["col_accepted_id"] == "7P94F"
    assert row["col_confidence"] == "medium"


def test_unknown_name(clb):
    row = resolve_genera(clb, ["Nothing"], {}, "Animalia")["Nothing"]
    assert row["col_match_type"] == "none" and row["col_id"] == "" and row["col_ncbi_path"] == ""


def test_answers_are_cached_per_release():
    fake = FakeCLB()
    cache = JsonCache(None)
    resolve_genera(ChecklistBank(RELEASE, cache, fetch=fake), ["Onchidella"], {}, "Animalia")
    n = len(fake.calls)
    resolve_genera(ChecklistBank(RELEASE, cache, fetch=fake), ["Onchidella"], {}, "Animalia")
    assert len(fake.calls) == n
    assert all(k.startswith("col:") for k in cache.data)


def test_alias_resolves_to_a_fixed_key():
    clb = ChecklistBank("3LXR", JsonCache(None), fetch=FakeCLB())
    assert clb.key == RELEASE and clb.info["alias"] == "COL26.10 XR"


def test_species(clb):
    rows = resolve_species(clb, {"Mitra mitra": "Mitra"}, {"Mitra": lineage("Mollusca")}, "Animalia",
                           {"Mitra mitra": 777})
    assert rows[0]["col_id"] == "SP1" and rows[0]["col_ncbi_id"] == "777" and rows[0]["col_ncbi_vs_ott"] == "same"


def test_compare_ids():
    assert compare_ids("36932", "262330;36932", "262330") == "merged"
    assert compare_ids("36932", "36932", "36932") == "same"
    assert compare_ids("1", "2;3", "2") == "different"
    assert compare_ids("", "2", "2") == "" and compare_ids("1", "", "") == ""


def test_join_ncbi_only_chooses_among_otts_ids():
    matches = {
        "Onchidella": {"ott_id": 162893, "ncbi_id": 262330, "status": "in_synthesis"},
        "Plain": {"ott_id": 2, "ncbi_id": 20, "status": "in_synthesis"},
        "Foreign": {"ott_id": 3, "ncbi_id": 30, "status": "in_synthesis"},
        "NoNcbi": {"ott_id": 4, "ncbi_id": None, "status": "attached_at_mrca"},
        "Weak": {"ott_id": 5, "ncbi_id": None, "status": "in_synthesis"},
    }
    info = {"Onchidella": {"ott_ncbi_ids": "262330;36932"}, "Plain": {"ott_ncbi_ids": "20"},
            "Foreign": {"ott_ncbi_ids": "30"}}
    col = {
        "Onchidella": {"col_ncbi_id": "36932", "col_confidence": "high"},
        "Plain": {"col_ncbi_id": "20", "col_confidence": "high"},
        "Foreign": {"col_ncbi_id": "99", "col_confidence": "high"},
        "NoNcbi": {"col_ncbi_id": "40", "col_confidence": "medium"},
        "Weak": {"col_ncbi_id": "50", "col_confidence": "low"},
    }
    ncbi = {"Onchidella": [262330, 1, 2], "Plain": [20], "Foreign": [30, 3], "NoNcbi": [0, 4], "Weak": [0]}
    out = join_ncbi(ncbi, matches, col, info, "col")
    assert ncbi == {"Onchidella": [36932, 1, 2], "Plain": [20], "Foreign": [30, 3], "NoNcbi": [40, 4], "Weak": [0]}
    assert out["Onchidella"] == {"join_ncbi_id": "36932", "join_ncbi_source": "col", "col_ncbi_vs_ott": "merged"}
    assert out["Foreign"]["col_ncbi_vs_ott"] == "different" and out["Foreign"]["join_ncbi_source"] == "ott"
    assert out["Weak"]["join_ncbi_source"] == ""
    before = {"Onchidella": [262330, 1, 2]}
    join_ncbi(before, {"Onchidella": matches["Onchidella"]}, col, info, "ott")
    assert before == {"Onchidella": [262330, 1, 2]}
    with pytest.raises(ValueError):
        join_ncbi({}, {}, {}, None, "nope")
