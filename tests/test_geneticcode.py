"""Genetic codes: tables, translation, resolution per record, and the taxdump builder."""

import pytest

from timegraft import filters, geneticcode
from timegraft.geneticcode import GeneticCodes, MIXED


# tables and translation --------------------------------------------------------

def test_tables_match_biopython():
    CodonTable = pytest.importorskip("Bio.Data.CodonTable")
    for code, aas in geneticcode.TABLES.items():
        table = CodonTable.unambiguous_dna_by_id[code]
        for codon, aa in zip(geneticcode._CODONS, aas):
            expected = "*" if codon in table.stop_codons else table.forward_table[codon]
            assert aa == expected, (code, codon)


def test_stop_codons():
    assert geneticcode.stop_codons(1) == {"TAA", "TAG", "TGA"}
    assert geneticcode.stop_codons(2) == {"TAA", "TAG", "AGA", "AGG"}
    assert geneticcode.stop_codons(5) == geneticcode.stop_codons(9) == geneticcode.stop_codons(4) == {"TAA", "TAG"}
    assert geneticcode.stop_codons(14) == {"TAG"}


def test_translate_differs_between_codes():
    assert geneticcode.translate("AAAAGATGA", 5) == "KSW"
    assert geneticcode.translate("AAAAGATGA", 9) == "NSW"
    assert geneticcode.translate("AAAAGATGA", 2) == "K*W"
    assert geneticcode.translate("AAAAGATGA", 13) == "KGW"
    assert geneticcode.translate("AAA-AGN", 1) == "KX"
    assert geneticcode.translate("AAAAGATGA", 1, frame=1) == "KD"


def test_best_frame():
    seq = "C" + "ATG" * 10 + "TAA" + "ATG"
    frame, n = geneticcode.best_frame(seq, 5)
    assert n == 0 and frame != 1
    assert geneticcode.frame_stops(seq, 5)[1] == 1


def test_unsupported_code():
    with pytest.raises(ValueError):
        geneticcode.translate("ATG", 6)


def test_vertebrate_and_flatworm_frames():
    vertebrate_stop = ("AGACC" * 30)[:150]  # AGA in every frame
    assert filters.has_open_frame(vertebrate_stop, genetic_code=5)
    assert not filters.has_open_frame(vertebrate_stop, genetic_code=2)
    taa = ("TAACC" * 30)[:150]  # TAA in every frame: Tyr under 14
    assert not filters.has_open_frame(taa, genetic_code=5)
    assert filters.has_open_frame(taa, genetic_code=14)


# resolution per record -----------------------------------------------------------

def test_shipped_table_resolves_common_groups():
    codes = GeneticCodes()
    assert codes.resolve({"phylum": "Arthropoda", "class": "Insecta"})[:2] == (5, "table")
    assert codes.resolve({"phylum": "Chordata", "class": "Actinopterygii"})[:2] == (2, "table")
    assert codes.resolve({"phylum": "Chordata", "class": "Ascidiacea"})[:2] == (13, "table")
    assert codes.resolve({"phylum": "Echinodermata"})[0] == 9
    assert codes.resolve({"phylum": "Cnidaria"})[0] == 4
    assert codes.resolve({"phylum": "Platyhelminthes", "class": "Trematoda"})[0] == 21


def test_most_specific_name_wins():
    codes = GeneticCodes()
    code, source, taxon = codes.resolve({"phylum": "Chordata", "class": "Reptilia", "order": "Squamata"})
    assert (code, source, taxon) == (2, "table", "Squamata")


def test_mixed_names_fall_through_to_default():
    codes = GeneticCodes(default=5)
    assert codes.resolve({"phylum": "Chordata"}) == (5, "default", "")
    assert codes.resolve({"genus": "Unheardof"}) == (5, "default", "")
    assert codes.resolve({}) == (5, "default", "")


def test_column_beats_everything():
    codes = GeneticCodes(overrides={"class": {"Actinopterygii": 5}})
    assert codes.resolve({"class": "Actinopterygii", "genetic_code": "1"}) == (1, "column", "")
    assert codes.resolve({"class": "Actinopterygii", "genetic_code": " "})[:2] == (5, "config")


def test_config_beats_table_at_the_same_rank_but_not_a_more_specific_one():
    codes = GeneticCodes(overrides={"phylum": {"Chordata": 2}})
    assert codes.resolve({"phylum": "Chordata"})[:2] == (2, "config")
    assert codes.resolve({"phylum": "Chordata", "class": "Ascidiacea"})[:2] == (13, "table")


def test_own_table_and_validation(tmp_path):
    path = tmp_path / "codes.tsv"
    path.write_text("# comment\nname\trank\tcode\tcodes\tncbi_taxid\nFoo\tphylum\t9\t\t\nBar\tphylum\tmixed\t\t\n")
    codes = GeneticCodes(table=path)
    assert codes.resolve({"phylum": "Foo"})[:2] == (9, "table")
    assert codes.resolve({"phylum": "Bar"})[1] == "default"
    assert codes.resolve({"phylum": "Arthropoda"})[1] == "default"
    assert GeneticCodes(table=False).resolve({"phylum": "Cnidaria"})[1] == "default"
    with pytest.raises(ValueError):
        GeneticCodes(overrides={"tribe": {"X": 5}})
    with pytest.raises(ValueError):
        GeneticCodes(overrides={"class": {"X": 6}})
    with pytest.raises(ValueError):
        GeneticCodes().resolve({"genetic_code": "99"})


# building the mapping from a taxdump -------------------------------------------------

def _write_taxdump(directory):
    # taxid, parent, rank, mitochondrial code
    nodes = [
        (1, 1, "no rank", 0),
        (2759, 1, "superkingdom", 1),
        (33208, 2759, "kingdom", 1),
        (6656, 33208, "phylum", 5),         # Arthropoda
        (50557, 6656, "class", 5),          # Insecta
        (7711, 33208, "phylum", 5),         # Chordata: mixed below
        (7742, 7711, "clade", 2),           # Vertebrata (not an emitted rank)
        (7898, 7742, "superclass", 2),      # Actinopterygii
        (7777, 7742, "class", 2),           # Chondrichthyes
        (7778, 7777, "subclass", 2),        # Elasmobranchii
        (7712, 7711, "subphylum", 13),      # Tunicata
        (7713, 7712, "class", 13),          # Ascidiacea
        (99, 6656, "genus", 5),
        (98, 99, "species", 1),             # a species code does not make its phylum mixed
        (2763, 2759, "phylum", 4),          # an algal phylum also called "Insecta" (homonym)
    ]
    names = [
        (6656, "Arthropoda", "scientific name"), (50557, "Insecta", "scientific name"),
        (7711, "Chordata", "scientific name"), (7742, "Vertebrata", "scientific name"),
        (7898, "Actinopterygii", "scientific name"), (7777, "Chondrichthyes", "scientific name"),
        (7778, "Elasmobranchii", "scientific name"), (7712, "Tunicata", "scientific name"),
        (7712, "Urochordata", "synonym"), (7713, "Ascidiacea", "scientific name"),
        (7713, "sea squirts", "common name"), (33208, "Metazoa", "scientific name"),
        (2763, "Rhodophyta", "scientific name"), (2763, "Insecta", "synonym"),
    ]
    with open(directory / "nodes.dmp", "w") as handle:
        for taxid, parent, rank, mgc in nodes:
            fields = [taxid, parent, rank, "", 0, 0, 1, 0, mgc, 0, 0, 0, ""]
            handle.write("\t|\t".join(map(str, fields)) + "\t|\n")
    with open(directory / "names.dmp", "w") as handle:
        for taxid, name, kind in names:
            handle.write("\t|\t".join([str(taxid), name, "", kind]) + "\t|\n")


def test_build_table(tmp_path):
    _write_taxdump(tmp_path)
    rows = {r["name"]: r for r in geneticcode.build_table(tmp_path)}
    assert rows["Arthropoda"]["code"] == "5"
    assert rows["Chordata"]["code"] == MIXED and rows["Chordata"]["codes"] == "2;5;13"
    assert rows["Elasmobranchii"]["code"] == "2" and rows["Elasmobranchii"]["rank"] == "subclass"
    assert rows["Actinopterygii"]["code"] == "2"
    assert rows["Urochordata"]["code"] == rows["Ascidiacea"]["code"] == "13"
    assert rows["Insecta"]["code"] == MIXED  # homonym with different codes
    assert rows["Metazoa"]["code"] == MIXED
    assert "Vertebrata" not in rows and "sea squirts" not in rows

    out = tmp_path / "codes.tsv"
    geneticcode.write_table(list(rows.values()), out, "test")
    codes = GeneticCodes(table=out)
    assert codes.resolve({"phylum": "Chordata", "class": "Elasmobranchii"})[:2] == (2, "table")
    assert codes.resolve({"phylum": "Arthropoda", "class": "Insecta"})[:2] == (5, "table")


def test_build_table_cli_on_a_tarball(tmp_path):
    import tarfile

    _write_taxdump(tmp_path)
    archive = tmp_path / "taxdump.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(tmp_path / "nodes.dmp", arcname="nodes.dmp")
        tar.add(tmp_path / "names.dmp", arcname="names.dmp")
    out = tmp_path / "out.tsv"
    geneticcode.main([str(archive), str(out), "--root", "33208"])
    names = {line.split("\t")[0] for line in out.read_text().splitlines() if not line.startswith("#")}
    assert "Rhodophyta" not in names and "Arthropoda" in names
