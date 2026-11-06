"""BCDM data package conversion."""

from timegraft.bcdm import convert, main
from timegraft.seqio import read_tsv

HEADER = ["processid", "sampleid", "bin_uri", "phylum", "class", "order", "family", "subfamily",
          "genus", "species", "country/ocean", "marker_code", "nuc"]
ROWS = [
    ["P1", "s1", "BOLD:AAA0001", "Mollusca", "Bivalvia", "Mytilida", "Mytilidae", "None", "Mytilus",
     "Mytilus edulis", "Portugal", "COI-5P", "acgt-acgt"],
    ["P2", "s2", "BOLD:AAA0002", "Annelida", "Polychaeta", "Phyllodocida", "Nereididae", "None", "Hediste",
     "None", "Spain", "COI-5P", "ACGTACGT"],
    ["P3", "s3", "BOLD:AAA0003", "Arthropoda", "Insecta", "Diptera", "Culicidae", "None", "Aedes",
     "None", "Portugal", "COI-5P", "ACGTACGT"],
    ["P4", "s4", "", "Mollusca", "Bivalvia", "Mytilida", "Mytilidae", "None", "Mytilus",
     "None", "Portugal", "COI-5P", "ACGTACGT"],
    ["P5", "s5", "BOLD:AAA0005", "Mollusca", "Bivalvia", "Mytilida", "Mytilidae", "None", "Mytilus",
     "None", "France", "COI-5P", "ACGTACGT"],
    ["P6", "s6", "BOLD:AAA0006", "Mollusca", "Bivalvia", "Mytilida", "Mytilidae", "None", "Mytilus",
     "None", "Portugal", "28S", "ACGTACGT"],
]


def test_convert_filters_and_maps(tmp_path):
    source = tmp_path / "bold.tsv"
    source.write_text("\n".join("\t".join(r) for r in [HEADER] + ROWS) + "\n")
    target = tmp_path / "records.tsv"
    counts = convert(source, target, taxa={"phylum": {"Mollusca", "Annelida"}}, countries={"Portugal", "Spain"})
    rows = read_tsv(target)
    assert [r["record_id"] for r in rows] == ["P1", "P2"]
    assert rows[0]["sequence"] == "ACGT-ACGT"
    assert rows[1]["species"] == ""  # 'None' is cleaned
    assert counts == {"read": 6, "kept": 2, "other_marker": 1, "no_unit": 1, "no_sequence": 0,
                      "taxon": 1, "country": 1}


def test_command_line(tmp_path, capsys):
    source = tmp_path / "bold.tsv"
    source.write_text("\n".join("\t".join(r) for r in [HEADER] + ROWS) + "\n")
    main([str(source), str(tmp_path / "out.tsv"), "--taxon", "family=Mytilidae", "--country", "Portugal"])
    assert "kept=1" in capsys.readouterr().out
