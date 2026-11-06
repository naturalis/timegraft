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
                      "taxon": 1, "country": 1, "unit_cap": 0}


def test_command_line(tmp_path, capsys):
    source = tmp_path / "bold.tsv"
    source.write_text("\n".join("\t".join(r) for r in [HEADER] + ROWS) + "\n")
    main([str(source), str(tmp_path / "out.tsv"), "--taxon", "family=Mytilidae", "--country", "Portugal"])
    assert "kept=1" in capsys.readouterr().out


def test_parquet_gives_the_same_records(tmp_path):
    duckdb = __import__("pytest").importorskip("duckdb")
    source = tmp_path / "bold.tsv"
    source.write_text("\n".join("\t".join(r) for r in [HEADER] + ROWS) + "\n")
    parquet = tmp_path / "bold.parquet"
    duckdb.execute(
        f"COPY (SELECT * FROM read_csv('{source}', delim='\t', header=true, all_varchar=true)) TO '{parquet}' (FORMAT parquet)"
    )
    filters = dict(taxa={"phylum": {"Mollusca", "Annelida"}}, countries={"Portugal", "Spain"})
    convert(source, tmp_path / "from_tsv.tsv", **filters)
    counts = convert(parquet, tmp_path / "from_parquet.tsv", **filters)
    assert read_tsv(tmp_path / "from_parquet.tsv") == read_tsv(tmp_path / "from_tsv.tsv")
    assert counts == {"read": 6, "kept": 2, "filtered_out": 4}


def test_exclude_and_cap_per_unit(tmp_path):
    rows = ROWS + [["P7", "s7", "BOLD:AAA0001", *ROWS[0][3:]], ["P8", "s8", "BOLD:AAA0001", *ROWS[0][3:]]]
    source = tmp_path / "bold.tsv"
    source.write_text("\n".join("\t".join(r) for r in [HEADER] + rows) + "\n")
    names = tmp_path / "phyla.txt"
    names.write_text("Arthropoda\n")
    main([str(source), str(tmp_path / "a.tsv"), "--exclude", f"phylum=@{names}", "--max-per-unit", "2"])
    kept = read_tsv(tmp_path / "a.tsv")
    assert [r["record_id"] for r in kept] == ["P1", "P2", "P5", "P7"]  # P3 is an insect, P8 is a third record of AAA0001
    duckdb = __import__("pytest").importorskip("duckdb")
    parquet = tmp_path / "bold.parquet"
    duckdb.execute(
        f"COPY (SELECT * FROM read_csv('{source}', delim='\t', header=true, all_varchar=true)) TO '{parquet}' (FORMAT parquet)"
    )
    main([str(parquet), str(tmp_path / "b.tsv"), "--exclude", f"phylum=@{names}", "--max-per-unit", "2"])
    assert sorted(r["record_id"] for r in read_tsv(tmp_path / "b.tsv")) == ["P1", "P2", "P5", "P7"]
