"""End-to-end run of the Snakemake workflow on the toy data."""

import shutil
import subprocess
from pathlib import Path

import dendropy
import pytest

from timegraft.dating import is_ultrametric, root_age
from timegraft.seqio import read_tsv

REPO = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(shutil.which("snakemake") is None or shutil.which("vsearch") is None,
                    reason="snakemake or vsearch not installed")
def test_toy_workflow(tmp_path):
    outdir = tmp_path / "out"
    subprocess.run(
        [
            "snakemake",
            "--cores", "2",
            "--quiet", "rules",
            "--snakefile", str(REPO / "workflow" / "Snakefile"),
            "--directory", str(REPO),
            "--config", f"outdir={outdir}",
        ],
        check=True,
    )
    tree = dendropy.Tree.get(
        path=str(outdir / "timegraft.tre"), schema="newick", preserve_underscores=True, rooting="force-rooted"
    )
    tips = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    # the ten units with a genus, plus the unit without a genus that was placed
    assert tips == {f"BOLD:AAA{i:04d}" for i in range(1, 11)} | {"BOLD:AAZ0002"}
    units_placed = {r["unit_id"]: r for r in read_tsv(outdir / "units_placed.tsv")}
    assert units_placed["BOLD:AAZ0002"]["status"] == "placed"
    assert units_placed["BOLD:AAZ0002"]["nearest_units"] == "BOLD:AAA0009"
    assert units_placed["BOLD:AAZ0001"]["status"] == "unplaced"
    assert is_ultrametric(tree, tolerance=1e-6)
    assert root_age(tree) == pytest.approx(50.0)

    status = {r["genus"]: r for r in read_tsv(outdir / "provenance.tsv")}
    assert status["Alpha"]["crown_source"] == "table"
    assert status["Alpha"]["crown_age"] == "8.0000"
    assert status["Delta"]["status"] == "single_unit"
    assert status["Epsilon"]["status"] == "pruned"
    for genus in ("Alpha", "Beta", "Gamma"):
        assert status[genus]["rooting"] == "outgroup"
        assert float(status[genus]["crown_age"]) < float(status[genus]["stem_age"])
        assert int(status[genus]["root_alternatives"]) >= 1
        assert float(status[genus]["pd_rel_range"]) >= 0.0
    # root uncertainty (#34): crown age ranges only where the age comes from the barcodes
    assert status["Alpha"]["crown_age_min"] == ""
    for genus in ("Beta", "Gamma"):
        if status[genus]["crown_source"].startswith("rate"):
            row = status[genus]
            assert float(row["crown_age_min"]) <= float(row["crown_age"]) <= float(row["crown_age_max"])
    assert status["Delta"]["root_alternatives"] == ""

    placements = {r["asv_id"]: r for r in read_tsv(outdir / "placement" / "placements.tsv")}
    assert {a: r["status"] for a, r in placements.items()} == {
        "asv_01": "matched", "asv_02": "matched", "asv_03": "placed", "asv_04": "unplaced", "asv_05": "placed",
    }
    assert placements["asv_03"]["unit_id"] == "BOLD:AAA0008"
    placed = dendropy.Tree.get(path=str(outdir / "placement" / "timegraft_placed.tre"), schema="newick",
                               preserve_underscores=True, rooting="force-rooted")
    assert is_ultrametric(placed, tolerance=1e-6)
    assert {"asv_03", "asv_05"} <= {l.taxon.label for l in placed.leaf_node_iter()}

    alpha = {r["sample"]: r for r in read_tsv(outdir / "diversity" / "alpha.tsv")}
    assert set(alpha) == {"estuary_1", "estuary_2", "marina_1"}
    assert alpha["estuary_2"]["richness"] == "3"
    assert len(read_tsv(outdir / "diversity" / "beta.tsv")) == 3

    summary = read_tsv(outdir / "validation" / "summary.tsv")
    assert {(r["method"], r["stratum"]) for r in summary} >= {("k2p", "all"), ("kmer", "genus")}
    assert (outdir / "validation" / "saturation.png").stat().st_size > 0

    rejected = {r["record_id"]: r for r in read_tsv(outdir / "units" / "rejected.tsv")}
    assert {k: r["reason"] for k, r in rejected.items()} == {
        "BAD_stop": "stop_codon", "BAD_short": "too_short",
        "BAD_vertebrate_stop": "stop_codon", "BAD_standard_stop": "stop_codon",
    }
    # genetic code per record: from the shipped mapping, the record's own column, or the default
    assert (rejected["BAD_vertebrate_stop"]["genetic_code"], rejected["BAD_vertebrate_stop"]["genetic_code_source"]) == ("2", "table")
    assert (rejected["BAD_standard_stop"]["genetic_code"], rejected["BAD_standard_stop"]["genetic_code_source"]) == ("1", "column")
    assert (rejected["BAD_stop"]["genetic_code"], rejected["BAD_stop"]["genetic_code_source"]) == ("5", "default")
    assert {r["genetic_code"] for r in read_tsv(outdir / "units" / "filtered.tsv")} == {"5"}
