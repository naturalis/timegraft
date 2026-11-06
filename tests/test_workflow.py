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
    assert tips == {f"BOLD:AAA{i:04d}" for i in range(1, 11)}
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

    rejected = {r["record_id"]: r["reason"] for r in read_tsv(outdir / "units" / "rejected.tsv")}
    assert rejected == {"BAD_stop": "stop_codon", "BAD_short": "too_short"}
