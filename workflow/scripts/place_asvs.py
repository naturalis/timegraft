"""Match ASVs to units, place the rest on the tree, and translate the ASV table to tips."""

import logging
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

import dendropy

from timegraft.placement import PLACEMENT_COLUMNS, classify, place_all, vsearch_hits
from timegraft.seqio import read_fasta, read_tsv, write_tsv

cfg = snakemake.params.settings
tree = dendropy.Tree.get(path=snakemake.input.tree, schema="newick", preserve_underscores=True,
                         rooting="force-rooted")
in_tree = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
records = read_tsv(snakemake.input.records)
record_unit = {r["record_id"]: r["unit_id"] for r in records}
unit_genus = {r["unit_id"]: r["genus"] for r in read_tsv(snakemake.input.reps)}
refs = [(r["record_id"], r["sequence"]) for r in records if r["unit_id"] in in_tree]
asvs = list(read_fasta(snakemake.input.asvs))

hits = vsearch_hits(
    asvs, refs,
    min_identity=cfg["place_identity"],
    query_cov=cfg.get("query_cov", 0.9),
    maxaccepts=cfg.get("maxaccepts", 20),
    threads=snakemake.threads,
)
rows = classify(hits, [a for a, _ in asvs], record_unit, unit_genus, in_tree, cfg["match_identity"])
placements = place_all(tree, rows, height=cfg.get("height", "rate"), rate=snakemake.params.rate,
                       terminal_fraction=cfg.get("terminal_fraction", 0.5))

tip_of = {p["asv_id"]: p["tip"] for p in placements if p["tip"]}
counts = defaultdict(float)
unplaced_reads = 0.0
for row in read_tsv(snakemake.input.asv_table):
    reads = float(row["reads"])
    tip = tip_of.get(row["asv_id"])
    if tip:
        counts[(row["sample"], tip)] += reads
    else:
        unplaced_reads += reads

Path(snakemake.output.tree).write_text(
    tree.as_string(schema="newick", suppress_rooting=True, preserve_spaces=True)
)
write_tsv(snakemake.output.placements, placements, PLACEMENT_COLUMNS)
write_tsv(snakemake.output.samples,
          [{"sample": s, "unit_id": t, "abundance": a} for (s, t), a in sorted(counts.items())],
          ["sample", "unit_id", "abundance"])
status = defaultdict(int)
for p in placements:
    status[p["status"]] += 1
logging.info("ASVs by status: %s; reads of unplaced ASVs: %d", dict(status), unplaced_reads)
