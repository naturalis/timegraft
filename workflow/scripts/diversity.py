"""Faith's PD, MPD, UniFrac and PhyloSor for the samples in the sample table."""

import logging
import sys
from collections import defaultdict

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

import dendropy

from timegraft.diversity import alpha, beta, edge_table
from timegraft.seqio import read_tsv, write_tsv

tree = dendropy.Tree.get(path=snakemake.input.tree, schema="newick", preserve_underscores=True,
                         rooting="force-rooted")
samples = defaultdict(dict)
for row in read_tsv(snakemake.input.samples):
    unit, abundance = row["unit_id"], float(row.get("abundance") or 1)
    samples[row["sample"]][unit] = samples[row["sample"]].get(unit, 0) + abundance

table, missing = edge_table(tree, samples)
alpha_rows = alpha(table)
write_tsv(snakemake.output.alpha, alpha_rows, ["sample", "richness", "faith_pd", "mpd", "mpd_abundance"])
write_tsv(snakemake.output.beta, beta(table),
          ["sample_a", "sample_b", "unifrac_unweighted", "unifrac_weighted", "phylosor"])
write_tsv(snakemake.output.missing,
          [{"sample": s, "unit_id": u} for s, units in sorted(missing.items()) for u in units],
          ["sample", "unit_id"])
logging.info("%d samples; %d unit occurrences not in the tree", len(samples),
             sum(len(u) for u in missing.values()))
