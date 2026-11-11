"""Faith's PD, MPD, UniFrac and PhyloSor for the samples in the sample table."""

import logging
import sys
from collections import defaultdict

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

import dendropy

from timegraft.diversity import alpha, beta, edge_table
from timegraft.graft import rescale_crowns
from timegraft.seqio import read_tsv, write_tsv

tree = dendropy.Tree.get(path=snakemake.input.tree, schema="newick", preserve_underscores=True,
                         rooting="force-rooted")
samples = defaultdict(dict)
for row in read_tsv(snakemake.input.samples):
    unit, abundance = row["unit_id"], float(row.get("abundance") or 1)
    samples[row["sample"]][unit] = samples[row["sample"]].get(unit, 0) + abundance

table, missing = edge_table(tree, samples)
alpha_rows = alpha(table)

# The same measures with every genus crown at its lower and at its upper bound:
# how much of a value rests on the crown ages (#36). Rescaling the subtrees
# gives the trees that runs with those crown ages would give.
genus_units = {}
for path in snakemake.input.subtrees:
    for row in read_tsv(path):
        sub = dendropy.Tree.get(data=row["newick"], schema="newick", preserve_underscores=True)
        genus_units[row["genus"]] = [leaf.taxon.label for leaf in sub.leaf_node_iter()]
provenance = [r for r in read_tsv(snakemake.input.provenance) if r["status"] == "subtree"]
for bound, column in (("min", "crown_lower"), ("max", "crown_upper")):
    # no lower bound (markers without a rate): the genus collapses to its stem
    ages = {r["genus"]: float(r[column] or 0) for r in provenance}
    moved = dendropy.Tree.get(path=snakemake.input.tree, schema="newick", preserve_underscores=True,
                              rooting="force-rooted")
    rescale_crowns(moved, genus_units, ages)
    rows = {r["sample"]: r for r in alpha(edge_table(moved, samples)[0])}
    for row in alpha_rows:
        row[f"faith_pd_crown_{bound}"] = rows[row["sample"]]["faith_pd"]
        row[f"mpd_crown_{bound}"] = rows[row["sample"]]["mpd"]
write_tsv(snakemake.output.alpha, alpha_rows,
          ["sample", "richness", "faith_pd", "faith_pd_crown_min", "faith_pd_crown_max",
           "mpd", "mpd_crown_min", "mpd_crown_max", "mpd_abundance"])
write_tsv(snakemake.output.beta, beta(table),
          ["sample_a", "sample_b", "unifrac_unweighted", "unifrac_weighted", "phylosor"])
write_tsv(snakemake.output.missing,
          [{"sample": s, "unit_id": u} for s, units in sorted(missing.items()) for u in units],
          ["sample", "unit_id"])
logging.info("%d samples; %d unit occurrences not in the tree", len(samples),
             sum(len(u) for u in missing.values()))
