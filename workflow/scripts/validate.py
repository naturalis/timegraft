"""Compare barcode distances with a reference tree, per taxonomic stratum."""

import logging
import sys
from collections import Counter, defaultdict

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

import dendropy

from timegraft.seqio import read_tsv, write_tsv
from timegraft.validation import compare, plot

cfg = snakemake.params.settings
ranks = list(cfg.get("ranks", ["genus", "family", "order"]))
methods = cfg.get("methods", {"kmer": {"method": "kmer", "k": 8}})

tree = dendropy.Tree.get(path=snakemake.input.reference, schema="newick", preserve_underscores=True,
                         rooting="force-rooted")

# rank values per unit: the majority among its records
votes = defaultdict(lambda: defaultdict(Counter))
for record in read_tsv(snakemake.input.records):
    for rank in ranks:
        if record.get(rank):
            votes[record["unit_id"]][rank][record[rank]] += 1
units = {}
for rep in read_tsv(snakemake.input.reps):
    info = {"sequence": rep["sequence"]}
    for rank in ranks:
        counter = votes[rep["unit_id"]][rank]
        info[rank] = counter.most_common(1)[0][0] if counter else ""
    units[rep["unit_id"]] = info

rows, summary, bins = compare(tree, units, methods, ranks, per_stratum=int(cfg.get("per_stratum", 200)),
                              seed=int(cfg.get("seed", 1)))
write_tsv(snakemake.output.pairs, rows, ["unit_a", "unit_b", "stratum", "reference", *methods])
write_tsv(snakemake.output.summary, summary,
          ["method", "stratum", "n_pairs", "spearman", "plateau_onset"])
write_tsv(snakemake.output.bins, [{"method": m, **b} for m, bs in bins.items() for b in bs],
          ["method", "reference_median", "median", "q25", "q75", "n"])
plot(rows, bins, ranks, snakemake.output.plot)
logging.info("%d pairs compared with %d methods", len(rows), len(methods))
