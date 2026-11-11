"""Group representatives by genus, pick sister-genus outgroups, and split into shards."""

import logging
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

from timegraft.backbone import pick_outgroup, read_backbone, shard_of, stem_ages
from timegraft.rates import RANKS, genus_taxonomy
from timegraft.seqio import read_tsv, write_tsv

reps = read_tsv(snakemake.input.reps)
backbone = read_backbone(snakemake.input.backbone)
stems = stem_ages(backbone)
n_shards = snakemake.params.n_shards
rooting = snakemake.params.rooting
n_outgroup = int(rooting.get("n_outgroup", 10))
choice = rooting.get("outgroup_choice", "spread")
levels = int(rooting.get("outgroup_levels", 3))

by_genus = defaultdict(list)
for rep in reps:
    by_genus[rep["genus"]].append(rep)

missing = [{"genus": g, "n_units": len(u)} for g, u in sorted(by_genus.items()) if g not in stems]
available = {g for g in by_genus if g in stems}

shards = defaultdict(list)
singles = []
for genus in sorted(available):
    units = sorted(by_genus[genus], key=lambda r: r["unit_id"])
    if len(units) == 1:
        singles.append({"genus": genus, "unit_id": units[0]["unit_id"]})
        continue
    rows = shards[shard_of(genus, n_shards)]
    for unit in units:
        rows.append({"genus": genus, "role": "in", "unit_id": unit["unit_id"],
                     "stem_age": stems[genus], "sequence": unit["sequence"]})
    picked = pick_outgroup(backbone, genus, by_genus, available, n_outgroup, choice, levels)
    for unit in picked:
        rows.append({"genus": genus, "role": "out", "unit_id": unit["unit_id"],
                     "stem_age": stems[genus], "sequence": unit["sequence"]})

columns = ["genus", "role", "unit_id", "stem_age", "sequence"]
for i, path in enumerate(snakemake.output.shards):
    write_tsv(path, shards.get(i, []), columns)
write_tsv(snakemake.output.singles, singles, ["genus", "unit_id"])
write_tsv(snakemake.output.missing, missing, ["genus", "n_units"])
# higher taxa per genus (the most common value among its records), for rates per taxon
taxonomy = genus_taxonomy(read_tsv(snakemake.input.records))
write_tsv(snakemake.output.taxonomy, [taxonomy[g] for g in sorted(taxonomy) if g in by_genus], RANKS)
logging.info("%d genera with subtrees, %d single-unit genera, %d genera missing from backbone",
             len(available) - len(singles), len(singles), len(missing))
