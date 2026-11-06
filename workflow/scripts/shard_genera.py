"""Group representatives by genus, pick sister-genus outgroups, and split into shards."""

import logging
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

from timegraft.backbone import nearest_genera, read_backbone, shard_of, stem_ages
from timegraft.seqio import read_tsv, write_tsv

reps = read_tsv(snakemake.input.reps)
backbone = read_backbone(snakemake.input.backbone)
stems = stem_ages(backbone)
n_shards = snakemake.params.n_shards
n_outgroup = snakemake.params.n_outgroup

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
    picked = []
    for sister in nearest_genera(backbone, genus, available):
        for unit in sorted(by_genus[sister], key=lambda r: r["unit_id"]):
            if len(picked) < n_outgroup:
                picked.append(unit)
    for unit in picked:
        rows.append({"genus": genus, "role": "out", "unit_id": unit["unit_id"],
                     "stem_age": stems[genus], "sequence": unit["sequence"]})

columns = ["genus", "role", "unit_id", "stem_age", "sequence"]
for i, path in enumerate(snakemake.output.shards):
    write_tsv(path, shards.get(i, []), columns)
write_tsv(snakemake.output.singles, singles, ["genus", "unit_id"])
write_tsv(snakemake.output.missing, missing, ["genus", "n_units"])
logging.info("%d genera with subtrees, %d single-unit genera, %d genera missing from backbone",
             len(available) - len(singles), len(singles), len(missing))
