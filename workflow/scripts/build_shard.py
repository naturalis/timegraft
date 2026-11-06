"""Build a dated subtree for every genus in one shard."""

import logging
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

from timegraft.build import build_genus, tree_to_newick
from timegraft.graft import PROVENANCE_COLUMNS
from timegraft.seqio import read_tsv, write_tsv

known = {}
if snakemake.input.crown_ages:
    known = {r["genus"]: float(r["crown_age"]) for r in read_tsv(snakemake.input.crown_ages)}

groups = defaultdict(lambda: {"in": [], "out": [], "stem": None})
for row in read_tsv(snakemake.input.shard):
    group = groups[row["genus"]]
    group[row["role"]].append((row["unit_id"], row["sequence"]))
    group["stem"] = float(row["stem_age"])

out = []
for genus in sorted(groups):
    group = groups[genus]
    tree, provenance = build_genus(genus, group["in"], group["out"], group["stem"],
                                   snakemake.params.marker_cfg, known)
    out.append({**provenance, "newick": tree_to_newick(tree)})
    logging.info("%s: %d units, crown %s Ma (%s), rooting %s", genus, provenance["n_units"],
                 provenance["crown_age"], provenance["crown_source"], provenance["rooting"])

write_tsv(snakemake.output[0], out, list(PROVENANCE_COLUMNS) + ["newick"])
