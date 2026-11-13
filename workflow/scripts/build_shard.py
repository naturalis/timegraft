"""Build a dated subtree for every genus in one shard."""

import logging
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

from timegraft.build import build_genus, tree_to_newick
from timegraft.calibrate import read_crown_table
from timegraft.graft import PROVENANCE_COLUMNS
from timegraft.rates import Rates
from timegraft.seqio import read_tsv, write_tsv

# crown ages: a point age in the user's table is used as is; its minimum and
# maximum, the TimeTree estimates and the stem age bound the age from the barcodes
known, timetree_crowns, crown_basis = {}, {}, {}
if snakemake.input.crown_ages:
    known = read_crown_table(snakemake.input.crown_ages)
if snakemake.input.timetree_crowns:
    for r in read_tsv(snakemake.input.timetree_crowns):
        if r.get("crown_age"):  # empty for a genus sampled with a single species (#54)
            timetree_crowns[r["genus"]] = float(r["crown_age"])
            crown_basis[r["genus"]] = r.get("crown_basis") or "global"

groups = defaultdict(lambda: {"in": [], "out": [], "stem": None})
for row in read_tsv(snakemake.input.shard):
    group = groups[row["genus"]]
    group[row["role"]].append((row["unit_id"], row["sequence"]))
    group["stem"] = float(row["stem_age"])

marker_cfg = dict(snakemake.params.marker_cfg)
crown_cfg = marker_cfg["crown_age"]
if "timetree_max_ratio" in crown_cfg:
    logging.warning("crown_age.timetree_max_ratio is no longer used: TimeTree crown ages are now an upper "
                    "bound on the age from the barcodes (crown_age.timetree_bound), not accepted or rejected")
# the substitution rate per genus, from its taxonomy (rate method only)
rates = Rates(crown_cfg.get("rate"), crown_cfg.get("rates"), snakemake.input.rate_table or None)
taxonomy = {r["genus"]: r for r in read_tsv(snakemake.input.taxonomy)}
if snakemake.input.lsd2:
    marker_cfg["dating"] = {**marker_cfg.get("dating", {}), "executable": str(Path(snakemake.input.lsd2).resolve())}

out = []
for genus in sorted(groups):
    group = groups[genus]
    tree, provenance = build_genus(genus, group["in"], group["out"], group["stem"],
                                   marker_cfg, known, timetree_crowns, rooting=snakemake.params.rooting,
                                   rate=rates.resolve(taxonomy.get(genus)), crown_basis=crown_basis)
    out.append({**provenance, "newick": tree_to_newick(tree)})
    logging.info("%s: %d units, crown %s Ma (%s), rooting %s (%s misplaced)", genus, provenance["n_units"],
                 provenance["crown_age"], provenance["crown_source"], provenance["rooting"],
                 provenance["violations"] or 0)

write_tsv(snakemake.output[0], out, list(PROVENANCE_COLUMNS) + ["newick"])
