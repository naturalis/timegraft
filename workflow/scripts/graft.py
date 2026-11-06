"""Graft all genus subtrees onto the backbone and write the final tree."""

import logging
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

from timegraft.backbone import read_backbone, stem_ages
from timegraft.dating import is_ultrametric, root_age
from timegraft.graft import PROVENANCE_COLUMNS, graft
from timegraft.seqio import read_tsv, write_tsv

backbone = read_backbone(snakemake.input.backbone)
stems = stem_ages(backbone)

subtrees, provenance = {}, {}
for path in snakemake.input.subtrees:
    for row in read_tsv(path):
        subtrees[row["genus"]] = (row["newick"], float(row["crown_age"]))
        provenance[row["genus"]] = {k: row[k] for k in PROVENANCE_COLUMNS}
singles = {r["genus"]: r["unit_id"] for r in read_tsv(snakemake.input.singles)}

tree, status = graft(backbone, subtrees, singles, prune_empty=snakemake.params.prune_empty)
if not is_ultrametric(tree, tolerance=1e-4):
    raise RuntimeError("grafted tree is not ultrametric")

rows = []
for genus, state in sorted(status.items()):
    row = provenance.get(genus, {"genus": genus, "stem_age": f"{stems[genus]:.4f}"})
    row["status"] = state
    if state == "single_unit":
        row["n_units"] = 1
    rows.append(row)

Path(snakemake.output.tree).write_text(
    tree.as_string(schema="newick", suppress_rooting=True, preserve_spaces=True)
)
write_tsv(snakemake.output.provenance, rows, PROVENANCE_COLUMNS)
logging.info("tree with %d tips, root age %.2f Ma", len(tree.leaf_nodes()), root_age(tree))
