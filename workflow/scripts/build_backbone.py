"""Fetch the genus topology from OToL, date it with TimeTree, and write the backbone."""

import logging
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

from timegraft.backbone_build import REPORT_COLUMNS, date_topology
from timegraft.opentree import OpenTree, genus_topology
from timegraft.seqio import read_tsv, write_tsv
from timegraft.timetree import TimeTree
from timegraft.webapi import JsonCache

settings = snakemake.params.settings
cache = JsonCache(settings.get("cache"))
genera = sorted({r["genus"] for r in read_tsv(snakemake.input.reps)})

try:
    otol = OpenTree(cache)
    matches = otol.match_genera(genera, context=settings.get("context", "All life"))
    ott_to_genus = {m["ott_id"]: g for g, m in matches.items() if m["status"] == "matched"}
    newick, missing, broken = otol.induced_subtree(ott_to_genus)
    for ott, reason in missing.items():
        matches[ott_to_genus.pop(ott)]["status"] = f"not_in_synthesis:{reason}"
    for ott in broken:
        matches[ott_to_genus[ott]]["status"] = "broken_attached_at_mrca"
    tree = genus_topology(newick, ott_to_genus, broken)

    ncbi = {g: m["ncbi_id"] for g, m in matches.items()}
    tree, report = date_topology(
        tree,
        ncbi,
        TimeTree(cache),
        per_side=int(settings.get("per_side", 2)),
        root_age=settings.get("root_age"),
    )
finally:
    cache.save()

Path(snakemake.output.tree).write_text(tree.as_string(schema="newick", suppress_rooting=True))
write_tsv(
    snakemake.output.genera,
    [{"genus": g, **m} for g, m in sorted(matches.items())],
    ["genus", "status", "ott_id", "ncbi_id"],
)
write_tsv(snakemake.output.nodes, report, REPORT_COLUMNS)
sources = {}
for row in report:
    sources[row["source"]] = sources.get(row["source"], 0) + 1
logging.info("%d genera in backbone; node ages by source: %s", len(tree.leaf_nodes()), sources)
