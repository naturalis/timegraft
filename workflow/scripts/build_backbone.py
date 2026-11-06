"""Genus topology from Open Tree of Life (static files or API), dated with TimeTree."""

import logging
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

import dendropy

from timegraft.backbone_build import REPORT_COLUMNS, date_topology
from timegraft.seqio import read_tsv, write_tsv
from timegraft.timetree import TimeTree
from timegraft.webapi import JsonCache

settings = snakemake.params.settings
cache = JsonCache(settings.get("cache"))
timetree = TimeTree(cache)
genera = sorted({r["genus"] for r in read_tsv(snakemake.input.reps)})
crowns = []

try:
    if settings.get("source", "static") == "static":
        from timegraft.otol_static import (SynthTree, Taxonomy, crown_ages, genus_topology_newick, lineage,
                                           lookup_ids)

        taxonomy = Taxonomy(snakemake.input.taxonomy, set(genera), context=settings.get("context"))
        synth = SynthTree(snakemake.input.synth)
        newick, placed = genus_topology_newick(synth, taxonomy, genera)
        matches = {g: {**taxonomy.genus[g], "status": placed.get(g, taxonomy.genus[g]["status"])} for g in genera}
        ncbi = {g: lookup_ids(taxonomy, g) for g in genera}
        higher = {g: lineage(taxonomy, g) for g in genera}
        tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True)
        if settings.get("crown_ages_from_timetree", True):
            crowns = crown_ages(taxonomy, timetree, [g for g in genera if g in placed and placed[g] != "unmatched"])
    else:
        from timegraft.opentree import OpenTree, genus_topology

        otol = OpenTree(cache)
        context = settings.get("context") or "All life"
        matches = otol.match_genera(genera, context=context)
        ott_to_genus = {m["ott_id"]: g for g, m in matches.items() if m["status"] == "matched"}
        newick, missing, broken = otol.induced_subtree(ott_to_genus)
        for ott, reason in missing.items():
            matches[ott_to_genus.pop(ott)]["status"] = f"not_in_synthesis:{reason}"
        for ott in broken:
            matches[ott_to_genus[ott]]["status"] = "attached_at_mrca"
        tree = genus_topology(newick, ott_to_genus, broken)
        ncbi = {g: m["ncbi_id"] for g, m in matches.items()}
        higher = None

    tree, report = date_topology(tree, ncbi, timetree, per_side=int(settings.get("per_side", 2)),
                                 root_age=settings.get("root_age"),
                                 per_genus=int(settings.get("species_per_genus", 3)), higher=higher)
finally:
    cache.save()

Path(snakemake.output.tree).write_text(tree.as_string(schema="newick", suppress_rooting=True))
write_tsv(snakemake.output.genera, [{"genus": g, **m} for g, m in sorted(matches.items())],
          ["genus", "status", "ott_id", "ncbi_id"])
write_tsv(snakemake.output.nodes, report, REPORT_COLUMNS)
write_tsv(snakemake.output.crowns, crowns, ["genus", "crown_age", "n_species"])
sources = {}
for row in report:
    sources[row["source"]] = sources.get(row["source"], 0) + 1
logging.info("%d genera in backbone; node ages by source: %s; %d crown ages from TimeTree",
             len(tree.leaf_nodes()), sources, len(crowns))
