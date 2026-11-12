"""Genus topology from Open Tree of Life (static files or API), dated with TimeTree."""

import json
import logging
from collections import Counter, defaultdict
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

import dendropy

from timegraft import __version__
from timegraft.backbone_build import (GENERA_COLUMNS, REPORT_COLUMNS, date_topology, genera_rows, join_ncbi,
                                      records_lineage, timetree_summary)
from timegraft.col import SPECIES_COLUMNS, compare_ids
from timegraft.otol_static import CROWN_COLUMNS
from timegraft.seqio import read_tsv, write_tsv
from timegraft.timetree import CACHE_VERSION, TimeTree, is_legacy
from timegraft.webapi import JsonCache

settings = snakemake.params.settings
cache = JsonCache(settings.get("cache"))
# version 1 cache entries are asked again for their evidence, keeping their ages (#47)
timetree = TimeTree(cache, upgrade_legacy=bool(settings.get("upgrade_cache", True)))
genera = sorted({r["genus"] for r in read_tsv(snakemake.input.reps)})
crowns = []
lineage_of = records_lineage(read_tsv(snakemake.input.records)) if snakemake.input.records else {}
extra = [lineage_of]
synth_nodes, support = None, None
concordance = []

# COL XR as identifier hub (#48): genera (and optionally species) resolved with the lineage of their records
col_cfg = settings.get("col") or {}
col_rows, col_species, clb = {}, [], None
if col_cfg.get("enabled", True):
    from timegraft.col import KINGDOMS, ChecklistBank, resolve_genera

    clb = ChecklistBank(col_cfg.get("release", "3LXR"), cache)
    kingdom = KINGDOMS.get(settings.get("context") or "", settings.get("context"))
    try:
        col_rows = resolve_genera(clb, genera, lineage_of, kingdom, workers=int(col_cfg.get("workers", 4)))
    finally:
        cache.save()
    extra.append(col_rows)

try:
    if settings.get("source", "static") == "static":
        from timegraft.otol_static import (SynthTree, Taxonomy, clade_richness, crown_ages,
                                           genus_anchor_info, genus_anchors, genus_topology_newick, lineage,
                                           lineage_names, lookup_ids, synth_nodes as synth_node_index,
                                           synthesis_support)

        # COL breaks ties between OTT homonyms, and only those (#48)
        prefer = {g: int(r["col_ott_id"]) for g, r in col_rows.items()
                  if r.get("col_ott_id") and r.get("col_confidence") in ("high", "medium")}
        taxonomy = Taxonomy(snakemake.input.taxonomy, set(genera), context=settings.get("context"), prefer=prefer)
        synth = SynthTree(snakemake.input.synth)
        newick, placed = genus_topology_newick(synth, taxonomy, genera)
        matches = {g: {**taxonomy.genus[g], "status": placed.get(g, taxonomy.genus[g]["status"])} for g in genera}
        ncbi = {g: lookup_ids(taxonomy, g) for g in genera}
        extra.append(join_ncbi(ncbi, matches, col_rows, taxonomy.match_info, col_cfg.get("join", "col")))
        extra.append({g: {"col_ott_vs_ott": compare_ids(r.get("col_ott_id", ""), "", str(matches[g]["ott_id"] or ""))}
                      for g, r in col_rows.items()})
        if col_cfg.get("species", False) and clb is not None:
            from timegraft.col import resolve_species

            names = {}
            for r in read_tsv(snakemake.input.records):
                if r.get("species") and r.get("genus") in matches:
                    names.setdefault(r["species"].strip(), r["genus"])
            try:
                col_species = resolve_species(clb, names, lineage_of, kingdom, taxonomy.species_ncbi,
                                              workers=int(col_cfg.get("workers", 4)))
            finally:
                cache.save()
        higher = {g: lineage(taxonomy, g) for g in genera}
        tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True)
        # genera in the full clade of every node, for the pure-birth interpolation (#16)
        anchors = genus_anchors(synth, taxonomy, genera)[0]
        richness = clade_richness(synth, anchors, tree, taxonomy)
        # evidence for the consistency layer (#47): how each genus was matched and anchored,
        # and which synthesis node stands behind each backbone node
        extra += [taxonomy.match_info, {g: lineage_names(taxonomy, g) for g in genera},
                  genus_anchor_info(synth, taxonomy, anchors)]
        synth_nodes = {n: synth.label[m] for n, m in synth_node_index(synth, anchors, tree).items()
                       if not n.is_leaf() and synth.label[m]}
        if snakemake.input.annotations:
            support = synthesis_support(snakemake.input.annotations, set(synth_nodes.values()))
        # placement concordance, report only (#49)
        conc_cfg = settings.get("concordance") or {}
        if conc_cfg.get("enabled", True):
            from timegraft import concordance as cc
            from timegraft.col import NCBI_DATASET, ChecklistBank

            info = taxonomy.match_info
            reference, phylum = {}, {}
            for g in set(genera) | set(lineage_of):
                col = col_rows.get(g) or {}
                rec = lineage_of.get(g) or {}
                if col.get("col_id"):
                    reference[g] = ("col", {r: col.get(f"col_{r}", "") for r in ("family", "order", "class")})
                elif rec:
                    reference[g] = ("records", {r: rec.get(f"records_{r}", "") for r in ("family", "order", "class")})
                phylum[g] = col.get("col_phylum") or rec.get("records_phylum") or ""
            names = cc.name_index(taxonomy)
            client = clb or ChecklistBank(col_cfg.get("release", 316691), cache)
            merged = {g: [i for i in (info.get(g, {}).get("ott_ncbi_ids") or "").split(";") if i] for g in genera}
            try:
                classifications = {i: client.classification(NCBI_DATASET, i)
                                   for ids in merged.values() if len(ids) > 1 for i in ids}
            finally:
                cache.save()
            hits = defaultdict(list)
            for r in read_tsv(snakemake.input.screened) if snakemake.input.screened else []:
                if r.get("best_other_genus"):
                    hits[r["genus"]].append(r["best_other_genus"])
            sisters = cc.backbone_sisters(tree)
            flag_level = cc.LEVELS[conc_cfg.get("flag_level", "class")]
            for g in genera:
                match = matches[g]
                row = {c: "" for c in cc.COLUMNS}
                row.update({"genus": g, "status": match["status"], "ott_id": match.get("ott_id") or "",
                            **{k: v for k, v in lineage_names(taxonomy, g).items() if k in cc.COLUMNS}})
                source, ref = reference.get(g, ("", {}))
                row["reference_source"] = source
                if match.get("ott_id"):
                    row.update(cc.lineage_distance(taxonomy, names, match["ott_id"], ref))
                    row["merged_ncbi_ids"] = ";".join(merged[g]) if len(merged[g]) > 1 else ""
                    row.update(cc.merged_id_distance(classifications, merged[g]))
                if g in sisters:
                    row.update(cc.neighbour_concordance(g, hits.get(g, []), sisters[g], phylum,
                                                        min_hits=int(conc_cfg.get("min_hits", 5))))
                concordance.append(cc.aggregate(row, flag_level, float(conc_cfg.get("max_share", 0.0)),
                                                int(conc_cfg.get("min_flags", 2))))
            concordance = cc.rank_rows(concordance)
            extra.append({r["genus"]: {k: r[k] for k in ("lineage_distance", "merged_id_distance",
                                                         "neighbour_concordance", "concordance_flag")}
                          for r in concordance})
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
        extra.append(join_ncbi(ncbi, matches, col_rows, None, col_cfg.get("join", "col")))
        higher = None
        richness = None

    tree, report = date_topology(tree, ncbi, timetree, per_side=int(settings.get("per_side", 2)),
                                 root_age=settings.get("root_age"),
                                 per_genus=int(settings.get("species_per_genus", 3)), higher=higher,
                                 interpolation=settings.get("interpolation", "pure_birth"), richness=richness,
                                 synth_nodes=synth_nodes, support=support)
finally:
    cache.save()

Path(snakemake.output.tree).write_text(tree.as_string(schema="newick", suppress_rooting=True))
write_tsv(snakemake.output.genera, genera_rows(matches, *extra), GENERA_COLUMNS)
write_tsv(snakemake.output.nodes, report, REPORT_COLUMNS)
write_tsv(snakemake.output.crowns, crowns, CROWN_COLUMNS)
write_tsv(snakemake.output.col_species, col_species, SPECIES_COLUMNS)
from timegraft.concordance import COLUMNS as CONCORDANCE_COLUMNS
write_tsv(snakemake.output.concordance, concordance, CONCORDANCE_COLUMNS)

# what the backbone was built from, and in what format the outputs are (#47)
entries = [v for k, v in cache.data.items() if k.startswith("mrca:") and isinstance(v, dict)]
metadata = {
    "output_format": {"genera.tsv": 4, "concordance.tsv": 1, "node_ages.tsv": 2, "crown_ages.tsv": 2, "timetree_cache": CACHE_VERSION},
    "timegraft_version": __version__,
    "source": settings.get("source", "static"),
    "taxonomy": settings.get("taxonomy", "ott3.7.3") if settings.get("source", "static") == "static" else "api",
    "synthesis": settings.get("synthesis", "opentree16.1") if settings.get("source", "static") == "static" else "api",
    "context": settings.get("context"),
    "synthesis_support": support is not None,
    "concordance": {
        "settings": settings.get("concordance") or {},
        "genera_by_flag": dict(Counter(r["concordance_flag"] or "concordant" for r in concordance)),
    },
    "col": None if clb is None else {
        **clb.info,
        "join": col_cfg.get("join", "col"),
        "genera_by_confidence": dict(Counter(r.get("col_confidence") or "unresolved" for r in col_rows.values())),
        "genera_by_ncbi_path": dict(Counter(r.get("col_ncbi_path") or "none" for r in col_rows.values())),
        "genera_joined_with_col_ncbi": sum(1 for r in genera_rows(matches, *extra) if r["join_ncbi_source"] == "col"),
        "species": len(col_species),
    },
    "timetree": {
        **timetree_summary(report, crowns),
        "cache_entries": len(entries),
        "cache_entries_legacy": sum(1 for v in entries if is_legacy(v)),
        "cache_entries_upgraded": sum(1 for v in entries if "upgraded_from" in v),
        "cache_entries_fresh_age_differs": sum(1 for v in entries if "fresh_age" in v),
    },
}
Path(snakemake.output.metadata).write_text(json.dumps(metadata, indent=1, sort_keys=True) + "\n")
sources = {}
for row in report:
    sources[row["source"]] = sources.get(row["source"], 0) + 1
logging.info("%d genera in backbone; node ages by source: %s (interpolation: %s); %d crown ages from TimeTree",
             len(tree.leaf_nodes()), sources, settings.get("interpolation", "pure_birth"), len(crowns))
