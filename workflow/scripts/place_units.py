"""Place reference units without a usable genus on the grafted tree."""

import logging
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

import dendropy

from collections import defaultdict

from timegraft.placement import attach_to_genus_stem, classify, place_all, vsearch_hits
from timegraft.seqio import read_tsv, write_tsv

COLUMNS = ["unit_id", "reason", "status", "nearest_units", "genus", "identity", "attach_age"]

tree = dendropy.Tree.get(path=snakemake.input.tree, schema="newick", preserve_underscores=True,
                         rooting="force-rooted")
unassigned = read_tsv(snakemake.input.unassigned)
# units with a long branch that the genus step left out of their subtree (#55)
for path in getattr(snakemake.input, "set_aside", None) or []:
    unassigned += read_tsv(path)
rows_out = []
if snakemake.params.enabled and unassigned:
    cfg = snakemake.params.settings
    in_tree = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    records = read_tsv(snakemake.input.records)
    record_unit = {r["record_id"]: r["unit_id"] for r in records}
    unit_genus = {r["unit_id"]: r["genus"] for r in read_tsv(snakemake.input.reps)}
    refs = [(r["record_id"], r["sequence"]) for r in records if r["unit_id"] in in_tree]
    queries = [(u["unit_id"], u["sequence"]) for u in unassigned]
    hits = vsearch_hits(queries, refs, min_identity=cfg["place_identity"],
                        query_cov=cfg.get("query_cov", 0.9), maxaccepts=cfg.get("maxaccepts", 20),
                        threads=snakemake.threads)
    # a different unit is never a match, so everything above place_identity is placed
    rows = classify(hits, [u for u, _ in queries], record_unit, unit_genus, in_tree, match_identity=2.0)
    placed = place_all(tree, rows, height=cfg.get("height", "rate"), rate=snakemake.params.rate,
                       terminal_fraction=cfg.get("terminal_fraction", 0.5))
    reason = {u["unit_id"]: u["reason"] for u in unassigned}
    # units set aside by screening (#38, #55) that match nothing well enough to be placed go on the stem
    # of their own genus: they stay with it, but at an age the genus's crown and stem set
    own_genus = {u["unit_id"]: u["genera"].split(":")[0] for u in unassigned if u["reason"].startswith("screen_")}
    genus_units = defaultdict(list)
    for unit, genus in unit_genus.items():
        genus_units[genus].append(unit)
    for p in placed:
        unit = p["asv_id"]
        if p["status"] == "unplaced" and unit in own_genus:
            age = attach_to_genus_stem(tree, unit, genus_units[own_genus[unit]], cfg.get("stem_fraction", 0.5))
            if age is not None:
                p.update(status="genus_stem", genus=own_genus[unit], attach_age=f"{age:.4f}")
        rows_out.append({"unit_id": unit, "reason": reason[unit], "status": p["status"],
                         "nearest_units": p["unit_id"], "genus": p["genus"], "identity": p["identity"],
                         "attach_age": p["attach_age"]})
else:
    rows_out = [{"unit_id": u["unit_id"], "reason": u["reason"], "status": "not_attempted"} for u in unassigned]

Path(snakemake.output.tree).write_text(
    tree.as_string(schema="newick", suppress_rooting=True, preserve_spaces=True)
)
write_tsv(snakemake.output.placements, rows_out, COLUMNS)
logging.info("%d unassigned units; placed %d; on the stem of their genus %d", len(unassigned),
             sum(1 for r in rows_out if r["status"] == "placed"), sum(1 for r in rows_out if r["status"] == "genus_stem"))
