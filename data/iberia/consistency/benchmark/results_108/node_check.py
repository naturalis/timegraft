"""Which backbone node did TimeTree date? (#108, step 2)

    python node_check.py BACKBONE_DIR LINEAGES OUT

Step 2.1 (ncbi_refinement.py) found that TimeTree's tree refines the NCBI taxonomy: a named answer is always the
NCBI common ancestor of the ids TimeTree recognised. So for every answer, named or not, the node TimeTree dated
lies within that common ancestor. Its genera in the run (those whose NCBI lineage passes through it) are compared
with the genera of the backbone node:

- ``match``: the same genera. TimeTree dated the intended node, or (unnamed) a node within it that spans it;
- ``wider``: the common ancestor also holds genera outside the backbone node. Named, TimeTree dated a deeper
  node; unnamed, it may have;
- ``narrower``: the common ancestor lacks genera of the backbone node. TimeTree dated a node within it;
- ``unknown``: an NCBI lineage is missing.

Genera without an NCBI id, or whose lineage ChecklistBank lacks, are left out of both sides. LINEAGES is the
lineage cache of ncbi_refinement.py; lineages missing from it are fetched and added.
"""

import csv
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import dendropy

from timegraft.col import NCBI_DATASET, request_json

backbone, lineage_path, out = sys.argv[1:4]
URL = "https://api.checklistbank.org/dataset/{}/taxon/{}/classification"
lineages = json.load(open(lineage_path)) if Path(lineage_path).exists() else {}


def lineage(taxon: str) -> list[str] | None:
    if taxon not in lineages:
        try:
            value = request_json(URL.format(NCBI_DATASET, taxon))
        except Exception:  # noqa: BLE001
            value = None
        lineages[taxon] = [str(c["id"]) for c in value if c.get("id")] if isinstance(value, list) and value else None
    known = lineages[taxon]
    return None if known is None else known + [taxon]


tree = dendropy.Tree.get(path=f"{backbone}/backbone.tre", schema="newick", preserve_underscores=True,
                         rooting="force-rooted")
genera_of = {f"n{i}": {leaf.taxon.label for leaf in node.leaf_iter()}
             for i, node in enumerate(tree.preorder_node_iter())}
ncbi = {r["genus"]: r["join_ncbi_id"] for r in csv.DictReader(open(f"{backbone}/genera.tsv"), delimiter="\t")
        if r.get("join_ncbi_id")}
rows = [r for r in csv.DictReader(open(f"{backbone}/node_ages.tsv"), delimiter="\t") if r["timetree_age"]]
todo = set(ncbi.values()) | {i for r in rows for i in r["timetree_found_ids"].split("+") if i}
with ThreadPoolExecutor(max_workers=8) as pool:
    list(pool.map(lineage, sorted(todo - set(lineages))))
json.dump(lineages, open(lineage_path, "w"))
genus_lineage = {g: set(lineage(i) or ()) for g, i in ncbi.items()}
known = {g for g, path in genus_lineage.items() if path}

with open(out, "w") as h:
    w = csv.writer(h, delimiter="\t", lineterminator="\n")
    w.writerow(["node", "n_genera", "example_genera", "timetree_mrca_ncbi", "ncbi_common_ancestor", "class",
                "n_inside", "n_outside", "n_missing", "outside_examples", "missing_examples", "timetree_age",
                "calibration_age", "age", "timetree_n_studies", "timetree_attempt"])
    for r in rows:
        paths = [lineage(i) for i in r["timetree_found_ids"].split("+") if i]
        node_genera = genera_of[r["node"]] & known
        if not paths or any(p is None for p in paths):
            cls, lca, inside, outside, missing = "unknown", "", set(), set(), set()
        else:
            prefix = []
            for step in zip(*paths):
                if len(set(step)) > 1:
                    break
                prefix.append(step[0])
            lca = prefix[-1] if prefix else ""
            in_lca = {g for g in known if lca in genus_lineage[g]}
            inside, outside, missing = in_lca & node_genera, in_lca - node_genera, node_genera - in_lca
            cls = "wider" if outside else ("narrower" if missing else "match")
        w.writerow([r["node"], r["n_genera"], r["example_genera"], r["timetree_mrca_ncbi"], lca, cls, len(inside),
                    len(outside), len(missing), ";".join(sorted(outside)[:5]), ";".join(sorted(missing)[:5]),
                    r["timetree_age"], r.get("calibration_age", ""), r["age"], r["timetree_n_studies"],
                    r["timetree_attempt"]])
