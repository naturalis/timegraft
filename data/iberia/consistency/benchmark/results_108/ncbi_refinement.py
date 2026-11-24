"""Does TimeTree's tree refine the NCBI taxonomy? (#108, step 2.1)

    python ncbi_refinement.py NODE_AGES CACHE LINEAGES OUT_NODES OUT_CACHE

If TimeTree's global tree only resolves polytomies of the NCBI taxonomy (every NCBI taxon a clade in it), the
node TimeTree reports for an ``mrca`` query, when it is named, must be the NCBI common ancestor of the ids it
recognised: a named node is an NCBI taxon, it contains those ids, and no smaller NCBI taxon does. This script
checks that for every named answer, both those calibrating the backbone (NODE_AGES, a node_ages.tsv) and all
in the TimeTree cache (CACHE, json or json.gz).

NCBI lineages come from ChecklistBank's copy of the NCBI taxonomy (dataset 2169), kept in LINEAGES (a json file,
created if missing, so that a rerun asks nothing again). Each answer is classed as:

- ``same``: TimeTree's node is the NCBI common ancestor of the recognised ids;
- ``deeper``: TimeTree's node is an NCBI ancestor of that common ancestor;
- ``shallower``: TimeTree's node lies inside the common ancestor, so in TimeTree's tree it contains ids that
  NCBI places outside it;
- ``elsewhere``: neither contains the other;
- ``unknown``: a lineage is missing in ChecklistBank.
"""

import csv
import gzip
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from timegraft.col import NCBI_DATASET, request_json

node_ages, cache_path, lineage_path, out_nodes, out_cache = sys.argv[1:6]
URL = "https://api.checklistbank.org/dataset/{}/taxon/{}/classification"


def load(path):
    path = Path(path)
    return json.load(gzip.open(path) if path.suffix == ".gz" else open(path))


lineages = json.load(open(lineage_path)) if Path(lineage_path).exists() else {}


def lineage(taxon: str) -> list[str] | None:
    """NCBI ids from the root to the taxon itself, or None where ChecklistBank does not know it."""
    if taxon not in lineages:
        try:
            value = request_json(URL.format(NCBI_DATASET, taxon))
        except Exception:  # noqa: BLE001 - an unknown id is answered with an error
            value = None
        lineages[taxon] = [str(c["id"]) for c in value if c.get("id")] if isinstance(value, list) and value else None
    known = lineages[taxon]
    return None if known is None else known + [taxon]


def classify(found: list[str], mrca: str) -> tuple[str, str]:
    paths = [lineage(i) for i in found]
    target = lineage(mrca)
    if target is None or any(p is None for p in paths):
        return "unknown", ""
    common = []
    for step in zip(*paths):
        if len(set(step)) > 1:
            break
        common.append(step[0])
    lca = common[-1] if common else ""
    if lca == mrca:
        return "same", lca
    if mrca in common:
        return "deeper", lca
    if lca and lca in target:
        return "shallower", lca
    return "elsewhere", lca


rows = [r for r in csv.DictReader(open(node_ages), delimiter="\t")
        if r["timetree_age"] and r["timetree_mrca_ncbi"] not in ("", "0") and "+" in r["timetree_found_ids"]]
cache = {k: v for k, v in load(cache_path).items()
         if k.startswith("mrca:") and isinstance(v, dict) and v.get("age") and v.get("mrca_ncbi")
         and len(v.get("found") or []) >= 2}
todo = {i for r in rows for i in r["timetree_found_ids"].split("+") + [r["timetree_mrca_ncbi"]]}
todo |= {str(i) for v in cache.values() for i in list(v["found"]) + [v["mrca_ncbi"]]}
todo -= set(lineages)
with ThreadPoolExecutor(max_workers=8) as pool:
    list(pool.map(lineage, sorted(todo)))
json.dump(lineages, open(lineage_path, "w"))

with open(out_nodes, "w") as h:
    w = csv.writer(h, delimiter="\t", lineterminator="\n")
    w.writerow(["node", "example_genera", "found_ids", "timetree_mrca_ncbi", "ncbi_common_ancestor", "class",
                "timetree_age", "timetree_n_studies"])
    for r in rows:
        cls, lca = classify(r["timetree_found_ids"].split("+"), r["timetree_mrca_ncbi"])
        w.writerow([r["node"], r["example_genera"], r["timetree_found_ids"], r["timetree_mrca_ncbi"], lca, cls,
                    r["timetree_age"], r["timetree_n_studies"]])
with open(out_cache, "w") as h:
    w = csv.writer(h, delimiter="\t", lineterminator="\n")
    w.writerow(["query", "found_ids", "timetree_mrca_ncbi", "ncbi_common_ancestor", "class", "n_studies"])
    for key, v in sorted(cache.items()):
        found = [str(i) for i in v["found"]]
        cls, lca = classify(found, str(v["mrca_ncbi"]))
        w.writerow([key, "+".join(found), v["mrca_ncbi"], lca, cls, len(v.get("studies") or [])])
