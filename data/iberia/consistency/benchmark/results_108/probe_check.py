"""Probe unnamed TimeTree answers with one more genus (#108, step 2; prototype of ttcheck).

    python probe_check.py BACKBONE_DIR NODE_CHECK LINEAGES CACHE OUT [PROBES] [named]

For an unnamed answer whose NCBI common ancestor also holds genera outside the backbone node (``wider`` in
NODE_CHECK, from node_check.py), TimeTree's dated node may or may not contain them. Every answer carries TimeTree's
own id of the node it dated (``mrca_ttid``). Asking again with the NCBI id of one such genus (inside the common
ancestor, outside the backbone node) added tells which:
the same ``mrca_ttid`` means the genus is inside TimeTree's node, so the answer dates a wider node than intended;
another means it is outside. Up to PROBES (default 3) genera are tried, nearest to the backbone node first (those
of its sister clade, then of the sister of its parent, ...). CACHE is a TimeTree cache (json), extended with the
probe queries; LINEAGES is the NCBI lineage cache of ncbi_refinement.py and node_check.py.

With ``named``, the named ``wider`` answers are probed instead. Their content is known (TimeTree's node is the
NCBI common ancestor), so every probe should come back inside: a check of the method.
"""

import csv
import json
import sys

import dendropy

from timegraft.webapi import JsonCache
from timegraft.timetree import TimeTree

backbone, node_check, lineage_path, cache_path, out = sys.argv[1:6]
n_probes = int(sys.argv[6]) if len(sys.argv) > 6 else 3
named = len(sys.argv) > 7 and sys.argv[7] == "named"
lineages = json.load(open(lineage_path))
tree = dendropy.Tree.get(path=f"{backbone}/backbone.tre", schema="newick", preserve_underscores=True,
                         rooting="force-rooted")
nodes = {f"n{i}": node for i, node in enumerate(tree.preorder_node_iter())}
ncbi = {r["genus"]: int(r["join_ncbi_id"]) for r in csv.DictReader(open(f"{backbone}/genera.tsv"), delimiter="\t")
        if r.get("join_ncbi_id")}
query = {r["node"]: [int(i) for i in r["query_ids"].split("+")]
         for r in csv.DictReader(open(f"{backbone}/node_ages.tsv"), delimiter="\t") if r["query_ids"]}
timetree = TimeTree(JsonCache(cache_path))


def nearest_outside(node, allowed: set[str]) -> list[str]:
    """Genera in ``allowed`` outside ``node``, nearest first: the sister clades, then those of the parent, ..."""
    inside = {leaf.taxon.label for leaf in node.leaf_iter()}
    order, current = [], node
    while current.parent_node is not None and len(order) < 50:
        for sister in current.parent_node.child_node_iter():
            if sister is current:
                continue
            order += sorted(leaf.taxon.label for leaf in sister.leaf_iter()
                            if leaf.taxon.label in allowed and leaf.taxon.label not in inside)
        current = current.parent_node
    return order


rows = list(csv.DictReader(open(node_check), delimiter="\t"))
genus_lineage = {g: set((lineages.get(str(i)) or []) + [str(i)]) for g, i in ncbi.items() if lineages.get(str(i))}
with open(out, "w") as h:
    w = csv.writer(h, delimiter="\t", lineterminator="\n")
    w.writerow(["node", "n_genera", "example_genera", "ttid", "probe_result", "probes", "probe_ttids"])
    for r in rows:
        if r["class"] != "wider" or (r["timetree_mrca_ncbi"] != "0") != named:
            continue
        ids = query[r["node"]]
        first = timetree._query(ids)
        ttid = first.get("mrca_ttid")
        found = set(first.get("found") or [])
        node = nodes[r["node"]]
        # candidate probes: genera outside the node with an NCBI id, nearest first
        in_lca = {g for g, path in genus_lineage.items() if r["ncbi_common_ancestor"] in path}
        candidates = [g for g in nearest_outside(node, in_lca) if ncbi[g] not in found]
        result, tried, ttids = "inconclusive", [], []
        for g in candidates:
            entry = timetree._query(sorted(set(ids) | {ncbi[g]}))
            if ncbi[g] not in (entry.get("found") or []):
                continue  # TimeTree does not know the probe: try the next
            tried.append(g)
            ttids.append(str(entry.get("mrca_ttid")))
            if entry.get("mrca_ttid") == ttid:
                result = "wider"
                break
            if len(tried) >= n_probes:
                break
        if result != "wider" and tried:
            result = "outside"
        w.writerow([r["node"], r["n_genera"], r["example_genera"], ttid, result, ";".join(tried), ";".join(ttids)])
timetree.cache.save()
