"""Rebuild the rooted, undated genus trees of the Iberian run (same code path as build_genus)
and write what the chronos experiments need.

Usage: export_trees.py RUN_DIR RECORDS_GZ OUT_DIR
"""
import csv
import gzip
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import dendropy
import numpy as np

from timegraft import distances
from timegraft.bionj import bionj
from timegraft.build import prepare_sequences, tree_to_newick
from timegraft.dating import mean_root_to_tip
from timegraft.rooting import root_genus

csv.field_size_limit(10**9)
run, records_gz, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)


def one(args):
    genus, ingroup, outgroup, stem = args
    items = prepare_sequences(ingroup + outgroup, "k2p", "mafft")
    labels = [n for n, _ in items]
    dist = distances.matrix([s for _, s in items], method="k2p")
    unrooted = bionj(labels, dist)
    in_names, out_names = [n for n, _ in ingroup], [n for n, _ in outgroup]
    tree, how, det = root_genus(unrooted, in_names, out_names, rule="tolerant", fallback="mad")
    neg = sum(1 for e in tree.postorder_edge_iter() if e.tail_node is not None and (e.length or 0) < 0)
    zero = sum(1 for e in tree.postorder_edge_iter() if e.tail_node is not None and (e.length or 0) == 0)
    idx = {n: i for i, n in enumerate(labels)}
    ii = [idx[n] for n in in_names]
    # mean K2P distance from the ingroup to each outgroup unit, and the largest distance within the genus
    d_out = {o: float(np.mean(dist[np.ix_(ii, [idx[o]])])) for o in out_names}
    d_in_max = float(dist[np.ix_(ii, ii)].max())
    return {
        "genus": genus, "n": len(ingroup), "stem_age": stem, "rooting": how,
        "violations": det["violations"], "mrt": mean_root_to_tip(tree),
        "n_neg": neg, "n_zero": zero, "d_in_max": d_in_max,
        "d_out": ";".join(f"{o}={d:.5f}" for o, d in d_out.items()),
        "newick": tree_to_newick(tree),
    }


if __name__ == "__main__":
    groups = defaultdict(lambda: {"in": [], "out": [], "stem": None})
    for shard in sorted((run / "genera").glob("shard_*.tsv")):
        with open(shard) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                g = groups[row["genus"]]
                g[row["role"]].append((row["unit_id"], row["sequence"]))
                g["stem"] = float(row["stem_age"])
    jobs = [(g, v["in"], v["out"], v["stem"]) for g, v in sorted(groups.items()) if len(v["in"]) >= 2]
    jobs.sort(key=lambda j: -len(j[1]))
    with Pool(2) as pool:
        rows = pool.map(one, jobs, chunksize=1)
    rows.sort(key=lambda r: r["genus"])
    with open(out / "trees.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # unit -> genus (as assigned in the run) and genus -> higher taxa (majority over records)
    unit_genus = {}
    for g, v in groups.items():
        for u, _ in v["in"]:
            unit_genus[u] = g
    with open(run / "genera" / "singles.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            unit_genus[row["unit_id"]] = row["genus"]
    ranks = defaultdict(lambda: defaultdict(Counter))
    with gzip.open(records_gz, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            for r in ("family", "order", "class", "phylum"):
                if row[r]:
                    ranks[row["genus"]][r][row[r]] += 1
    with open(out / "taxonomy.tsv", "w") as fh:
        fh.write("genus\tfamily\torder\tclass\tphylum\n")
        for g in sorted(ranks):
            fh.write(g + "\t" + "\t".join(
                ranks[g][r].most_common(1)[0][0] if ranks[g][r] else "" for r in ("family", "order", "class", "phylum")) + "\n")

    # backbone divergence age between each genus and the genus of each of its outgroup units
    bb = dendropy.Tree.get(path=str(run / "backbone" / "backbone.tre"), schema="newick", preserve_underscores=True)
    bb.calc_node_ages(ultrametricity_precision=1e-3)
    tip = {l.taxon.label: l for l in bb.leaf_node_iter()}
    anc = {}
    for label, leaf in tip.items():
        chain, node = [], leaf
        while node is not None:
            chain.append(node)
            node = node.parent_node
        anc[label] = chain
    with open(out / "outgroup_ages.tsv", "w") as fh:
        fh.write("genus\tunit\tout_genus\tdiv_age\tk2p\n")
        for r in rows:
            for item in filter(None, r["d_out"].split(";")):
                u, d = item.rsplit("=", 1)
                og = unit_genus.get(u, "")
                age = ""
                if og in anc and r["genus"] in anc:
                    mine = set(map(id, anc[r["genus"]]))
                    age = f"{next(n for n in anc[og] if id(n) in mine).age:.4f}"
                fh.write(f"{r['genus']}\t{u}\t{og}\t{age}\t{d}\n")
    print(len(rows), "genus trees written")
