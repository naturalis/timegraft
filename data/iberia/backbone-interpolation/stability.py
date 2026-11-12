"""How much the age of an undated backbone node depends on which other genera are sampled (#16).

Usage: python stability.py <run directory with backbone/> <richness.tsv>

Drops a random share of the genera, prunes the backbone, keeps the TimeTree
ages of nodes that survive, interpolates again and compares the age of every
surviving undated node with its age in the full backbone. A method whose ages
are a property of the clade gives the same age whatever else is sampled; one
that depends on the sampled topology does not. Repeated 50 times at each
share.
"""

import csv
import math
import random
import statistics
import sys
from pathlib import Path

import dendropy

from timegraft.backbone_build import interpolate

run, rich_path = Path(sys.argv[1]), Path(sys.argv[2])
tree = dendropy.Tree.get(path=run / "backbone" / "backbone.tre", schema="newick", rooting="force-rooted",
                         preserve_underscores=True)
rows = {r["node"]: r for r in csv.DictReader(open(run / "backbone" / "node_ages.tsv"), delimiter="\t")}
names = {f"n{i}": node for i, node in enumerate(tree.preorder_node_iter())}
rich = list(csv.DictReader(open(rich_path), delimiter="\t"))
dated = {names[n]: float(r["age"]) for n, r in rows.items() if r["source"] in ("timetree", "timetree_higher_taxa")}
dated[tree.seed_node] = float(rows["n0"]["age"])
clade = {node: frozenset(l.taxon.label for l in node.leaf_iter()) for node in tree.internal_nodes()}
by_clade = {v: k for k, v in clade.items()}
METHODS = {
    "even": dict(method="even"),
    "pb_sampled": dict(method="pure_birth"),
    "pb_genera": dict(method="pure_birth", richness={names[r["node"]]: int(r["synth_genera"]) for r in rich if r["synth_genera"]}),
    "pb_species": dict(method="pure_birth", richness={names[r["node"]]: int(r["synth_tips"]) for r in rich if r["synth_tips"]}),
}
full = {m: {**dated, **interpolate(tree, dated, **kw)} for m, kw in METHODS.items()}
genera = sorted(l.taxon.label for l in tree.leaf_node_iter())


def full_node(leaves: frozenset):
    """The node of the full backbone that is the MRCA of ``leaves``."""
    node = tree.find_node_with_taxon_label(next(iter(leaves)))
    while not leaves <= clade.get(node, frozenset()):
        node = node.parent_node
    return node


rng = random.Random(16)
print("share kept, method: median and 90th percentile of |log2(age after subsampling / age in full backbone)|"
      " over undated nodes")
for keep in (0.75, 0.5, 0.25):
    errs = {m: [] for m in METHODS}
    for rep in range(50):
        kept = set(rng.sample(genera, int(len(genera) * keep)))
        sub = dendropy.Tree(tree)
        sub.prune_taxa_with_labels([g for g in genera if g not in kept])
        sub.suppress_unifurcations()
        mapping = {node: full_node(frozenset(l.taxon.label for l in node.leaf_iter())) for node in sub.internal_nodes()}
        anchors = {node: dated[f] for node, f in mapping.items() if f in dated}
        anchors[sub.seed_node] = dated[tree.seed_node] if mapping[sub.seed_node] is tree.seed_node else anchors.get(
            sub.seed_node, full["even"][mapping[sub.seed_node]])
        for m, kw in METHODS.items():
            kw = dict(kw)
            if "richness" in kw:
                kw["richness"] = {node: kw["richness"][f] for node, f in mapping.items() if f in kw["richness"]}
            ages = interpolate(sub, anchors, **kw)
            errs[m].extend(abs(math.log2(a / full[m][mapping[node]])) for node, a in ages.items())
    for m, e in errs.items():
        print(f"  {keep:.2f} {m:>10}: median {statistics.median(e):.2f}, 90th percentile {statistics.quantiles(e, n=10)[-1]:.2f}"
              f" (n={len(e)})")
