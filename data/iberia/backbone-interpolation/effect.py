"""Effect of the interpolation method on the Iberian backbone as it was built (#16).

Usage: python effect.py <run directory with backbone/> <richness.tsv>

Keeps the run's TimeTree ages, interpolates the other nodes with each method
and reports how undated nodes and genus stem ages change against even spacing
(which reproduces the run), the total branch length, and how often a genus
stem comes out younger than the genus crown age from TimeTree, which cannot be.
"""

import csv
import math
import statistics
import sys
from pathlib import Path

import dendropy

from timegraft.backbone_build import interpolate

run, rich_path = Path(sys.argv[1]), Path(sys.argv[2])
tree = dendropy.Tree.get(path=run / "backbone" / "backbone.tre", schema="newick", rooting="force-rooted",
                         preserve_underscores=True)
rows = {r["node"]: r for r in csv.DictReader(open(run / "backbone" / "node_ages.tsv"), delimiter="\t")}
names = {f"n{i}": n for i, n in enumerate(tree.preorder_node_iter())}
dated = {names[k]: float(r["age"]) for k, r in rows.items() if r["source"] in ("timetree", "timetree_higher_taxa")}
dated[tree.seed_node] = float(rows["n0"]["age"])
rich = list(csv.DictReader(open(rich_path), delimiter="\t"))
METHODS = {
    "even": dict(method="even"),
    "pb_sampled": dict(method="pure_birth"),
    "pb_genera": dict(method="pure_birth", richness={names[r["node"]]: int(r["synth_genera"]) for r in rich if r["synth_genera"]}),
    "pb_species": dict(method="pure_birth", richness={names[r["node"]]: int(r["synth_tips"]) for r in rich if r["synth_tips"]}),
}
ages = {m: {**dated, **interpolate(tree, dated, **kw)} for m, kw in METHODS.items()}
for node in tree.internal_nodes():
    assert abs(ages["even"][node] - float(rows[[k for k, v in names.items() if v is node][0]]["age"])) < 1e-3
undated = [n for n in tree.internal_nodes() if n not in dated]
leaves = {l.taxon.label: l for l in tree.leaf_node_iter()}
crowns = {r["genus"]: float(r["crown_age"]) for r in csv.DictReader(open(run / "backbone" / "crown_ages.tsv"), delimiter="\t")
          if r["crown_age"]}
interpolated_stems = [g for g in leaves if leaves[g].parent_node not in dated]
with_crown = [g for g in interpolated_stems if g in crowns]
print(f"{len(undated)} undated nodes; {len(interpolated_stems)} genera with an interpolated stem, "
      f"{len(with_crown)} of them with a TimeTree crown age")


def length(m):
    return sum(ages[m][n.parent_node] - ages[m].get(n, 0.0) for n in tree.preorder_node_iter() if n.parent_node)


for m in METHODS:
    shift = [math.log2(ages[m][n] / ages["even"][n]) for n in undated]
    q = statistics.quantiles(shift, n=4)
    stems = [ages[m][leaves[g].parent_node] for g in interpolated_stems]
    conflicts = sum(ages[m][leaves[g].parent_node] < crowns[g] for g in with_crown)
    print(f"  {m:>10}: undated nodes log2 vs even median {q[1]:+.2f} (quartiles {q[0]:+.2f}, {q[2]:+.2f}); "
          f"median interpolated stem {statistics.median(stems):.1f} Ma; total length {length(m):,.0f} Ma; "
          f"stem younger than TimeTree crown: {conflicts}")
