"""Tips of the synthetic tree below every backbone node of a run (#16).

Usage: python richness.py <run directory with backbone/> <otol_dir> <out.tsv> [context]

``otol_dir`` holds opentree16.1/labelled_supertree.tre and ott3.7.3/taxonomy.tsv,
as downloaded by the workflow. Checks that the topology induced again from
these files is that of the run's backbone, then writes per node (as named in
node_ages.tsv) the genera in the backbone, and the genera and tips below it
in the synthetic tree.
"""

import csv
import sys
from pathlib import Path

import dendropy

from timegraft.otol_static import SynthTree, Taxonomy, clade_richness, genus_anchors, genus_topology_newick

run, otol, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
context = sys.argv[4] if len(sys.argv) > 4 else "Metazoa"
tree = dendropy.Tree.get(path=run / "backbone" / "backbone.tre", schema="newick", rooting="force-rooted",
                         preserve_underscores=True)
genera = sorted(r["genus"] for r in csv.DictReader(open(run / "backbone" / "genera.tsv"), delimiter="\t"))
taxonomy = Taxonomy(otol / "ott3.7.3" / "taxonomy.tsv", set(genera), context=context)
synth = SynthTree(otol / "opentree16.1" / "labelled_supertree.tre")
anchor, _ = genus_anchors(synth, taxonomy, genera)
newick, _ = genus_topology_newick(synth, taxonomy, genera)
again = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", preserve_underscores=True,
                          taxon_namespace=tree.taxon_namespace)


def clades(t):
    return {frozenset(l.taxon.label for l in n.leaf_iter()) for n in t.internal_nodes()}


assert clades(again) == clades(tree), "the synthesis files do not give the run's topology"
tips = clade_richness(synth, anchor, tree)
genera_rich = clade_richness(synth, anchor, tree, taxonomy)
names = {node: f"n{i}" for i, node in enumerate(tree.preorder_node_iter())}
with open(out, "w") as fh:
    fh.write("node\tn_genera\tsynth_genera\tsynth_tips\n")
    for node in tree.preorder_node_iter():
        if not node.is_leaf():
            fh.write(f"{names[node]}\t{len(node.leaf_nodes())}\t{genera_rich.get(node, '')}\t{tips.get(node, '')}\n")
