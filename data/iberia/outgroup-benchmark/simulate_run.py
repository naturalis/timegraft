"""Simulated timegraft run with known genus roots, for testing the rooting diagnostics."""
import json
import random
import sys
from pathlib import Path

import dendropy
import numpy as np
import yaml

from timegraft.backbone import nearest_genera, stem_ages
from timegraft.build import build_genus
from timegraft.graft import PROVENANCE_COLUMNS
from timegraft.seqio import write_tsv
from timegraft.build import tree_to_newick

out = Path(sys.argv[1])
rng = random.Random(7)
nrng = np.random.default_rng(7)
N_GENERA, L, RATE, KAPPA, ALPHA = 60, 658, 0.0115, 10.0, 0.3

# backbone: pure birth, scaled to 150 Ma
bb = dendropy.simulate.treesim.birth_death_tree(1.0, 0.0, num_extant_tips=N_GENERA, rng=rng)
for i, leaf in enumerate(bb.leaf_node_iter()):
    leaf.taxon.label = f"G{i:02d}"
bb.calc_node_ages()
scale = 150.0 / bb.seed_node.age
for e in bb.postorder_edge_iter():
    if e.length:
        e.length *= scale
bb.calc_node_ages()

# genus subtrees grafted on the tips
truth = {}
full = dendropy.Tree.get(data=bb.as_string(schema="newick"), schema="newick", rooting="force-rooted")
for leaf in list(full.leaf_node_iter()):
    g = leaf.taxon.label
    stem = leaf.edge.length
    crown = stem * rng.uniform(0.3, 0.8)
    n = rng.randint(3, 25)
    sub = dendropy.simulate.treesim.birth_death_tree(1.0, 0.0, num_extant_tips=n, rng=rng)
    sub.calc_node_ages()
    f = crown / sub.seed_node.age
    for e in sub.postorder_edge_iter():
        if e.length:
            e.length *= f
    for i, l in enumerate(sub.leaf_node_iter()):
        l.taxon = full.taxon_namespace.require_taxon(label=f"{g}_u{i:02d}")
    kids = [frozenset(x.taxon.label for x in c.leaf_iter()) for c in sub.seed_node.child_node_iter()]
    truth[g] = sorted(sorted(k) for k in kids)
    leaf.edge.length = stem - crown
    leaf.taxon = None
    for c in list(sub.seed_node.child_node_iter()):
        leaf.add_child(c)

# sequences: K2P + gamma, lognormal rate per branch (non-clocklike)
site_rates = nrng.gamma(ALPHA, 1.0 / ALPHA, L)
beta = 1.0 / (KAPPA + 2.0)
alpha_ts = KAPPA * beta
ts_of = np.array([2, 3, 0, 1])  # A<->G, C<->T with A0 C1 G2 T3
seq = {full.seed_node: nrng.integers(0, 4, L)}
for node in full.preorder_node_iter():
    if node.parent_node is None:
        continue
    t = node.edge.length * RATE * float(nrng.lognormal(0, 0.6)) * site_rates
    e1 = np.exp(-4 * beta * t)
    e2 = np.exp(-2 * (alpha_ts + beta) * t)
    p_ts = 0.25 + 0.25 * e1 - 0.5 * e2
    p_tv = 0.25 - 0.25 * e1
    parent = seq[node.parent_node]
    u = nrng.random(L)
    child = parent.copy()
    is_ts = u < p_ts
    is_tv = (u >= p_ts) & (u < p_ts + 2 * p_tv)
    child[is_ts] = ts_of[parent[is_ts]]
    tv_choice = nrng.integers(0, 2, L)
    tv_targets = np.where(parent % 2 == 0, np.array([1, 3])[tv_choice], np.array([0, 2])[tv_choice])
    child[is_tv] = tv_targets[is_tv]
    seq[node] = child

reps = []
for leaf in full.leaf_node_iter():
    s = "".join("ACGT"[b] for b in seq[leaf])
    g = leaf.taxon.label.split("_")[0]
    if rng.random() < 0.08:  # fragment
        a = rng.randint(0, L - 350)
        s = "-" * a + s[a:a + 350] + "-" * (L - a - 350)
    reps.append({"unit_id": leaf.taxon.label, "genus": g, "sequence": s})
# misassigned units: relabel a unit as belonging to another genus
for _ in range(8):
    r = rng.choice(reps)
    r["genus"] = f"G{rng.randrange(N_GENERA):02d}"

(out / "units").mkdir(parents=True, exist_ok=True)
(out / "backbone").mkdir(exist_ok=True)
write_tsv(out / "units" / "representatives.tsv", reps, ["unit_id", "genus", "sequence"])
bb.write(path=str(out / "backbone" / "backbone.tre"), schema="newick", suppress_rooting=True)
json.dump(truth, open(out / "truth.json", "w"))

# shard and build exactly as the pipeline does
cfg = yaml.safe_load(open("config/config.yaml"))
marker_cfg = cfg["markers"]["coi"]
stems = stem_ages(bb)
by_genus = {}
for r in reps:
    by_genus.setdefault(r["genus"], []).append(r)
available = {g for g in by_genus if g in stems}
rows, prov = [], []
for g in sorted(available):
    units = sorted(by_genus[g], key=lambda r: r["unit_id"])
    if len(units) < 2:
        continue
    picked = []
    for sis in nearest_genera(bb, g, available):
        for u in sorted(by_genus[sis], key=lambda r: r["unit_id"]):
            if len(picked) < cfg["n_outgroup"]:
                picked.append(u)
    rows += [{"genus": g, "role": "in", "unit_id": u["unit_id"], "stem_age": stems[g], "sequence": u["sequence"]} for u in units]
    rows += [{"genus": g, "role": "out", "unit_id": u["unit_id"], "stem_age": stems[g], "sequence": u["sequence"]} for u in picked]
    tree, p = build_genus(g, [(u["unit_id"], u["sequence"]) for u in units],
                          [(u["unit_id"], u["sequence"]) for u in picked], stems[g], marker_cfg, {}, {})
    prov.append({**p, "newick": tree_to_newick(tree)})
write_tsv(out / "genera" / "shard_0.tsv", rows, ["genus", "role", "unit_id", "stem_age", "sequence"])
write_tsv(out / "subtrees" / "shard_0.tsv", prov, list(PROVENANCE_COLUMNS) + ["newick"])
yaml.safe_dump({"outdir": str(out), "backbone": str(out / "backbone" / "backbone.tre")}, open(out / "run.yaml", "w"))
from collections import Counter
print(Counter(p["rooting"] for p in prov))
