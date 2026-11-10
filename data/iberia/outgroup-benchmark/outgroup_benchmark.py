"""Benchmark of outgroup choices for rooting genus subtrees (issue #19).

Per genus, the ingroup and every candidate outgroup unit are aligned once and
a K2P matrix is computed; each method then builds BioNJ on its part of that
matrix, so differences between methods are not alignment differences.

Methods
  ref            the run's own two outgroup units (reference root on real data)
  sister_random  two other random units of the first sister genus (noise floor)
  far_data       the two closest units, from up to four sister clades, that are
                 farther from every ingroup unit than the two most distant ingroup
                 units are from each other
  family_random  two random units from other families of the same order
  family_near    the two units from other families of the same order closest to the ingroup
  tolerant10     ten units spread over the sister clades, rooted at the edge with
                 the fewest units on the wrong side (minimum violation)
  mad, midpoint  no outgroup

Strict methods (all with an outgroup but tolerant10) fail when the ingroup is
not a split; they then give no root.
"""
import argparse
import csv
import glob
import json
import random
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np

from timegraft import distances
from timegraft.backbone import read_backbone, stem_ages
from timegraft.bionj import bionj
import subprocess, tempfile
from pathlib import Path
from timegraft.seqio import read_fasta, write_fasta


def align_mafft(items):
    """FFT-NS-2: fast, and every method shares the same alignment."""
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "in.fasta"
        write_fasta(f, [(f"s{i}", s.replace("-", "")) for i, (_, s) in enumerate(items)])
        res = subprocess.run(["mafft", "--retree", "2", "--maxiterate", "0", "--quiet", "--thread", "1", str(f)],
                             check=True, capture_output=True, text=True)
        (Path(tmp) / "out.fasta").write_text(res.stdout)
        al = dict(read_fasta(Path(tmp) / "out.fasta"))
    return [(n, al[f"s{i}"]) for i, (n, _) in enumerate(items)]
from timegraft.dating import mean_root_to_tip
from timegraft.rooting import midpoint, root_on_outgroup
from timegraft.rooting_diagnostics import (copy_tree, describe, leaf_sets, mad_root, pick_spread,
                                           reroot, rogue_counts, root_split, sister_clades)

RATE = 0.0115
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
STRICT = ["ref", "sister_random", "far_data", "family_random", "family_near"]
METHODS = STRICT + ["tolerant10", "mad", "midpoint"]


def mv_rooted(tree, ingroup, outgroup):
    """Root at the side with the fewest misplaced units, then drop the outgroup."""
    labels = ingroup | outgroup
    best = None
    for node, below in leaf_sets(tree).items():
        if node.parent_node is None:
            continue
        for side, child_side in ((below, True), (labels - below, False)):
            if not ingroup & side:
                continue
            cost = (len(ingroup - side) + len(outgroup & side), -len(ingroup & side))
            if best is None or cost < best[0]:
                best = (cost, below, node.edge.length or 0.0)
    _, below, length = best
    t = reroot(copy_tree(tree), below, length / 2.0)
    t.prune_taxa_with_labels(list(outgroup))
    t.suppress_unifurcations()
    return t, best[0][0]


def evaluate(rooted, keep):
    if rooted is None:
        return {}
    d = describe(rooted, RATE)
    return {"split": root_split(rooted, keep), "crown": d["crown"], "pd_rel": d["pd_rel"]}


def work(task):
    rng = random.Random(task["genus"])
    I = task["in"]
    ingroup = frozenset(n for n, _ in I)
    cands = {}
    for key in ("base", "sister", "far", "family", "near10"):
        for n, s in task[key]:
            if n not in ingroup:
                cands[n] = s
    items = align_mafft(I + list(cands.items()))
    labels = [n for n, _ in items]
    idx = {n: i for i, n in enumerate(labels)}
    D = distances.matrix([s for _, s in items], method="k2p")
    iin = [idx[n] for n, _ in I]
    sub = D[np.ix_(iin, iin)]
    dmax = sub.max()
    dmin = {n: D[idx[n], iin].min() for n in cands}

    def strict(outs):
        outs = [o for o in outs if o in idx]
        if len(outs) < 1:
            return None, None
        sel = list(iin) + [idx[o] for o in outs]
        tree = bionj([labels[i] for i in sel], D[np.ix_(sel, sel)])
        out = frozenset(outs)
        if rogue_counts(tree, ingroup, out)["rogue_in"] != 0:
            return "fail", float(np.mean([dmin[o] for o in outs]))
        rooted, _ = root_on_outgroup(copy_tree(tree), ingroup, out)
        return rooted, float(np.mean([dmin[o] for o in outs]))

    choice = {}
    choice["ref"] = [n for n, _ in task["base"]]
    sis = [n for n, _ in task["sister"] if n not in choice["ref"]]
    choice["sister_random"] = rng.sample(sis, 2) if len(sis) >= 2 else []
    far = sorted((n for n, _ in task["far"] if n in dmin and dmin[n] > dmax), key=lambda n: dmin[n])
    choice["far_data"] = far[:2]
    fam = [n for n, _ in task["family"]]
    choice["family_random"] = rng.sample(fam, 2) if len(fam) >= 2 else []
    choice["family_near"] = sorted(fam, key=lambda n: dmin[n])[:2] if len(fam) >= 2 else []

    rows = {}
    for m in STRICT:
        if len(choice[m]) < 2:
            rows[m] = {"status": "no candidates"}
            continue
        rooted, dist = strict(choice[m])
        if rooted == "fail":
            rows[m] = {"status": "fail", "dist": dist}
        else:
            rows[m] = {"status": "ok", "dist": dist, **evaluate(rooted, ingroup)}
    near = [n for n, _ in task["near10"]]
    if len(near) >= 3:
        sel = list(iin) + [idx[o] for o in near]
        tree = bionj([labels[i] for i in sel], D[np.ix_(sel, sel)])
        rooted, viol = mv_rooted(tree, ingroup, frozenset(near))
        rows["tolerant10"] = {"status": "ok", "violations": viol,
                              "dist": float(np.mean([dmin[o] for o in near])), **evaluate(rooted, ingroup)}
    else:
        rows["tolerant10"] = {"status": "no candidates"}
    itree = bionj([labels[i] for i in iin], sub)
    side, to_child, _, _ = mad_root(copy_tree(itree))
    rows["mad"] = {"status": "ok", **evaluate(reroot(copy_tree(itree), side, to_child), ingroup)}
    rows["midpoint"] = {"status": "ok", **evaluate(midpoint(copy_tree(itree)), ingroup)}

    truth = task.get("truth")
    if truth is not None:
        ref = frozenset(frozenset(k) & ingroup for k in truth)
        if any(not k for k in ref):
            ref = None
    else:
        ref = rows["ref"].get("split")
    out = []
    for m in METHODS:
        r = rows[m]
        row = {"genus": task["genus"], "set": task["set"], "n_in": len(I), "method": m, "status": r["status"],
               "dist": r.get("dist", ""), "violations": r.get("violations", "")}
        if r["status"] == "ok" and ref is not None:
            row["same_root"] = r["split"] == ref
            if rows["ref"].get("crown"):
                row["crown_ratio"] = r["crown"] / rows["ref"]["crown"]
                row["pd_ratio"] = r["pd_rel"] / rows["ref"]["pd_rel"]
        out.append(row)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--taxonomy", required=True, help="genus, family, order")
    ap.add_argument("--truth", help="truth.json from the simulation")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()

    run = a.run
    pg = {r["genus"]: r for r in rd(f"{run}/diagnostics/rooting/per_genus.tsv")}
    rp = {}
    try:
        rp = {r["genus"]: float(r["random_pair_ok"]) for r in rd(f"{run}/diagnostics/rooting/random_pairs.tsv")}
    except FileNotFoundError:
        pass
    truth = json.load(open(a.truth)) if a.truth else None
    tax = {r["genus"]: r for r in rd(a.taxonomy)}
    reps = rd(f"{run}/units/representatives.tsv")
    by_genus = defaultdict(list)
    for r in reps:
        by_genus[r["genus"]].append(r)
    bb = read_backbone(f"{run}/backbone/backbone.tre")
    avail = {g for g in by_genus if g in stem_ages(bb)}
    shard = defaultdict(lambda: {"in": [], "out": []})
    for p in glob.glob(f"{run}/genera/shard_*.tsv"):
        for r in rd(p):
            shard[r["genus"]][r["role"]].append((r["unit_id"], r["sequence"]))

    tasks = []
    for g, d in sorted(shard.items()):
        if len(d["in"]) < 3 or len(d["out"]) < 2:
            continue
        if truth is not None:
            sset = "simulated"
        elif rp.get(g) == 1.0:
            sset = "stable"
        elif pg[g]["split_ok"] == "False":
            sset = "failed"
        else:
            continue
        rng = random.Random(g + "pool")
        clades = sister_clades(bb, g, avail, levels=4)
        if not clades:
            continue
        sister = [(r["unit_id"], r["sequence"]) for r in by_genus[clades[0][0]]]
        sister = rng.sample(sister, min(60, len(sister)))
        far = [(r["unit_id"], r["sequence"]) for c in clades for x in c for r in by_genus[x]]
        far = rng.sample(far, min(150, len(far)))
        t = tax.get(g, {})
        fam = [(r["unit_id"], r["sequence"]) for x in avail if x != g and tax.get(x, {}).get("order")
               and tax[x]["order"] == t.get("order") and tax[x]["family"] != t.get("family") for r in by_genus[x]]
        fam = rng.sample(sorted(fam), min(100, len(fam)))
        near10 = pick_spread(sister_clades(bb, g, avail, levels=3), by_genus, 10)
        tasks.append({"genus": g, "set": sset, "in": d["in"], "base": d["out"], "sister": sister,
                      "far": far, "family": fam, "near10": near10,
                      "truth": truth.get(g) if truth else None})
    print(Counter(t["set"] for t in tasks), flush=True)
    with Pool(a.threads) as p:
        rows = [r for rs in p.imap_unordered(work, tasks, chunksize=1) for r in rs]
    cols = ["genus", "set", "n_in", "method", "status", "dist", "violations", "same_root", "crown_ratio", "pd_ratio"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["set"], r["genus"], METHODS.index(r["method"]))))


if __name__ == "__main__":
    main()
