"""Evaluate the release gate of the consistency layer (#56): every metric of GATE.md, with its status.

    python evaluate.py BENCH OUT

BENCH holds ``main/results/iberia`` and ``branch/results/iberia`` (the two full runs), ``main_run.json`` and
``branch_run.json`` (time and memory), ``tests.txt`` (the last line of pytest) and ``exp/`` as
``run_experiments.sh`` leaves it. OUT gets ``metrics.tsv`` (the gate) and a table per experiment.
"""

import csv
import filecmp
import json
import random
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from benchlib import Backbone, PDTree, factor, quantile, read, spearman, write

BENCH, OUT = Path(sys.argv[1]), Path(sys.argv[2])
MAIN, BRANCH, EXP = BENCH / "main/results/iberia", BENCH / "branch/results/iberia", BENCH / "exp"
SEEDS = [1, 2, 3, 4, 5]
HERE = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)
metrics = []


def gate(number, name, legacy, default, threshold, passed, note=""):
    status = "report" if passed is None else ("pass" if passed else "FAIL")
    metrics.append({"metric": number, "name": name, "legacy": legacy, "default": default, "threshold": threshold,
                    "status": status, "note": note})
    print(f"{number:>3} {status:6} {name}: legacy {legacy}; default {default}; threshold {threshold} {note}")


def fmt(x, digits=3):
    return "" if x is None else f"{x:.{digits}f}"


# 0. Soundness: the legacy settings on consistency-layer against the run from main ----------------------
main_bb, legacy_bb = Backbone(MAIN / "backbone/backbone.tre"), Backbone(EXP / "legacy_clean/backbone.tre")
a, b = main_bb.node_ages(), legacy_bb.node_ages()
same_tree = set(a) == set(b)
age_diff = max((abs(a[k] - b[k]) for k in a if k in b), default=0.0)
main_dec = {r["unit_id"]: r["decision"] for r in read(MAIN / "units/screened.tsv")}
legacy_dec = {r["unit_id"]: r["decision"] for r in read(EXP / "legacy_units/screened.tsv")}
main_crowns = {r["genus"]: r["crown_age"] for r in read(MAIN / "backbone/crown_ages.tsv") if r["crown_age"]}
legacy_crowns = {r["genus"]: r["crown_age"] for r in read(EXP / "legacy_clean/crown_ages.tsv") if r["crown_age"]}
crown_diff = sum(1 for g in set(main_crowns) | set(legacy_crowns)
                 if g not in main_crowns or g not in legacy_crowns
                 or abs(float(main_crowns[g]) - float(legacy_crowns[g])) > 1e-3)
gate("0", "legacy settings reproduce main", "",
     f"clades {'same' if same_tree else 'differ'}; largest age difference {age_diff:.4f} Ma; "
     f"{sum(main_dec.get(u) != d for u, d in legacy_dec.items())} of {len(legacy_dec)} screening decisions differ; "
     f"{crown_diff} of {len(main_crowns)} TimeTree crowns differ",
     "identical", same_tree and age_diff < 1e-3 and main_dec == legacy_dec and crown_diff == 0,
     f"backbone.tre byte-identical: {filecmp.cmp(MAIN / 'backbone/backbone.tre', EXP / 'legacy_clean/backbone.tre', False)}")

# the harness against the workflow: the default backbone of the full run and of backbone.py
default_bb, e2e_bb = Backbone(EXP / "default_clean/backbone.tre"), Backbone(BRANCH / "backbone/backbone.tre")
harness_same = filecmp.cmp(EXP / "default_clean/backbone.tre", BRANCH / "backbone/backbone.tre", False)

genera = {r["genus"]: r for r in read(EXP / "default_clean/genera.tsv")}
klass = {g: r["records_class"] for g, r in genera.items()}
phylum = {g: r["records_phylum"] for g, r in genera.items()}

# 1a. Onchidella ------------------------------------------------------------------------------------------


def company_phyla(bb, genus):
    return Counter(phylum.get(g, "?") for g in bb.smallest_clade_with_company(genus))


def unit_company(tree_path, reps_path, genus):
    """Phyla of the genera of the other units in the smallest clade of the final tree around the genus's units."""
    import dendropy

    unit_genus = {r["unit_id"]: r["genus"] for r in read(reps_path)}
    tree = dendropy.Tree.get(path=str(tree_path), schema="newick", preserve_underscores=True, rooting="force-rooted")
    own = {u for u, g in unit_genus.items() if g == genus}
    leaves = [leaf for leaf in tree.leaf_node_iter() if leaf.taxon.label in own]
    if not leaves:
        return Counter()
    node = leaves[0].parent_node
    while node is not None:
        others = {unit_genus.get(leaf.taxon.label, "") for leaf in node.leaf_iter()} - {genus, ""}
        if others:
            return Counter(phylum.get(g, "?") for g in others)
        node = node.parent_node
    return Counter()


legacy_on = unit_company(MAIN / "timegraft.tre", MAIN / "units/representatives.tsv", "Onchidella")
default_on = unit_company(BRANCH / "timegraft.tre", BRANCH / "units/representatives.tsv", "Onchidella")
gate("1a", "Onchidella's nearest genera in the final tree", dict(legacy_on), dict(default_on), "all Mollusca",
     bool(default_on) and set(default_on) == {"Mollusca"},
     f"backbone: legacy {dict(company_phyla(main_bb, 'Onchidella'))}, default {dict(company_phyla(e2e_bb, 'Onchidella'))}")

# 1b-1g. Injected misplacements -----------------------------------------------------------------------------
def injection(name):
    """Every injected genus of the replicates under EXP/name, and the side effects per replicate."""
    rows, collateral = [], []
    for seed in SEEDS:
        folder = EXP / f"{name}/s{seed}"
        manifest = read(folder / "manifest.tsv")
        injected = {m["genus"] for m in manifest}
        conc = {r["genus"]: r for r in read(folder / "default/concordance.tsv")}
        gen = {r["genus"]: r for r in read(folder / "default/genera.tsv")}
        d_bb, l_bb = Backbone(folder / "default/backbone.tre"), Backbone(folder / "legacy/backbone.tre")
        for m in manifest:
            g = m["genus"]
            c = conc.get(g, {})

            def own_class(bb):
                if g not in bb.leaf:
                    return None
                others = bb.smallest_clade_with_company(g)
                return bool(others) and all(klass.get(o) == klass[g] for o in others)

            def stem_factor(bb, clean):
                return factor(bb.stem_age(g), clean.stem_age(g)) if g in bb.leaf and g in clean.leaf else None

            rows.append({
                **m, "n_evaluable": c.get("n_evaluable", ""), "n_flags": c.get("n_flags", ""),
                "lineage_flag": c.get("lineage_flag", ""), "merged_flag": c.get("merged_flag", ""),
                "neighbour_flag": c.get("neighbour_flag", ""), "neighbour_hits": c.get("neighbour_hits", ""),
                "concordance_flag": c.get("concordance_flag", ""), "placement_action": c.get("placement_action", ""),
                "reanchor_target": gen.get(g, {}).get("reanchor_target", ""),
                "default_in_backbone": g in d_bb.leaf, "default_own_class": own_class(d_bb),
                "legacy_in_backbone": g in l_bb.leaf, "legacy_own_class": own_class(l_bb),
                "default_stem_factor": fmt(stem_factor(d_bb, default_bb)),
                "legacy_stem_factor": fmt(stem_factor(l_bb, legacy_bb)),
            })
        false_flags = [g for g, c in conc.items() if c["concordance_flag"] == "flag" and g not in injected
                       and g != "Onchidella"]
        moved = {}
        for arm, bb, clean in (("legacy", l_bb, legacy_bb), ("default", d_bb, default_bb)):
            shared = [g for g in bb.leaf if g in clean.leaf and g not in injected]
            moved[arm] = sum(1 for g in shared if factor(bb.stem_age(g), clean.stem_age(g)) > 1.10)
            moved[arm + "_n"] = len(shared)
        collateral.append({"seed": seed, "injected": len(manifest), "false_flags": len(false_flags),
                           "false_flag_genera": ";".join(false_flags),
                           "legacy_other_stems_moved": moved["legacy"], "legacy_other_genera": moved["legacy_n"],
                           "default_other_stems_moved": moved["default"], "default_other_genera": moved["default_n"]})
    return rows, collateral


def injection_summary(rows, collateral):
    flagged = [r for r in rows if r["concordance_flag"] == "flag"]
    corrected = [r for r in flagged if r["default_own_class"]]
    near = [r for r in corrected if r["default_stem_factor"] and float(r["default_stem_factor"]) <= 1.5]
    return {
        "injected": len(rows),
        "two_checks_evaluable": sum(1 for r in rows if r["n_evaluable"] and int(r["n_evaluable"]) >= 2),
        "flagged": len(flagged), "review": sum(1 for r in rows if r["concordance_flag"] == "review"),
        "back_in_own_class": len(corrected), "stem_within_1.5": len(near),
        "legacy_in_own_class": sum(1 for r in rows if r["legacy_own_class"]),
        "legacy_stem_within_1.5": sum(1 for r in rows if r["legacy_stem_factor"] and float(r["legacy_stem_factor"]) <= 1.5),
        "false_flags": sum(c["false_flags"] for c in collateral),
        "other_stems_moved_default": sum(c["default_other_stems_moved"] for c in collateral),
        "other_stems_moved_legacy": sum(c["legacy_other_stems_moved"] for c in collateral),
    }


rows, collateral = injection("inject")
write(OUT / "injection_genera.tsv", rows)
write(OUT / "injection_replicates.tsv", collateral)
flagged = [r for r in rows if r["concordance_flag"] == "flag"]
evaluable = [r for r in rows if r["n_evaluable"] and int(r["n_evaluable"]) >= 2]
corrected = [r for r in flagged if r["default_own_class"]]
for kind in ("phylum", "class"):
    sub = [r for r in rows if r["kind"] == kind]
    print(f"    injected to another {kind}: {len(sub)}; flagged {sum(r['concordance_flag'] == 'flag' for r in sub)}; "
          f"with two checks {sum(1 for r in sub if r in evaluable)}; "
          f"back in own class {sum(1 for r in sub if r['default_own_class'])}; "
          f"legacy in own class {sum(1 for r in sub if r['legacy_own_class'])}")
share = len([r for r in evaluable if r["concordance_flag"] == "flag"]) / len(evaluable) if evaluable else 0.0
gate("1b", "injected misplacements with two checks: flagged", "0 (no check)",
     f"{share:.3f} ({sum(r['concordance_flag'] == 'flag' for r in evaluable)} of {len(evaluable)})", ">= 0.90",
     share >= 0.90)
gate("1c", "injected misplacements, all: flagged", "0 (no check)",
     f"{len(flagged) / len(rows):.3f} ({len(flagged)} of {len(rows)})", ">= 0.50", len(flagged) / len(rows) >= 0.50)
share = len(corrected) / len(flagged) if flagged else 0.0
gate("1d", "flagged injected genera back in their own class",
     f"{sum(1 for r in rows if r['legacy_own_class'])} of {len(rows)} in own class",
     f"{share:.3f} ({len(corrected)} of {len(flagged)})", ">= 0.90", share >= 0.90)
mean_false = statistics.mean(c["false_flags"] for c in collateral)
gate("1e", "genera without an injected error flagged, per replicate", "",
     f"{mean_false:.1f} ({[c['false_flags'] for c in collateral]})", "<= 1", mean_false <= 1)
near = [r for r in corrected if r["default_stem_factor"] and float(r["default_stem_factor"]) <= 1.5]
share = len(near) / len(corrected) if corrected else 0.0
legacy_far = [r for r in rows if r["legacy_stem_factor"] and float(r["legacy_stem_factor"]) > 1.5]
gate("1f", "corrected genera with a stem age within 1.5 of the one without the error",
     f"{len(rows) - len(legacy_far)} of {len(rows)} injected within 1.5",
     f"{share:.3f} ({len(near)} of {len(corrected)})", ">= 0.80", share >= 0.80)
lm, dm = sum(c["legacy_other_stems_moved"] for c in collateral), sum(c["default_other_stems_moved"] for c in collateral)
gate("1g", "other genera whose stem age moves by more than 10%, five replicates together", lm, dm,
     "<= half of legacy", dm <= 0.5 * lm)

# The two variants added after the first results (not part of the gate): see run_posthoc.sh
variants = [{"injection": "merge (the gate)", **injection_summary(rows, collateral)}]
for name, label in (("inject_own", "merge, host's own species left out (as Onchidella)"),
                    ("inject_name", "name only: no merged NCBI ids, no species moved")):
    if (EXP / name / "s5/default/metadata.json").exists():
        v_rows, v_collateral = injection(name)
        write(OUT / f"{name}_genera.tsv", v_rows)
        variants.append({"injection": label, **injection_summary(v_rows, v_collateral)})
write(OUT / "injection_variants.tsv", variants)
for v in variants:
    print("   ", v)

# 2. Backbone conflicts ----------------------------------------------------------------------------------
lmeta = json.load(open(EXP / "legacy_clean/metadata.json"))["ages"]
dmeta = json.load(open(EXP / "default_clean/metadata.json"))["ages"]
lf, df = lmeta["final"], dmeta["final"]
gate("2a", "TimeTree calibrations off from the final age by more than 1.5", f"{lf['severe']} of {lf['n']}",
     f"{df['severe']} of {df['n']}", "<= half of legacy", df["severe"] <= 0.5 * lf["severe"],
     f"largest factor: legacy {lf['max_factor']}, default {df['max_factor']}; 90%: {lf['p90_factor']}, {df['p90_factor']}")
gate("2b", "RMS of calibration residuals over their uncertainty", lf["weighted_rms"], df["weighted_rms"],
     "lower than legacy", df["weighted_rms"] < lf["weighted_rms"])
final_tree = Backbone(BRANCH / "timegraft.tre")
gate("2c", "nodes older than their parent", f"backbone {main_bb.negative}",
     f"backbone {e2e_bb.negative}, final tree {final_tree.negative}", "0",
     e2e_bb.negative == 0 and final_tree.negative == 0)

# 3. Crown tails -----------------------------------------------------------------------------------------


def crown_tails(path):
    out = Counter()
    for r in read(path):
        if r.get("status", "subtree") != "subtree" or not r["crown_age"]:
            continue
        crown = float(r["crown_age"])
        out["genera"] += 1
        if r.get("barcode_age"):
            ratio = crown / float(r["barcode_age"])
            out["over20"] += ratio > 20
            out["over10"] += ratio > 10
            out["young"] += r["crown_source"].startswith("timetree") and ratio < 1 / 3
        if r.get("root_to_tip") and float(r["root_to_tip"]) < 0.02 and crown > 20:
            out["shallow_old"] += 1
    return out


lt, dt = crown_tails(MAIN / "provenance.tsv"), crown_tails(BRANCH / "provenance.tsv")
gate("3a", "crown age more than 20 times the barcode age", lt["over20"], dt["over20"], "0", dt["over20"] == 0,
     f"of {lt['genera']} and {dt['genera']} genera with a subtree")
gate("3b", "crown age more than 10 times the barcode age", lt["over10"], dt["over10"], "fewer than legacy",
     dt["over10"] < lt["over10"])
gate("3c", "TimeTree crown used that is less than a third of the barcode age", lt["young"], dt["young"], "0",
     dt["young"] == 0)
gate("3d", "root-to-tip below 0.02 and crown above 20 Ma", lt["shallow_old"], dt["shallow_old"], "none", None)

# 4. Hold-out --------------------------------------------------------------------------------------------
holdout = OUT / "holdout.tsv"  # kept gzipped in the repository
subprocess.run([sys.executable, str(HERE.parent / "agefit/holdout.py"), str(EXP / "default_clean"), "50", str(holdout)],
               check=True, stdout=subprocess.DEVNULL)
errors = defaultdict(list)
for r in read(holdout):
    errors[r["method"]].append(factor(float(r["age"]), float(r["timetree_age"])))
lmed, dmed = statistics.median(errors["precedence"]), statistics.median(errors["fit"])
l90, d90 = quantile(errors["precedence"], 0.9), quantile(errors["fit"], 0.9)
gate("4a", "held-out TimeTree ages: median error factor", fmt(lmed), fmt(dmed), "<= legacy x 1.02", dmed <= lmed * 1.02,
     f"{len(errors['fit'])} held-out ages")
gate("4b", "held-out TimeTree ages: 90% quantile of the error factor", fmt(l90), fmt(d90), "<= legacy x 1.02",
     d90 <= l90 * 1.02)

# 5. Perturbations ----------------------------------------------------------------------------------------
norequery = Backbone(EXP / "default_norequery/backbone.tre")
nr = norequery.node_ages()
nr_diff = max(abs(nr[k] - v) for k, v in default_bb.node_ages().items())
moves = {arm: [] for arm in ("legacy", "default")}
stems = {arm: [] for arm in ("legacy", "default")}
noise_rows, pd_rho = [], {arm: [] for arm in ("legacy", "default")}
rng = random.Random(56)
genus_samples = None
for seed in SEEDS:
    for arm, clean in (("legacy", legacy_bb), ("default", norequery)):
        bb = Backbone(EXP / f"noise/s{seed}/{arm}/backbone.tre")
        ages, ref = bb.node_ages(), clean.node_ages()
        f = [factor(ages[k], ref[k]) for k in ages if k in ref and ref[k] > 0 and ages[k] > 0]
        moves[arm] += f
        s = [factor(bb.stem_age(g), clean.stem_age(g)) for g in bb.leaf if g in clean.leaf]
        stems[arm] += s
        # PD of genus samples on the backbone: is the ranking of samples the same with and without the noise?
        if genus_samples is None:
            pool = sorted(set(legacy_bb.leaf) & set(default_bb.leaf))
            genus_samples = [rng.sample(pool, 100) for _ in range(200)]
        t0, t1 = PDTree(EXP / ("legacy_clean" if arm == "legacy" else "default_norequery") / "backbone.tre"), \
            PDTree(EXP / f"noise/s{seed}/{arm}/backbone.tre")
        rho = spearman([t0.pd(s_) for s_ in genus_samples], [t1.pd(s_) for s_ in genus_samples])
        pd_rho[arm].append(rho)
        noise_rows.append({"seed": seed, "arm": arm, "nodes": len(f), "same_clades": len(f) == len(ages),
                           "median_factor": fmt(statistics.median(f)), "p90_factor": fmt(quantile(f, 0.9)),
                           "share_over_1.5": fmt(sum(x > 1.5 for x in f) / len(f)),
                           "stem_p90_factor": fmt(quantile(s, 0.9)), "pd_rank_correlation": fmt(rho, 4)})
write(OUT / "noise.tsv", noise_rows)
l, d = quantile(moves["legacy"], 0.9), quantile(moves["default"], 0.9)
gate("5a", "ages perturbed: 90% quantile of the factor by which node ages move", fmt(l), fmt(d), "<= legacy x 1.05",
     d <= l * 1.05, f"medians {fmt(statistics.median(moves['legacy']))}, {fmt(statistics.median(moves['default']))}; "
                    f"default without re-queries differs from the default by at most {nr_diff:.4f} Ma")
l, d = (sum(x > 1.5 for x in moves[arm]) / len(moves[arm]) for arm in ("legacy", "default"))
gate("5b", "ages perturbed: share of node ages that move by more than 1.5", fmt(l), fmt(d), "<= legacy + 0.01",
     d <= l + 0.01, f"stems, 90% quantile: {fmt(quantile(stems['legacy'], 0.9))}, {fmt(quantile(stems['default'], 0.9))}")

moves = {arm: [] for arm in ("legacy", "default")}
stems = {arm: [] for arm in ("legacy", "default")}
sub_rows = []
for seed in SEEDS:
    for arm, clean in (("legacy", legacy_bb), ("default", default_bb)):
        bb = Backbone(EXP / f"sub/s{seed}/{arm}/backbone.tre")
        f = [factor(age, clean.mrca_age(tips)) for tips, age in bb.node_ages().items()
             if tips <= set(clean.leaf) and age > 0]
        s = []
        for g in bb.leaf:
            node = bb.parent[bb.leaf[g]]
            if bb.tips[node] <= set(clean.leaf) and bb.age[node] > 0:
                s.append(factor(bb.age[node], clean.mrca_age(bb.tips[node])))
        moves[arm] += f
        stems[arm] += s
        sub_rows.append({"seed": seed, "arm": arm, "genera": len(bb.leaf), "nodes": len(f),
                         "median_factor": fmt(statistics.median(f)), "p90_factor": fmt(quantile(f, 0.9)),
                         "share_over_1.5": fmt(sum(x > 1.5 for x in f) / len(f)),
                         "stem_p90_factor": fmt(quantile(s, 0.9)),
                         "stem_share_over_1.5": fmt(sum(x > 1.5 for x in s) / len(s))})
write(OUT / "subsample.tsv", sub_rows)
l, d = quantile(moves["legacy"], 0.9), quantile(moves["default"], 0.9)
gate("5c", "half the genera dropped: 90% quantile of the factor by which clade ages move", fmt(l), fmt(d),
     "<= legacy x 1.05", d <= l * 1.05,
     f"share over 1.5: {fmt(sum(x > 1.5 for x in moves['legacy']) / len(moves['legacy']))}, "
     f"{fmt(sum(x > 1.5 for x in moves['default']) / len(moves['default']))}")
l, d = quantile(stems["legacy"], 0.9), quantile(stems["default"], 0.9)
gate("5d", "half the genera dropped: the same for genus stem ages", fmt(l), fmt(d), "<= legacy x 1.05", d <= l * 1.05)

# 6. Diversity ---------------------------------------------------------------------------------------------
old, new = PDTree(MAIN / "timegraft.tre"), PDTree(BRANCH / "timegraft.tre")
unit_genus = {r["unit_id"]: r["genus"] for r in read(BRANCH / "units/representatives.tsv")}
shared = sorted(set(old.tip) & set(new.tip))
by_genus = defaultdict(list)
for u in shared:
    if u in unit_genus:
        by_genus[unit_genus[u]].append(u)
rng = random.Random(6)
all_genera, rich = sorted(by_genus), sorted(g for g, us in by_genus.items() if len(us) >= 5)
designs = {
    "150 units at random": [rng.sample(shared, 150) for _ in range(200)],
    "100 genera, one unit each": [[rng.choice(by_genus[g]) for g in rng.sample(all_genera, 100)] for _ in range(200)],
    "100 genera, five units each": [[u for g in rng.sample(rich, 100) for u in rng.sample(by_genus[g], 5)]
                                    for _ in range(200)],
}
pd_rows, rhos = [], []
for name, samples in designs.items():
    x, y = [old.pd(s) for s in samples], [new.pd(s) for s in samples]
    rho = spearman(x, y)
    rhos.append(rho)
    pd_rows.append({"design": name, "samples": len(samples), "rank_correlation": fmt(rho, 4),
                    "median_pd_legacy": fmt(statistics.median(x), 1), "median_pd_default": fmt(statistics.median(y), 1),
                    "median_ratio": fmt(statistics.median(b / a for a, b in zip(x, y)), 4),
                    "ratio_p05": fmt(quantile([b / a for a, b in zip(x, y)], 0.05), 4),
                    "ratio_p95": fmt(quantile([b / a for a, b in zip(x, y)], 0.95), 4)})
write(OUT / "pd.tsv", pd_rows)
gate("6a", "PD of the same samples on the two final trees: rank correlation", "",
     "; ".join(f"{r['design']}: {r['rank_correlation']}" for r in pd_rows), ">= 0.90 in each design", min(rhos) >= 0.90,
     f"{len(shared)} units in both trees; {len(old.tip)} and {len(new.tip)} tips")

# 6b. stem ages that differ between the runs, and whether the outputs say why
main_gen = {r["genus"]: r for r in read(MAIN / "backbone/genera.tsv")}
e2e_gen = {r["genus"]: r for r in read(BRANCH / "backbone/genera.tsv")}
node_rows = read(BRANCH / "backbone/node_ages.tsv")
order = {f"n{i}": n for i, n in enumerate(e2e_bb.tree.preorder_node_iter())}
row_of = {order[r["node"]]: r for r in node_rows if r["node"] in order}
legacy_anchor = {r["genus"]: r.get("anchor_node", "") for r in read(EXP / "legacy_clean/genera.tsv")}
cause_rows = []
for g in sorted(set(main_bb.leaf) & set(e2e_bb.leaf)):
    f = factor(e2e_bb.stem_age(g), main_bb.stem_age(g))
    if f <= 1.10:
        continue
    r = e2e_gen[g]
    causes = []
    if r.get("placement_action") in ("reanchored",):
        causes.append("reanchored")
    if r.get("anchor_node", "") != legacy_anchor.get(g, ""):
        causes.append("anchor")
    if r.get("join_ncbi_source") == "col" and r.get("join_ncbi_id") != main_gen.get(g, {}).get("ncbi_id"):
        causes.append("identifier")
    fitted, nearest = False, None
    for node in e2e_bb.path_to_root(g):
        row = row_of.get(node)
        if not row or not row.get("precedence_age") or not row.get("age"):
            continue
        departs = factor(float(row["age"]), float(row["precedence_age"])) > 1.01 if float(row["precedence_age"]) > 0 else False
        if nearest is None and row.get("timetree_age"):
            nearest = departs
        fitted = fitted or departs
    if fitted:
        causes.append("fit_on_path")
    cause_rows.append({"genus": g, "stem_legacy": fmt(main_bb.stem_age(g), 2), "stem_default": fmt(e2e_bb.stem_age(g), 2),
                       "factor": fmt(f), "causes": ";".join(causes), "fit_at_nearest_dated_ancestor": nearest})
write(OUT / "stem_changes.tsv", cause_rows)
n_changed = len(cause_rows)
with_cause = sum(1 for r in cause_rows if r["causes"])
strict = sum(1 for r in cause_rows if r["causes"].replace("fit_on_path", "").strip(";") or r["fit_at_nearest_dated_ancestor"])
share = with_cause / n_changed if n_changed else 1.0
gate("6b", "stem ages that differ by more than 10%: share with a recorded cause", "",
     f"{share:.3f} ({with_cause} of {n_changed})", ">= 0.90", share >= 0.90,
     f"counting a fitted age only at the nearest dated ancestor: {strict} of {n_changed}; "
     f"causes: {dict(Counter(c for r in cause_rows for c in r['causes'].split(';') if c))}")
l, d = statistics.mean(pd_rho["legacy"]), statistics.mean(pd_rho["default"])
gate("6c", "backbone PD of genus samples: rank correlation with and without perturbed ages", fmt(l, 4), fmt(d, 4),
     ">= legacy - 0.01", d >= l - 0.01)
gate("6d", "PD ratio default over legacy, and total tree length", fmt(old.total, 0),
     f"{fmt(new.total, 0)} ({new.total / old.total - 1:+.1%}); "
     + "; ".join(f"{r['design']}: {r['median_ratio']}" for r in pd_rows), "none", None)

# The two full runs side by side (report only)
prov = {arm: {r["genus"]: r for r in read(root / "provenance.tsv")} for arm, root in (("legacy", MAIN), ("default", BRANCH))}
placed = {arm: Counter(r["status"] for r in read(root / "units_placed.tsv")) for arm, root in (("legacy", MAIN), ("default", BRANCH))}
decisions = {"legacy": Counter(main_dec.values()),
             "default": Counter(r["decision"] for r in read(BRANCH / "units/screened.tsv"))}
both = [g for g in prov["legacy"] if g in prov["default"] and prov["legacy"][g]["status"] == "subtree"
        and prov["default"][g]["status"] == "subtree"]
crown_factor = [factor(float(prov["legacy"][g]["crown_age"]), float(prov["default"][g]["crown_age"])) for g in both]
side = []


def add(quantity, legacy, default):
    side.append({"quantity": quantity, "legacy": legacy, "default": default})


add("genera in the backbone", len(main_bb.leaf), len(e2e_bb.leaf))
add("backbone nodes", len(main_bb.by_tips), len(e2e_bb.by_tips))
for arm_stat in ("subtree", "single_unit"):
    add(f"genera: {arm_stat}", *(sum(1 for r in prov[arm].values() if r["status"] == arm_stat) for arm in ("legacy", "default")))
add("units set aside by screening", *(decisions[arm]["other_genus"] + decisions[arm]["isolated"] for arm in ("legacy", "default")))
add("units in genus subtrees", *(sum(int(r["n_units"]) for r in prov[arm].values() if r["status"] == "subtree")
                                 for arm in ("legacy", "default")))
add("genera with units on the wrong side of the root",
    *(sum(1 for r in prov[arm].values() if r.get("violations") and int(r["violations"]) > 0) for arm in ("legacy", "default")))
add("units of a genus among its outgroup",
    *(sum(len(r["misplaced"].split(";")) for r in prov[arm].values() if r.get("misplaced")) for arm in ("legacy", "default")))
add("crown age from TimeTree", *(sum(1 for r in prov[arm].values() if r.get("crown_source") == "timetree") for arm in ("legacy", "default")))
add("TimeTree crown demoted to a bound", "", sum(1 for r in prov["default"].values() if r.get("crown_source") == "midpoint_timetree_demoted"))
add("sum of crown ages, genera with a subtree in both (Ma)", *(fmt(sum(float(prov[arm][g]["crown_age"]) for g in both), 0)
                                                               for arm in ("legacy", "default")))
add("median crown age (Ma)", *(fmt(statistics.median(float(prov[arm][g]["crown_age"]) for g in both), 2) for arm in ("legacy", "default")))
add("crown ages that differ by more than 10%", "", sum(f > 1.10 for f in crown_factor))
add("crown ages that differ by more than a factor 1.5", "", sum(f > 1.5 for f in crown_factor))
add("genus stems that differ by more than 10%", "", len(cause_rows))
add("genus stems that differ by more than a factor 1.5", "", sum(float(r["factor"]) > 1.5 for r in cause_rows))
add("units placed next to their nearest unit", placed["legacy"]["placed"], placed["default"]["placed"])
add("units on the stem of their genus", placed["legacy"]["genus_stem"], placed["default"]["genus_stem"])
add("units not placed", placed["legacy"]["unplaced"], placed["default"]["unplaced"])
add("tips of the final tree", len(old.tip), len(new.tip))
add("total branch length of the final tree (Ma)", fmt(old.total, 0), fmt(new.total, 0))
add("units flagged numt_like", "", sum(1 for r in read(BRANCH / "units/screened.tsv") if r.get("translation") == "numt_like"))
add("units listed as long branches", "", sum(len(r["rate_outliers"].split(";")) for r in prov["default"].values() if r.get("rate_outliers")))
write(OUT / "end_to_end.tsv", side)

# Ablation: one feature switched back at a time, against the default backbone (report only)
ablation = []
d_crowns = {r["genus"]: r["crown_age"] for r in read(EXP / "default_clean/crown_ages.tsv")}
for name in ("identifiers_from_ott", "no_reanchoring", "mrca_anchoring", "precedence_ages", "global_crowns"):
    folder = EXP / "ablate" / name
    if not (folder / "metadata.json").exists():
        continue
    bb = Backbone(folder / "backbone.tre")
    shared_genera = [g for g in bb.leaf if g in default_bb.leaf]
    f = [factor(bb.stem_age(g), default_bb.stem_age(g)) for g in shared_genera]
    crowns = {r["genus"]: r["crown_age"] for r in read(folder / "crown_ages.tsv")}
    meta = json.load(open(folder / "metadata.json"))["ages"]["final"]
    ablation.append({"switched_back": name, "genera": len(bb.leaf),
                     "clades_not_in_default": len(set(bb.node_ages()) - set(default_bb.node_ages())),
                     "stems_moved_over_10pct": sum(x > 1.10 for x in f), "stems_moved_over_1.5": sum(x > 1.5 for x in f),
                     "calibrations_off_by_1.5": meta["severe"], "weighted_rms": meta["weighted_rms"],
                     "timetree_crowns_changed": sum(1 for g in set(crowns) | set(d_crowns)
                                                    if crowns.get(g, "") != d_crowns.get(g, ""))})
write(OUT / "ablation.tsv", ablation)

# 7. Operations ------------------------------------------------------------------------------------------
lrun, drun = json.load(open(BENCH / "main_run.json")), json.load(open(BENCH / "branch_run.json"))


def rule_seconds(root):
    out = {}
    for path in (root / ".snakemake/metadata").iterdir():
        try:
            meta = json.load(open(path))
        except (OSError, ValueError):
            continue
        if meta.get("rule") and meta.get("endtime"):
            out.setdefault(meta["rule"], set()).add(round(meta["endtime"] - meta["starttime"], 1))
    return {rule: round(sum(times), 1) for rule, times in out.items()}


lsec, dsec = rule_seconds(BENCH / "main"), rule_seconds(BENCH / "branch")
write(OUT / "runtime.tsv", [{"rule": rule, "legacy_seconds": lsec.get(rule, ""), "default_seconds": dsec.get(rule, "")}
                            for rule in sorted(set(lsec) | set(dsec))]
      + [{"rule": "whole run (wall)", "legacy_seconds": lrun["seconds"], "default_seconds": drun["seconds"]},
         {"rule": "peak memory (MB)", "legacy_seconds": lrun["peak_rss_mb"], "default_seconds": drun["peak_rss_mb"]}])
ratio = drun["seconds"] / lrun["seconds"]
gate("7a", "wall time of the full run (minutes)", fmt(lrun["seconds"] / 60, 1), fmt(drun["seconds"] / 60, 1),
     "<= 1.5 x legacy", ratio <= 1.5 and drun["returncode"] == 0 and lrun["returncode"] == 0, f"ratio {ratio:.2f}")
gate("7b", "peak memory of the full run (MB)", lrun["peak_rss_mb"], drun["peak_rss_mb"], "<= 1.5 x legacy and < 7000",
     drun["peak_rss_mb"] <= 1.5 * lrun["peak_rss_mb"] and drun["peak_rss_mb"] < 7000)
FILES = ("backbone.tre", "node_ages.tsv", "genera.tsv", "crown_ages.tsv", "concordance.tsv")
again = {f: filecmp.cmp(EXP / "default_clean" / f, EXP / "repro/again" / f, False) for f in FILES}
workflow = {f: filecmp.cmp(EXP / "default_clean" / f, BRANCH / "backbone" / f, False) for f in FILES}
screen_same = filecmp.cmp(EXP / "default_units_again/screened.tsv", BRANCH / "units/screened.tsv", False)
gate("7c", "backbone and screening rebuilt: byte-identical", "",
     f"backbone again: {sum(again.values())} of {len(FILES)} files; harness against workflow: "
     f"{sum(workflow.values())} of {len(FILES)}; screening: {'identical' if screen_same else 'differs'}", "identical",
     all(again.values()) and all(workflow.values()) and screen_same,
     "; ".join(f for f, same in {**again, **{k + ' (workflow)': v for k, v in workflow.items()}}.items() if not same))
offline_ok = (EXP / "repro/offline/metadata.json").exists()
offline = {f: offline_ok and filecmp.cmp(EXP / "default_clean" / f, EXP / "repro/offline" / f, False) for f in FILES}
gate("7d", "backbone rebuilt with the network blocked", "",
     f"{'built' if offline_ok else 'failed'}; {sum(offline.values())} of {len(FILES)} files identical",
     "succeeds, identical", offline_ok and all(offline.values()))
tests = (BENCH / "tests.txt").read_text().strip() if (BENCH / "tests.txt").exists() else "not run"
gate("7e", "test suite", "", tests, "passes", "passed" in tests and "failed" not in tests and "error" not in tests)

write(OUT / "metrics.tsv", metrics)
must = [m for m in metrics if m["status"] != "report"]
failed = [m["metric"] for m in must if m["status"] == "FAIL"]
print(f"\n{len(must) - len(failed)} of {len(must)} must metrics pass; failed: {failed or 'none'}")
