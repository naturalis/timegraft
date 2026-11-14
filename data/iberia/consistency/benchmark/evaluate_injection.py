"""Metrics 1d, 1f and 1g of the gate (#56) for the injected merges, with and without the identifier check (#69).

    python evaluate_injection.py EXP OUT

EXP is ``exp/`` as ``run_injection.sh`` leaves it. Three arms, each against its own clean build, with the
logic of ``injection()`` in ``evaluate.py``:

- legacy: ``inject/sN/legacy`` against ``legacy_clean``;
- default: ``inject/sN/default`` against ``default_clean`` (``reanchor.check_ids: true``);
- control: ``inject_control/sN/default`` against ``default_clean_control`` (the same code, ``check_ids: false``),
  which reproduces the default arm of the gate.

OUT gets ``injection_69.tsv`` (one row per arm), ``injection_69_genera.tsv`` (every injected genus in the
default arm, with the ids dropped) and ``injection_69_collateral.tsv``: why the stem of each other genus moves
by more than 10% in the default arm, from its stem node in ``node_ages.tsv`` against the clean build.
"""

import sys
from collections import Counter
from pathlib import Path

from benchlib import Backbone, factor, read, write

EXP, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
SEEDS = [1, 2, 3, 4, 5]
klass = {r["genus"]: r["records_class"] for r in read(EXP / "default_clean/genera.tsv")}


def load(folder):
    """The backbone and the node_ages.tsv row of every internal node (node ids number all nodes in preorder)."""
    bb = Backbone(folder / "backbone.tre")
    rows = {r["node"]: r for r in read(folder / "node_ages.tsv")}
    by = {n: rows[f"n{i}"] for i, n in enumerate(bb.tree.preorder_node_iter()) if not n.is_leaf()}
    assert len(by) == len(rows), folder
    return bb, by


def arm(folder_of, clean_folder):
    clean = Backbone(clean_folder / "backbone.tre")
    rows, collateral = [], []
    for seed in SEEDS:
        manifest = read(EXP / f"inject/s{seed}/manifest.tsv")
        injected = {m["genus"] for m in manifest}
        folder = folder_of(seed)
        conc = {r["genus"]: r for r in read(folder / "concordance.tsv")}
        gen = {r["genus"]: r for r in read(folder / "genera.tsv")}
        bb = Backbone(folder / "backbone.tre")
        for m in manifest:
            g = m["genus"]
            c = conc.get(g, {})
            own = None
            if g in bb.leaf:
                others = bb.smallest_clade_with_company(g)
                own = bool(others) and all(klass.get(o) == klass[g] for o in others)
            f = factor(bb.stem_age(g), clean.stem_age(g)) if g in bb.leaf and g in clean.leaf else None
            rows.append({**m, "concordance_flag": c.get("concordance_flag", ""),
                         "placement_action": c.get("placement_action", ""), "own_class": own,
                         "stem_factor": "" if f is None else f"{f:.3f}",
                         "reanchor_ids_checked": gen.get(g, {}).get("reanchor_ids_checked", ""),
                         "reanchor_ids_dropped": gen.get(g, {}).get("reanchor_ids_dropped", ""),
                         "reanchor_dropped_ids": c.get("reanchor_dropped_ids", "")})
        shared = [g for g in bb.leaf if g in clean.leaf and g not in injected]
        collateral.append(sum(1 for g in shared if factor(bb.stem_age(g), clean.stem_age(g)) > 1.10))
    return rows, collateral


def summary(name, rows, collateral):
    flagged = [r for r in rows if r["concordance_flag"] == "flag"]
    corrected = [r for r in flagged if r["own_class"]]
    near = [r for r in corrected if r["stem_factor"] and float(r["stem_factor"]) <= 1.5]
    return {"arm": name, "injected": len(rows), "flagged": len(flagged), "back_in_own_class": len(corrected),
            "stem_within_1.5": len(near), "share_within_1.5": f"{len(near) / len(corrected):.3f}" if corrected else "",
            "other_stems_moved": sum(collateral), "other_stems_moved_per_seed": ";".join(map(str, collateral))}


def causes():
    """Why each other stem moves by more than 10% in the default arm."""
    cbb, crow = load(EXP / "default_clean")
    out = []
    for seed in SEEDS:
        manifest = {m["genus"] for m in read(EXP / f"inject/s{seed}/manifest.tsv")}
        bb, row = load(EXP / f"inject/s{seed}/default")
        gen = {r["genus"]: r for r in read(EXP / f"inject/s{seed}/default/genera.tsv")}
        moved = {g for g in manifest if gen.get(g, {}).get("placement_action") == "reanchored"}
        for g in sorted(bb.leaf):
            if g in manifest or g not in cbb.leaf or factor(bb.stem_age(g), cbb.stem_age(g)) <= 1.10:
                continue
            node, cnode = bb.parent[bb.leaf[g]], cbb.parent[cbb.leaf[g]]
            a, b = row[node], crow[cnode]
            if bb.tips[node] - manifest != cbb.tips[cnode] - manifest:
                cause = "stem clade differs"
            elif a["query_ids"] != b["query_ids"]:
                cause = "other TimeTree query"
            elif a["source"].startswith("interpolated"):
                cause = "interpolated"
            else:
                cause = "same calibration, refitted"
            # the nearest dated ancestor, and the re-anchored injected genera in its clade
            anc = node
            while anc in bb.parent and row[anc]["source"].startswith("interpolated"):
                anc = bb.parent[anc]
            out.append({"seed": seed, "genus": g, "factor": f"{factor(bb.stem_age(g), cbb.stem_age(g)):.3f}",
                        "cause": cause, "source": a["source"], "clean_source": b["source"],
                        "dated_ancestor": row[anc]["node"], "dated_ancestor_genera": len(bb.tips[anc]),
                        "dated_ancestor_age": f"{bb.age[anc]:.1f}",
                        "reanchored_in_dated_ancestor": ";".join(sorted(set(bb.tips[anc]) & moved))})
    return out


results = []
rows, collateral = arm(lambda s: EXP / f"inject/s{s}/legacy", EXP / "legacy_clean")
results.append(summary("legacy", rows, collateral))
rows, collateral = arm(lambda s: EXP / f"inject/s{s}/default", EXP / "default_clean")
write(OUT / "injection_69_genera.tsv", rows)
results.append(summary("default (check_ids: true)", rows, collateral))
rows, collateral = arm(lambda s: EXP / f"inject_control/s{s}/default", EXP / "default_clean_control")
results.append(summary("control (check_ids: false)", rows, collateral))
legacy_moved = results[0]["other_stems_moved"]
for r in results:
    r["other_stems_moved_vs_legacy"] = f"{r['other_stems_moved'] / legacy_moved:.3f}" if legacy_moved else ""
    print(r)
write(OUT / "injection_69.tsv", results)
collat = causes()
write(OUT / "injection_69_collateral.tsv", collat)
print(dict(Counter(c["cause"] for c in collat)))
print("in a dated clade that holds a re-anchored genus:", sum(1 for c in collat if c["reanchored_in_dated_ancestor"]))
