"""Detection metrics of the gate (#56) with and without the lineage-alone rule (#71).

    python evaluate_lineage_alone.py EXP OUT

EXP is ``exp/`` as ``run_lineage_alone.sh`` leaves it. Two arms, the default (``lineage_alone_level: class``)
and the control (``off``, the rule of the gate), each on the clean run and on both injection variants:
``merge`` (the gate's injection: OTT merges the genus, with its NCBI ids, into a taxon elsewhere) and ``name``
(only the name moves). Per variant and arm, with the definitions of ``evaluate.py``:

- flagged (1b, 1c), by rule, and listed for review;
- placed back in their own class, and placed in another class (1d: wrongly);
- genera without an injected error that are flagged (1e; *Onchidella* not counted);
- stems of the corrected genera within a factor 1.5 of the clean build of the same arm;
- stems of other genera that move by more than 10% (as in 1g), against the clean build of the same arm.

OUT gets ``lineage_alone.tsv`` (one row per variant and arm), ``lineage_alone_genera.tsv`` (every injected
genus, both arms) and ``lineage_alone_clean.tsv`` (the genera flagged in the clean run, both arms).
"""

import sys
from collections import Counter
from pathlib import Path

from benchlib import Backbone, factor, read, write

EXP, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
SEEDS = [1, 2, 3, 4, 5]
ARMS = {"default": "lineage_alone_level: class", "control": "lineage_alone_level: off"}
klass = {r["genus"]: r["records_class"] for r in read(EXP / "default_clean/genera.tsv")}

summary, genera, clean_rows = [], [], []
for arm in ARMS:
    conc = read(EXP / f"{arm}_clean/concordance.tsv")
    clean_rows += [{"arm": arm, "genus": r["genus"], "flag_rule": r["flag_rule"], "severity": r["severity"],
                    "n_flags": r["n_flags"], "lineage_distance": r["lineage_distance"]}
                   for r in conc if r["concordance_flag"] == "flag"]
for variant in ("merge", "name"):
    for arm, setting in ARMS.items():
        clean = Backbone(EXP / f"{arm}_clean/backbone.tre")
        rows, false_flags, moved = [], [], []
        for seed in SEEDS:
            manifest = read(EXP / f"{variant}/s{seed}/manifest.tsv")
            injected = {m["genus"] for m in manifest}
            folder = EXP / f"{variant}/s{seed}/{arm}"
            conc = {r["genus"]: r for r in read(folder / "concordance.tsv")}
            bb = Backbone(folder / "backbone.tre")
            for m in manifest:
                g = m["genus"]
                c = conc.get(g, {})
                own = None
                if g in bb.leaf:
                    others = bb.smallest_clade_with_company(g)
                    own = bool(others) and all(klass.get(o) == klass[g] for o in others)
                f = factor(bb.stem_age(g), clean.stem_age(g)) if g in bb.leaf and g in clean.leaf else None
                rows.append({"variant": variant, "arm": arm, "seed": seed, **m,
                             "lineage_distance": c.get("lineage_distance", ""), "n_flags": c.get("n_flags", ""),
                             "concordance_flag": c.get("concordance_flag", ""), "flag_rule": c.get("flag_rule", ""),
                             "placement_action": c.get("placement_action", ""), "own_class": own,
                             "stem_factor": "" if f is None else f"{f:.3f}"})
            moved.append(sum(1 for g in bb.leaf if g in clean.leaf and g not in injected
                             and factor(bb.stem_age(g), clean.stem_age(g)) > 1.10))
            false_flags.append(sorted(g for g, c in conc.items() if c["concordance_flag"] == "flag"
                                      and g not in injected and g != "Onchidella"))
        flagged = [r for r in rows if r["concordance_flag"] == "flag"]
        corrected = [r for r in flagged if r["own_class"]]
        near = [r for r in corrected if r["stem_factor"] and float(r["stem_factor"]) <= 1.5]
        rules = Counter(r["flag_rule"] for r in flagged)
        summary.append({
            "variant": variant, "arm": arm, "setting": setting, "injected": len(rows), "flagged": len(flagged),
            "flagged_share": f"{len(flagged) / len(rows):.3f}",
            "by_two_checks": rules.get("two_checks", 0), "by_lineage_alone": rules.get("lineage_alone", 0),
            "review": sum(1 for r in rows if r["concordance_flag"] == "review"),
            "back_in_own_class": len(corrected),
            "in_another_class": sum(1 for r in flagged if r["placement_action"] == "reanchored" and r["own_class"] is False),
            "false_flags": sum(map(len, false_flags)),
            "false_flags_per_seed": ";".join(str(len(f)) for f in false_flags),
            "false_flag_genera": ";".join(g for f in false_flags for g in f),
            "corrected_stem_within_1.5": len(near),
            "other_stems_moved": sum(moved), "other_stems_moved_per_seed": ";".join(map(str, moved)),
        })
        genera += rows
for row in summary:
    print(row)
write(OUT / "lineage_alone.tsv", summary)
write(OUT / "lineage_alone_genera.tsv", genera)
write(OUT / "lineage_alone_clean.tsv", clean_rows)
