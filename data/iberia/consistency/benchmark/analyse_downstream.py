"""The downstream runs of #81 in more detail than the gate metrics (report only).

    python analyse_downstream.py EXP INJECTION_GENERA OUT

EXP is ``exp/`` with ``downstream/`` as ``run_downstream.sh`` leaves it; INJECTION_GENERA is the
``injection_genera.tsv`` that ``evaluate.py`` wrote for the same benchmark (which injected genera were flagged and
moved back). OUT gets:

- ``downstream_runs.tsv``: per run, against the clean run of the same setting: genera with a crown age in both,
  how far crown ages move (median, 90% quantile, share above 1.1 and 1.5), tips and total length of the final
  tree, and units placed;
- ``downstream_injected.tsv``: per injected genus of replicates 1 and 2, its crown age in the clean run and with
  the error, for both settings, and whether the default flagged and corrected it;
- ``downstream_collateral.tsv``: per injection run, the genera without an injected error whose crown age moves by
  more than 10%, and how many of those have a stem that moved too.
"""

import statistics
import sys
from collections import Counter
from pathlib import Path

from benchlib import PDTree, factor, quantile, read, write

EXP, GENERA, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
D = EXP / "downstream"
SEEDS, ARMS = (1, 2), ("legacy", "default")


def provenance(folder):
    return {r["genus"]: r for r in read(folder / "provenance.tsv")}


def crowns(rows):
    return {g: float(r["crown_age"]) for g, r in rows.items() if r.get("status") == "subtree" and r.get("crown_age")}


def fmt(x, digits=3):
    return "" if x is None else f"{x:.{digits}f}"


clean = {arm: provenance(D / "clean" / arm) for arm in ARMS}
clean_crowns = {arm: crowns(clean[arm]) for arm in ARMS}
clean_length = {arm: PDTree(D / "clean" / arm / "timegraft.tre").total for arm in ARMS}
runs = []
for kind in ("inject", "noise", "sub"):
    for seed in SEEDS:
        for arm in ARMS:
            folder = D / kind / f"s{seed}" / arm
            rows = provenance(folder)
            c = crowns(rows)
            f = [factor(c[g], clean_crowns[arm][g]) for g in c if g in clean_crowns[arm] and c[g] > 0
                 and clean_crowns[arm][g] > 0]
            tree = PDTree(folder / "timegraft.tre")
            placed = Counter(r["status"] for r in read(folder / "units_placed.tsv"))
            runs.append({"experiment": kind, "seed": seed, "arm": arm, "genera_compared": len(f),
                         "crown_median_factor": fmt(statistics.median(f)), "crown_p90_factor": fmt(quantile(f, 0.9)),
                         "crown_share_over_1.1": fmt(sum(x > 1.1 for x in f) / len(f)),
                         "crown_share_over_1.5": fmt(sum(x > 1.5 for x in f) / len(f)),
                         "tips": len(tree.tip), "total_length_ratio": fmt(tree.total / clean_length[arm], 4),
                         "units_placed": placed.get("placed", 0), "units_on_genus_stem": placed.get("genus_stem", 0)})
write(OUT / "downstream_runs.tsv", runs)

injected = [r for r in read(GENERA) if int(r["seed"]) in SEEDS]
rows_out, collateral = [], []
for seed in SEEDS:
    mine = {r["genus"]: r for r in injected if int(r["seed"]) == seed}
    for arm in ARMS:
        run = provenance(D / "inject" / f"s{seed}" / arm)
        c = crowns(run)
        stems_moved = 0
        moved = []
        for g in sorted(set(c) & set(clean_crowns[arm]) - set(mine)):
            if factor(c[g], clean_crowns[arm][g]) > 1.10:
                moved.append(g)
                a, b = run[g].get("stem_age"), clean[arm][g].get("stem_age")
                stems_moved += bool(a and b and factor(float(a), float(b)) > 1.10)
        collateral.append({"seed": seed, "arm": arm, "other_genera_compared": len(set(c) & set(clean_crowns[arm]) - set(mine)),
                           "other_crowns_moved_over_10pct": len(moved), "of_which_stem_moved_over_10pct": stems_moved,
                           "genera": ";".join(moved[:30]) + (";..." if len(moved) > 30 else "")})
    for g, r in sorted(mine.items()):
        row = {"seed": seed, "genus": g, "kind": r["kind"], "flagged": r["concordance_flag"] == "flag",
               "corrected": r["concordance_flag"] == "flag" and r["default_own_class"] == "True"}
        for arm in ARMS:
            c = crowns(provenance(D / "inject" / f"s{seed}" / arm))
            before, after = clean_crowns[arm].get(g), c.get(g)
            row[f"{arm}_crown_clean"] = fmt(before, 2)
            row[f"{arm}_crown_injected"] = fmt(after, 2)
            row[f"{arm}_crown_factor"] = fmt(factor(after, before)) if before and after else ""
        rows_out.append(row)
write(OUT / "downstream_injected.tsv", rows_out)
write(OUT / "downstream_collateral.tsv", collateral)
for arm in ARMS:
    f = [float(r[f"{arm}_crown_factor"]) for r in rows_out if r[f"{arm}_crown_factor"]]
    print(f"{arm}: injected genera with a crown in both {len(f)} of {len(rows_out)}; within 1.5: "
          f"{sum(x <= 1.5 for x in f)}; median factor {fmt(statistics.median(f)) if f else ''}")
for r in collateral:
    print(r["seed"], r["arm"], r["other_crowns_moved_over_10pct"], "of", r["other_genera_compared"],
          "other crowns moved; stems moved too:", r["of_which_stem_moved_over_10pct"])
