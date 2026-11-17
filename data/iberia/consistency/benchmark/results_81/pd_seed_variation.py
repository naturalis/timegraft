"""How much metric 6f varies with the seed of its samples alone (#81).

    python pd_seed_variation.py DOWNSTREAM OUT [FIRST LAST]

DOWNSTREAM is ``exp/downstream/`` as ``run_downstream.sh`` leaves it. For each sample seed from FIRST to LAST
(default 6 to 25; the gate uses 6), computes 6f exactly as ``evaluate.py`` does (the three designs of 6a, 200
samples each, drawn per setting from the units in its clean and subsampled final trees, mean over the two
subsampling replicates), and writes per seed and design the value for each setting, the difference default minus
legacy, and whether 6f would pass (the difference at least -0.01). The final rows summarise each design: mean and
standard deviation of the difference, and the number of seeds where it fails; and the number of seeds where at
least one design fails, which is when 6f as declared fails.
"""

import random
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchlib import PDTree, read, spearman, write  # noqa: E402

D, OUT = Path(sys.argv[1]), Path(sys.argv[2])
FIRST, LAST = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (6, 25)
DESIGNS = ("150 units at random", "100 genera, one unit each", "100 genera, five units each")


def designs(ta, tb, unit_genus, seed):
    """The three designs of 6a with ``seed``, as evaluate.py draws them for 6f."""
    both = sorted(set(ta.tip) & set(tb.tip))
    by_genus = defaultdict(list)
    for u in both:
        if u in unit_genus:
            by_genus[unit_genus[u]].append(u)
    r = random.Random(seed)
    names, many = sorted(by_genus), sorted(g for g, us in by_genus.items() if len(us) >= 5)
    draws = [[r.sample(both, 150) for _ in range(200)],
             [[r.choice(by_genus[g]) for g in r.sample(names, 100)] for _ in range(200)],
             [[u for g in r.sample(many, 100) for u in r.sample(by_genus[g], 5)] for _ in range(200)]]
    return {d: spearman([ta.pd(x) for x in s], [tb.pd(x) for x in s]) for d, s in zip(DESIGNS, draws)}


trees = {}
for arm in ("legacy", "default"):
    unit_genus = {r["unit_id"]: r["genus"] for r in read(D / "clean" / arm / "units/representatives.tsv")}
    trees[arm] = (unit_genus, PDTree(D / "clean" / arm / "timegraft.tre"),
                  [PDTree(D / "sub" / f"s{s}" / arm / "timegraft.tre") for s in (1, 2)])
rows, diffs, any_fail = [], defaultdict(list), 0
for seed in range(FIRST, LAST + 1):
    value = {}
    for arm, (unit_genus, clean, subs) in trees.items():
        per = defaultdict(list)
        for sub in subs:
            for design, rho in designs(clean, sub, unit_genus, seed).items():
                per[design].append(rho)
        value[arm] = {d: statistics.mean(v) for d, v in per.items()}
    failed = False
    for d in DESIGNS:
        diff = value["default"][d] - value["legacy"][d]
        diffs[d].append(diff)
        failed |= diff < -0.01
        rows.append({"seed": seed, "design": d, "legacy": f"{value['legacy'][d]:.4f}",
                     "default": f"{value['default'][d]:.4f}", "difference": f"{diff:.4f}",
                     "passes": "yes" if diff >= -0.01 else "no"})
    any_fail += failed
    print(seed, {d: rows[-3 + i]["difference"] for i, d in enumerate(DESIGNS)}, flush=True)
for d in DESIGNS:
    v = diffs[d]
    rows.append({"seed": "all", "design": d, "legacy": "", "default": "",
                 "difference": f"mean {statistics.mean(v):.4f}, sd {statistics.stdev(v):.4f}",
                 "passes": f"fails in {sum(x < -0.01 for x in v)} of {len(v)} seeds"})
rows.append({"seed": "all", "design": "any design (6f as declared)", "legacy": "", "default": "", "difference": "",
             "passes": f"fails in {any_fail} of {len(diffs[DESIGNS[0]])} seeds"})
write(OUT, rows)
for r in rows[-4:]:
    print(r)
