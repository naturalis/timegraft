"""Why metric 6f fails (#81): PD rank correlation under subsampling with and without the placed units.

    python pd_without_placed.py DOWNSTREAM OUT

DOWNSTREAM is ``exp/downstream/`` as ``run_downstream.sh`` leaves it; run from the benchmark folder (it imports
``benchlib``). For each setting, the PD of the three sample designs of metric 6a (seed 6) on the clean final tree
and on the final tree of each subsampling replicate (1 and 2), as ``evaluate.py`` computes 6f, and again with the
samples drawn only from units that neither run placed by nearest neighbour or on their genus's stem
(``units_placed.tsv`` status ``placed`` or ``genus_stem``). Writes OUT as a table of mean rank correlations.
"""

import csv
import random
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchlib import PDTree, read, spearman  # noqa: E402

D, OUT = Path(sys.argv[1]), Path(sys.argv[2])


def designs(ta, tb, unit_genus, exclude):
    both = sorted((set(ta.tip) & set(tb.tip)) - exclude)
    by_genus = defaultdict(list)
    for u in both:
        if u in unit_genus:
            by_genus[unit_genus[u]].append(u)
    r = random.Random(6)
    names, many = sorted(by_genus), sorted(g for g, us in by_genus.items() if len(us) >= 5)
    draws = {"150 units at random": [r.sample(both, 150) for _ in range(200)],
             "100 genera, one unit each": [[r.choice(by_genus[g]) for g in r.sample(names, 100)] for _ in range(200)],
             "100 genera, five units each": [[u for g in r.sample(many, 100) for u in r.sample(by_genus[g], 5)]
                                             for _ in range(200)]}
    return {k: spearman([ta.pd(x) for x in s], [tb.pd(x) for x in s]) for k, s in draws.items()}


rows = []
for arm in ("legacy", "default"):
    unit_genus = {r["unit_id"]: r["genus"] for r in read(D / "clean" / arm / "units/representatives.tsv")}
    clean = PDTree(D / "clean" / arm / "timegraft.tre")
    for units in ("all units", "without placed units"):
        rhos = defaultdict(list)
        for seed in (1, 2):
            run = D / "sub" / f"s{seed}" / arm
            placed = set()
            if units != "all units":
                for path in (D / "clean" / arm / "units_placed.tsv", run / "units_placed.tsv"):
                    placed |= {r["unit_id"] for r in read(path) if r["status"] in ("placed", "genus_stem")}
            for design, rho in designs(clean, PDTree(run / "timegraft.tre"), unit_genus, placed).items():
                rhos[design].append(rho)
        rows += [{"arm": arm, "units": units, "design": d, "pd_rank_correlation": f"{statistics.mean(v):.4f}"}
                 for d, v in rhos.items()]
with open(OUT, "w", newline="") as handle:
    w = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
for r in rows:
    print(r)
