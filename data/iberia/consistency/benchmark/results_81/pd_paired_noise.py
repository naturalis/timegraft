"""How much the paired PD metrics of gate version 3 vary with the sample seed (#104).

    python pd_paired_noise.py DOWNSTREAM OUT DRAWS FIRST LAST

DOWNSTREAM is ``exp/downstream/`` as ``run_downstream.sh`` leaves it. For each sample seed from FIRST to LAST,
computes 1o, 6e and 6f as gate version 3 does (``benchlib.paired_pd``: per replicate, the same DRAWS samples per
design for both settings, drawn from the units in all trees compared), and writes per experiment, seed and
design the value for each setting and the difference default minus legacy. The last rows give, per experiment
and design, the mean and standard deviation of the difference across seeds, and in how many seeds it is below
-0.01 (where the metric would fail).
"""

import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchlib import DESIGNS, PDTree, paired_pd, read, write  # noqa: E402

D, OUT, DRAWS, FIRST, LAST = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
ARMS = ("legacy", "default")
unit_genus = {}
for arm in ARMS:
    unit_genus.update({r["unit_id"]: r["genus"] for r in read(D / "clean" / arm / "units/representatives.tsv")})
clean = {arm: PDTree(D / "clean" / arm / "timegraft.tre") for arm in ARMS}
runs = {kind: [{arm: PDTree(D / kind / f"s{s}" / arm / "timegraft.tre") for arm in ARMS} for s in (1, 2)]
        for kind in ("inject", "noise", "sub")}
metric = {"inject": "1o", "noise": "6e", "sub": "6f"}
rows, diffs = [], {(k, d): [] for k in runs for d in DESIGNS}
for seed in range(FIRST, LAST + 1):
    for kind, reps in runs.items():
        value = paired_pd(clean, reps, unit_genus, seed, DRAWS)
        for d in DESIGNS:
            diff = value["default"][d] - value["legacy"][d]
            diffs[(kind, d)].append(diff)
            rows.append({"metric": metric[kind], "seed": seed, "design": d, "legacy": f"{value['legacy'][d]:.4f}",
                         "default": f"{value['default'][d]:.4f}", "difference": f"{diff:.4f}"})
    print(seed, {metric[k]: [r["difference"] for r in rows[-9:] if r["metric"] == metric[k]] for k in runs}, flush=True)
for (kind, d), v in diffs.items():
    rows.append({"metric": metric[kind], "seed": "all", "design": d, "legacy": "", "default": "",
                 "difference": f"mean {statistics.mean(v):.4f}, sd {statistics.stdev(v):.4f}, "
                               f"below -0.01 in {sum(x < -0.01 for x in v)} of {len(v)}"})
write(OUT, rows)
for r in rows[-9:]:
    print(r["metric"], r["design"], r["difference"])
