"""The options left at their conservative defaults, each against the default, end to end (#73).

    PYTHONPATH=src python data/iberia/consistency/options/compare.py OUT DEFAULT NAME=RESULTS [NAME=RESULTS ...]

DEFAULT and each RESULTS are result folders (``results/iberia``) of full runs of the workflow on
``data/iberia/records.tsv.gz`` with ``config/config.yaml`` and ``config/iberia.yaml``; each alternative adds an
overlay that changes one setting. Writes to OUT:

- ``options.tsv``: one row per run, the default first, with the measures of #73: rooting (genera with units on
  the wrong side of the root, units among the outgroup, crowns that depend strongly on the root, as in #75),
  crown ages (count and size of changes against the default; the tails of gate metrics 3a to 3d), PD of the
  gate's three sample designs on the final trees against the default (metric 6a), and units set aside and
  placed;
- ``crown_changes.tsv``: every genus whose crown age differs by more than 10% from the default, per run.

Crown ages are compared for genera with a subtree in both runs. The PD designs and their seed (6) are those
of ``benchmark/evaluate.py``, drawn from the units in both trees.
"""

import csv
import random
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmark"))
from benchlib import PDTree, quantile, spearman  # noqa: E402

csv.field_size_limit(sys.maxsize)
SET_ASIDE = ("other_genus", "isolated", "numt_like")
ROOT_FACTOR = 1.25  # crown_age_max / crown_age_min at which a crown depends strongly on the root (#75)


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows):
    columns = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as handle:
        out = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n", restval="")
        out.writeheader()
        out.writerows(rows)


class Run:
    def __init__(self, folder):
        folder = Path(folder)
        self.prov = {r["genus"]: r for r in read(folder / "provenance.tsv")}
        self.subtrees = {g: r for g, r in self.prov.items() if r.get("status") == "subtree"}
        self.screened = read(folder / "units/screened.tsv")
        self.placed = read(folder / "units_placed.tsv")
        self.aside_genus = [r for f in sorted((folder / "subtrees").glob("set_aside_*.tsv")) for r in read(f)]
        self.tree = PDTree(folder / "timegraft.tre")
        self.unit_genus = {r["unit_id"]: r["genus"] for r in read(folder / "units/representatives.tsv")}

    def crowns(self):
        return {g: float(r["crown_age"]) for g, r in self.subtrees.items() if r.get("crown_age")}


def tails(run):
    out = Counter()
    for r in run.subtrees.values():
        if not r["crown_age"]:
            continue
        crown = float(r["crown_age"])
        if r.get("barcode_age"):
            ratio = crown / float(r["barcode_age"])
            out["3a"] += ratio > 20
            out["3b"] += ratio > 10
            out["3c"] += r["crown_source"].startswith("timetree") and ratio < 1 / 3
        if r.get("root_to_tip") and float(r["root_to_tip"]) < 0.02 and crown > 20:
            out["3d"] += 1
    return out


def strongly(r):
    lo, hi = r.get("crown_age_min", ""), r.get("crown_age_max", "")
    return bool(lo and hi and float(lo) > 0 and float(hi) / float(lo) >= ROOT_FACTOR)


def pd_designs(a, b, unit_genus):
    """Rank correlation of PD on the two final trees, for the three designs of gate metric 6a."""
    shared = sorted(set(a.tree.tip) & set(b.tree.tip))
    by_genus = defaultdict(list)
    for u in shared:
        if u in unit_genus:
            by_genus[unit_genus[u]].append(u)
    rng = random.Random(6)
    genera, rich = sorted(by_genus), sorted(g for g, us in by_genus.items() if len(us) >= 5)
    designs = {
        "150 units": [rng.sample(shared, 150) for _ in range(200)],
        "100 genera x 1": [[rng.choice(by_genus[g]) for g in rng.sample(genera, 100)] for _ in range(200)],
        "100 genera x 5": [[u for g in rng.sample(rich, 100) for u in rng.sample(by_genus[g], 5)] for _ in range(200)],
    }
    out = {}
    for name, samples in designs.items():
        x, y = [a.tree.pd(s) for s in samples], [b.tree.pd(s) for s in samples]
        ratio = [q / p for p, q in zip(x, y)]
        out[name] = (spearman(x, y), statistics.median(ratio), quantile(ratio, 0.05), quantile(ratio, 0.95))
    return out, len(shared)


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    runs = {"default": Run(sys.argv[2])}
    for arg in sys.argv[3:]:
        name, folder = arg.split("=", 1)
        runs[name] = Run(folder)
    base = runs["default"]
    base_crowns = base.crowns()
    rows, changes = [], []
    for name, run in runs.items():
        crowns = run.crowns()
        both = sorted(set(crowns) & set(base_crowns))
        rel = {g: crowns[g] / base_crowns[g] - 1 for g in both}
        moved = {g: x for g, x in rel.items() if abs(x) > 0.10}
        t = tails(run)
        row = {
            "run": name,
            "units set aside by screening": sum(1 for r in run.screened if r["decision"] in SET_ASIDE),
            "of which numt_like": sum(1 for r in run.screened if r["decision"] == "numt_like"),
            "units set aside from subtrees": len(run.aside_genus),
            "units placed next to their nearest units": sum(1 for r in run.placed if r["status"] == "placed"),
            "units joined to their genus's stem": sum(1 for r in run.placed if r["status"] == "genus_stem"),
            "units left out (unplaced)": sum(1 for r in run.placed if r["status"] == "unplaced"),
            "tips in the final tree": len(run.tree.tip),
            "genus subtrees": len(run.subtrees),
            "genera with units on the wrong side of the root": sum(
                1 for r in run.subtrees.values() if r.get("violations") and int(r["violations"]) > 0),
            "units of a genus among its outgroup": sum(
                len(r["misplaced"].split(";")) for r in run.subtrees.values() if r.get("misplaced")),
            "crowns that depend strongly on the root": sum(1 for r in run.subtrees.values() if strongly(r)),
            "crowns compared": len(both),
            "crowns changed by more than 10%": len(moved),
            "older / younger": f"{sum(1 for x in moved.values() if x > 0)} / {sum(1 for x in moved.values() if x < 0)}",
            "crowns changed by more than a quarter": sum(1 for x in rel.values() if abs(x) > 0.25),
            "median change of the changed crowns": f"{statistics.median(moved.values()):+.3f}" if moved else "",
            "sum of crown ages against the default": f"{sum(crowns[g] for g in both) / sum(base_crowns[g] for g in both) - 1:+.4f}" if both else "",
            "3a crown > 20 x barcode age": t["3a"], "3b crown > 10 x barcode age": t["3b"],
            "3c TimeTree crown < barcode age / 3": t["3c"], "3d root-to-tip < 0.02, crown > 20 Ma": t["3d"],
        }
        if name != "default":
            pd, shared = pd_designs(base, run, base.unit_genus | run.unit_genus)
            row["units in both final trees"] = shared
            for design, (rho, med, lo, hi) in pd.items():
                row[f"6a PD rank correlation, {design}"] = f"{rho:.4f}"
                row[f"PD ratio median [5%, 95%], {design}"] = f"{med:.4f} [{lo:.4f}, {hi:.4f}]"
            for g in sorted(moved, key=lambda g: -abs(moved[g])):
                d, r = base.subtrees[g], run.subtrees[g]
                changes.append({"run": name, "genus": g, "default_crown_age": f"{base_crowns[g]:.4f}",
                                "crown_age": f"{crowns[g]:.4f}", "change": f"{moved[g]:+.3f}",
                                "default_crown_source": d["crown_source"], "crown_source": r["crown_source"],
                                "default_n_units": d["n_units"], "n_units": r["n_units"],
                                "default_misplaced": d.get("misplaced", ""), "misplaced": r.get("misplaced", ""),
                                "rate_outliers": r.get("rate_outliers", ""), "misplaced_set_aside": r.get("misplaced_set_aside", "")})
        rows.append(row)
    write(out / "options.tsv", rows)
    write(out / "crown_changes.tsv", changes)
    keys = list(rows[0])
    for k in dict.fromkeys(k for r in rows for k in r):
        print(f"{k:55}", " | ".join(str(r.get(k, ""))[:22] for r in rows))


if __name__ == "__main__":
    main()
