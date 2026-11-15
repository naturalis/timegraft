"""The relative-rate figures of the screening section, from two genus-stage runs (#55, #76).

    python data/iberia/consistency/screening/rate_outliers.py FLAG SET_ASIDE OUT_DIR

FLAG and SET_ASIDE are result folders (``results/iberia``) of the workflow at the defaults
(``config/config.yaml``, ``config/iberia.yaml``) on ``data/iberia/records.tsv.gz``, through the genus stage:
FLAG with ``screening.rate_outlier_action: flag`` (the default), SET_ASIDE the same genus stage again on the
same ``units/``, ``backbone/`` and ``genera/`` with ``rate_outlier_action: set_aside``. Writes
``rate_outliers.tsv`` and ``rate_outliers_set_aside.tsv`` (the provenance columns of every genus subtree) to
OUT_DIR, and prints the figures quoted in ``data/iberia/consistency/README.md`` and
``docs/technical/configuration.md``.
"""

import csv
import statistics
import sys
from pathlib import Path

csv.field_size_limit(sys.maxsize)
COLUMNS = ["genus", "n_units", "root_to_tip", "rate_outliers", "rate_outlier_z", "crown_age", "crown_source"]


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def subtrees(folder):
    rows = [r for f in sorted((Path(folder) / "subtrees").glob("shard_*.tsv")) for r in read(f)]
    return {r["genus"]: r for r in sorted(rows, key=lambda r: r["genus"]) if r["status"] == "subtree"}


def write(path, rows, columns):
    with open(path, "w", newline="") as handle:
        out = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        out.writeheader()
        out.writerows(rows)


def units(row):
    return row["rate_outliers"].split(";") if row.get("rate_outliers") else []


def main():
    flag_dir, aside_dir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    flag, aside = subtrees(flag_dir), subtrees(aside_dir)
    write(out / "rate_outliers.tsv", flag.values(), COLUMNS)
    write(out / "rate_outliers_set_aside.tsv", aside.values(), COLUMNS[:5] + ["rate_outlier_action"] + COLUMNS[5:])

    assessed = [r for r in flag.values() if r["rate_outlier_z"]]
    outliers = [u for r in flag.values() for u in units(r)]
    genera = [g for g, r in flag.items() if units(r)]
    z = {g: float(r["rate_outlier_z"]) for g, r in flag.items() if r["rate_outlier_z"]}
    print(f"assessed: {len(assessed)} genera (at least four units), of {len(flag)} subtrees")
    print(f"long branches: {len(outliers)} units in {len(genera)} genera; genera above 5: "
          f"{sum(1 for g in genera if z[g] > 5)}, above 10: {sum(1 for g in genera if z[g] > 10)}")
    print("largest scores:", [(g, flag[g]["n_units"], round(z[g], 1)) for g in sorted(genera, key=lambda g: -z[g])[:6]])

    screened = read(flag_dir / "units/screened.tsv")
    translation = {r["unit_id"]: r["translation"] for r in screened}
    numt = [u for u in outliers if translation.get(u) == "numt_like"]
    share_all = sum(1 for r in screened if r["translation"] == "numt_like") / len(screened)
    print(f"numt_like: {len(numt)} of {len(outliers)} long-branch units ({len(numt) / len(outliers):.1%}), "
          f"against {share_all:.1%} of all {len(screened)} units")

    # with set_aside: the same genera, built without their long branches
    moved = [g for g in genera if g in aside]
    rtt = [float(aside[g]["root_to_tip"]) / float(flag[g]["root_to_tip"]) - 1 for g in moved
           if flag[g]["root_to_tip"] and aside[g]["root_to_tip"] and float(flag[g]["root_to_tip"]) > 0]
    set_aside = sum(len(units(aside[g])) for g in moved if aside[g].get("rate_outlier_action") == "set_aside")
    kept_two = [g for g in moved if aside[g].get("rate_outlier_action") not in ("set_aside", "")]
    print(f"set aside: {set_aside} units from {sum(1 for g in moved if aside[g].get('rate_outlier_action') == 'set_aside')} "
          f"subtrees; other actions: {sorted(set(aside[g].get('rate_outlier_action') for g in kept_two))}")
    print(f"mean root-to-tip in those genera: median change {statistics.median(rtt):+.1%}")
    changes = {}
    for g in moved:
        a, b = flag[g]["crown_age"], aside[g]["crown_age"]
        if a and b and abs(float(a) - float(b)) > 1e-4:
            changes[g] = float(b) / float(a) - 1
    rel = list(changes.values())
    print(f"crown ages that change: {len(changes)} of {len(moved)} genera; median {statistics.median(rel):+.1%}; "
          f"older {sum(1 for x in rel if x > 0)}, younger {sum(1 for x in rel if x < 0)}; "
          f"by more than a quarter: {sum(1 for x in rel if abs(x) > 0.25)}")
    big = sorted((g for g in changes if abs(changes[g]) > 0.25), key=lambda g: changes[g])
    print("  more than a quarter:", [(g, round(float(flag[g]["crown_age"]), 1), round(float(aside[g]["crown_age"]), 1)) for g in big])
    same = [g for g in moved if g not in changes]
    from collections import Counter
    print("  unchanged crowns by source:", dict(Counter(flag[g]["crown_source"] for g in same)))
    total_a = sum(float(r["crown_age"]) for r in flag.values() if r["crown_age"])
    total_b = sum(float(r["crown_age"]) for r in aside.values() if r["crown_age"])
    print(f"sum of all crown ages: {total_b / total_a - 1:+.2%}")
    other = [g for g in aside if g not in moved and flag.get(g, {}).get("crown_age") != aside[g]["crown_age"]]
    print("crowns that change in genera without long branches (should be none):", other)


if __name__ == "__main__":
    main()
