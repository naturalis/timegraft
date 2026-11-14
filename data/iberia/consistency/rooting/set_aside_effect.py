"""What setting aside the misplaced units of a genus would do to its crown age (#75).

    PYTHONPATH=src python data/iberia/consistency/rooting/set_aside_effect.py RUN OUT.tsv [GENUS ...]

RUN is a result folder (``results/iberia``) of a run through the genus stage at the default settings of
``config/config.yaml`` and ``config/iberia.yaml``. For every genus with units among its outgroup
(``misplaced`` in ``subtrees/shard_*.tsv``), or only the genera named, the subtree is built again as
``workflow/scripts/build_shard.py`` builds it: once with all its units, which must reproduce the recorded
crown age, and once without the misplaced units, as an option ``set_aside`` for them would build it (like
``rate_outlier_action: set_aside``). Genera that would be left with fewer than two units are reported, not
built. Needs MAFFT on the path.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import yaml

from timegraft.build import build_genus
from timegraft.rates import Rates

csv.field_size_limit(sys.maxsize)
REPO = Path(__file__).resolve().parents[4]


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def merge(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def main():
    run, target, only = Path(sys.argv[1]), Path(sys.argv[2]), set(sys.argv[3:])
    config = merge(yaml.safe_load(open(REPO / "config/config.yaml")), yaml.safe_load(open(REPO / "config/iberia.yaml")))
    marker_cfg = config["markers"][config["marker"]]
    crown_cfg = marker_cfg["crown_age"]
    rates = Rates(crown_cfg.get("rate"), crown_cfg.get("rates"), crown_cfg.get("rate_table") or None)
    timetree, basis = {}, {}
    for r in read(run / "backbone/crown_ages.tsv"):
        if r.get("crown_age"):
            timetree[r["genus"]] = float(r["crown_age"])
            basis[r["genus"]] = r.get("crown_basis") or "global"
    known = {}
    if config.get("crown_ages"):
        from timegraft.calibrate import read_crown_table
        known = read_crown_table(REPO / config["crown_ages"])
    taxonomy = {r["genus"]: r for r in read(run / "genera/taxonomy.tsv")}
    groups = defaultdict(lambda: {"in": [], "out": [], "stem": None})
    for f in sorted((run / "genera").glob("shard_*.tsv")):
        for r in read(f):
            groups[r["genus"]][r["role"]].append((r["unit_id"], r["sequence"]))
            groups[r["genus"]]["stem"] = float(r["stem_age"])
    prov = {r["genus"]: r for f in sorted((run / "subtrees").glob("shard_*.tsv")) for r in read(f)}

    rows = []
    for genus in sorted(g for g, r in prov.items() if r.get("misplaced") and (not only or g in only)):
        group, recorded = groups[genus], prov[genus]
        misplaced = set(recorded["misplaced"].split(";"))

        def build(ingroup):
            return build_genus(genus, ingroup, group["out"], group["stem"], marker_cfg, known, timetree,
                               rooting=config.get("rooting"), rate=rates.resolve(taxonomy.get(genus)),
                               crown_basis=basis)[1]

        full = build(group["in"])
        if abs(float(full["crown_age"]) - float(recorded["crown_age"])) > 1e-3:
            raise AssertionError(f"{genus}: rebuilt crown {full['crown_age']}, recorded {recorded['crown_age']}")
        kept = [(n, s) for n, s in group["in"] if n not in misplaced]
        row = {"genus": genus, "n_units": len(group["in"]), "n_misplaced": len(misplaced),
               "crown_age": recorded["crown_age"], "crown_source": recorded["crown_source"],
               "crown_age_min": recorded["crown_age_min"], "crown_age_max": recorded["crown_age_max"]}
        if len(kept) < 2:
            row.update(set_aside="fewer than two units left")
        else:
            alt = build(kept)
            row.update(set_aside="built", crown_age_without=alt["crown_age"], crown_source_without=alt["crown_source"],
                       violations_without=alt["violations"], misplaced_without=alt["misplaced"],
                       crown_factor=f"{max(float(alt['crown_age']), float(full['crown_age'])) / min(float(alt['crown_age']), float(full['crown_age'])):.3f}",
                       crown_age_min_without=alt["crown_age_min"], crown_age_max_without=alt["crown_age_max"])
        rows.append(row)
        print(genus, row.get("crown_age"), "->", row.get("crown_age_without", row["set_aside"]), flush=True)
    columns = list(dict.fromkeys(k for r in rows for k in r))
    with open(target, "w", newline="") as handle:
        out = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n", restval="")
        out.writeheader()
        out.writerows(rows)


if __name__ == "__main__":
    main()
