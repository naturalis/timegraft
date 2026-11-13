"""Keep a random half of the genera of a run (#56).

    python subsample.py UNITS_DIR SEED OUT_DIR [SHARE]

Writes ``representatives.tsv``, ``filtered.tsv`` and ``screened.tsv`` with the rows of the genera kept;
in ``screened.tsv`` a nearest unit of a dropped genus is blanked.
The genera are those of the representatives, sorted, and drawn with the seed.
"""

import csv
import random
import sys
from pathlib import Path

units, seed, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
share = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
csv.field_size_limit(10**9)
out.mkdir(parents=True, exist_ok=True)
genera = sorted({r["genus"] for r in csv.DictReader(open(units / "representatives.tsv"), delimiter="\t")})
keep = set(random.Random(int(seed)).sample(genera, round(share * len(genera))))
for name in ("representatives.tsv", "filtered.tsv", "screened.tsv"):
    with open(units / name, newline="") as src, open(out / name, "w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(dst, reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for r in reader:
            if r["genus"] not in keep:
                continue
            if name == "screened.tsv" and r.get("best_other_genus") and r["best_other_genus"] not in keep:
                # a nearest unit of a genus that was dropped is no evidence in the smaller run
                for field in ("best_other", "best_other_genus", "best_other_identity", "best_other_scope"):
                    if field in r:
                        r[field] = ""
            writer.writerow(r)
(out / "genera.txt").write_text("\n".join(sorted(keep)) + "\n")
print(f"seed {seed}: {len(keep)} of {len(genera)} genera kept")
