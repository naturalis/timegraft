"""Crown ages of every genus under the four combinations of crown basis and plausibility test (#54).

Usage: python crown_compare.py <subtrees dir> <global crown_ages.tsv> <sampled crown_ages.tsv> [out.tsv]
Recomputes the crown age decision from the provenance of a genus-stage run (root-to-tip, rate, stem),
which is all the decision depends on, so the trees need not be rebuilt for each combination.
"""
import csv, glob, math, statistics, sys
import dendropy
from timegraft.dating import crown_age

folder, global_path, sampled_path = sys.argv[1:4]
out = sys.argv[4] if len(sys.argv) > 4 else None
prov = {}
for path in sorted(glob.glob(f"{folder}/subtrees_*.tsv")):
    for r in csv.DictReader(open(path), delimiter="\t"):
        prov[r["genus"]] = r
def crowns(path):
    return {r["genus"]: float(r["crown_age"]) for r in csv.DictReader(open(path), delimiter="\t") if r["crown_age"]}
tables = {"global": crowns(global_path), "sampled": crowns(sampled_path)}
plaus = {"enabled": True, "rate_factor": 10, "saturation": 0.15, "young_factor": 3}
rows = []
for genus, r in sorted(prov.items()):
    d, rate, stem = float(r["root_to_tip"]), float(r["rate"]), float(r["stem_age"])
    tree = dendropy.Tree.get(data=f"(A:{d},B:{d});", schema="newick", rooting="force-rooted")
    row = {"genus": genus, "root_to_tip": d, "stem_age": stem, "barcode_age": d / rate}
    for basis in ("global", "sampled"):
        for test in (False, True):
            c = crown_age(genus, stem, tree, {}, method="rate", rate=rate, suggested=tables[basis],
                          estimate="timetree_midpoint", plausibility=plaus if test else None)
            key = f"{basis}{'+test' if test else ''}"
            row[key] = c.age
            row[key + "_source"] = c.source
    rows.append(row)
# check: the run itself used sampled+test
mismatch = [r["genus"] for r in rows if abs(r["sampled+test"] - float(prov[r["genus"]]["crown_age"])) > 1e-3]
print(f"{len(rows)} genera; recomputed differs from the run for {len(mismatch)}: {mismatch[:5]}")
for key in ("global", "global+test", "sampled", "sampled+test"):
    ratio = [r[key] / r["barcode_age"] for r in rows]
    tt = sum(1 for r in rows if r[key + "_source"] == "timetree")
    print(f"{key:>13}: TimeTree used for {tt:4d}; crown over barcode age > 10: {sum(x > 10 for x in ratio):3d}, "
          f"> 20: {sum(x > 20 for x in ratio):3d}; median {statistics.median(ratio):.2f}; "
          f"low divergence (< 0.03) and > 10: {sum(1 for r, x in zip(rows, ratio) if r['root_to_tip'] < 0.03 and x > 10)}")
changed = [r for r in rows if abs(math.log(r["sampled+test"] / r["global"])) > 0.01]
print(f"crown ages that change (new against before #54): {len(changed)}; older {sum(r['sampled+test'] > r['global'] for r in changed)}")
if out:
    with open(out, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
