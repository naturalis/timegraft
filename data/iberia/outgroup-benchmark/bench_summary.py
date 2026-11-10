import csv, statistics as st, sys
from collections import defaultdict

METHODS = ["ref", "sister_random", "far_data", "family_random", "family_near", "tolerant10", "mad", "midpoint"]
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
for s in sorted({r["set"] for r in rows}):
    sub = [r for r in rows if r["set"] == s]
    n = len({r["genus"] for r in sub})
    print(f"\n## {s} ({n} genera)\n")
    print("| method | usable | passed test | same root as reference | crown off >25% | PD off >25% | median dist. to genus |")
    print("| --- | --- | --- | --- | --- | --- | --- |")
    for m in METHODS:
        rs = [r for r in sub if r["method"] == m]
        usable = [r for r in rs if r["status"] != "no candidates"]
        ok = [r for r in usable if r["status"] == "ok"]
        same = [r for r in ok if r.get("same_root") in ("True", "False")]
        sr = sum(r["same_root"] == "True" for r in same)
        cr = [float(r["crown_ratio"]) for r in ok if r.get("crown_ratio")]
        pr = [float(r["pd_ratio"]) for r in ok if r.get("pd_ratio")]
        off = lambda v: sum(x < 0.8 or x > 1.25 for x in v)
        d = [float(r["dist"]) for r in usable if r.get("dist")]
        print(f"| {m} | {len(usable)} | {len(ok)}/{len(usable)} | "
              f"{sr}/{len(same)} ({100 * sr / len(same):.0f}%) | " if same else f"| {m} | {len(usable)} | {len(ok)}/{len(usable)} | n/a | ", end="")
        print(f"{off(cr)}/{len(cr)} | {off(pr)}/{len(pr)} | {st.median(d):.3f} |" if d else f"{off(cr)}/{len(cr)} | {off(pr)}/{len(pr)} | |")
