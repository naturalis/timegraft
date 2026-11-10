"""Compare genus rooting and ages between the strict and the tolerant Iberian runs."""
import csv, glob, statistics as st, sys
from collections import Counter

import dendropy

old, new = sys.argv[1], sys.argv[2]
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
sub = lambda run: {r["genus"]: r for p in sorted(glob.glob(f"{run}/subtrees/shard_*.tsv")) for r in rd(p)}
A, B = sub(old), sub(new)
print("genus subtrees:", len(A), len(B))
print("rooting, strict run:  ", dict(Counter(r["rooting"] for r in A.values())))
print("rooting, tolerant run:", dict(Counter(r["rooting"] for r in B.values())))
v = [int(r["violations"]) for r in B.values() if r["violations"] not in ("", "0")]
print("units on the wrong side where >0: n", len(v), "median", st.median(v), "max", max(v),
      "| 1:", sum(x == 1 for x in v), "2:", sum(x == 2 for x in v), "3-5:", sum(3 <= x <= 5 for x in v), ">5:", sum(x > 5 for x in v))
mis = [r["misplaced"].split(";") for r in B.values() if r["misplaced"]]
print("genera with genus units among the outgroup:", len(mis), "units:", sum(map(len, mis)))
share = [int(r["violations"]) / (int(r["n_units"]) + len(r["outgroup"].split(";"))) for r in B.values() if r["violations"] not in ("", "0")]
print("units on the wrong side as share of genus + outgroup: median", f"{st.median(share):.0%}", "max", f"{max(share):.0%}")
both = [g for g in A if g in B]
ratio = [float(B[g]["crown_age"]) / float(A[g]["crown_age"]) for g in both]
off = lambda xs: sum(x < 0.8 or x > 1.25 for x in xs)
print(f"crown age tolerant/strict: median {st.median(ratio):.2f}, off by >25%: {off(ratio)}/{len(ratio)}")
for label, pred in [("strict rooted on outgroup", lambda g: A[g]["rooting"] == "outgroup"),
                    ("strict fell back to midpoint", lambda g: A[g]["rooting"] != "outgroup")]:
    r = [float(B[g]["crown_age"]) / float(A[g]["crown_age"]) for g in both if pred(g)]
    print(f"  {label}: n {len(r)}, median {st.median(r):.2f}, off by >25%: {off(r)}")
print("crown sources strict:  ", dict(Counter(r["crown_source"] for r in A.values())))
print("crown sources tolerant:", dict(Counter(r["crown_source"] for r in B.values())))


def pd_total(path):
    t = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True, rooting="force-rooted")
    return sum(e.length or 0 for e in t.postorder_edge_iter() if e.tail_node is not None), len(t.leaf_nodes())


pa, na = pd_total(f"{old}/timegraft.tre")
pb, nb = pd_total(f"{new}/timegraft.tre")
print(f"total PD (Myr): strict {pa:,.0f} ({na} tips), tolerant {pb:,.0f} ({nb} tips), ratio {pb / pa:.3f}")
