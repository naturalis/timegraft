"""Exploratory: nearest relative of every unit representative among all others, before screening existed.

Usage: unit_matches.py HITS INDEX TAXONOMY OUT
  HITS: blast6 output of vsearch --usearch_global reps.fa --db reps.fa --self --id 0.80 --maxaccepts 10 --maxrejects 64 --query_cov 0.7
  INDEX: tab-separated neutral name, unit_id, genus, n_records, length (as written alongside reps.fa)
  TAXONOMY: genera/taxonomy.tsv of the run
"""
import csv, sys, statistics as st
from collections import Counter, defaultdict
hits, index, taxonomy, out = sys.argv[1:5]
unit, genus, nrec = {}, {}, {}
for l in open(index):
    k, u, g, n, _ = l.rstrip("\n").split("\t"); unit[k] = u; genus[k] = g; nrec[k] = int(n or 1)
tax = {r["genus"]: r for r in csv.DictReader(open(taxonomy), delimiter="\t")}
RANKS = ("family", "order", "class", "phylum")
def shared(a, b):
    if a == b: return "genus"
    for r in RANKS:
        x, y = tax.get(a, {}).get(r, ""), tax.get(b, {}).get(r, "")
        if x and x == y: return r
    return "none"
best_same, best_other = {}, {}
for l in open(hits):
    f = l.split("\t"); q, t, ident = f[0], f[1], float(f[2]) / 100
    if genus[q] == genus[t]:
        if ident > best_same.get(q, (0,))[0]: best_same[q] = (ident, t)
    elif ident > best_other.get(q, (0,))[0]: best_other[q] = (ident, t)
size = Counter(genus.values())
rows = []
for k in unit:
    s, o = best_same.get(k, (0, ""))[0], best_other.get(k, (0, ""))
    level = shared(genus[k], genus[o[1]]) if o[1] else ""
    rows.append({"unit_id": unit[k], "genus": genus[k], "genus_units": size[genus[k]], "n_records": nrec[k],
                 "best_same_genus": f"{s:.3f}" if s else "", "best_other": f"{o[0]:.3f}" if o[0] else "",
                 "other_genus": genus[o[1]] if o[1] else "", "other_unit": unit[o[1]] if o[1] else "", "shared_rank": level})
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
multi = [r for r in rows if r["genus_units"] >= 2]
print(f"{len(rows)} units; {len(multi)} in genera with at least two units")
def f(x): return float(x) if x else 0.0
print("\nNearest relative in the data, units of genera with at least two units")
cat = Counter()
for r in multi:
    s, o = f(r["best_same_genus"]), f(r["best_other"])
    if not s and not o: cat["no match at 80% anywhere"] += 1
    elif s >= o: cat["nearest is a congener"] += 1
    elif o - s < 0.02: cat["another genus, by less than 2%"] += 1
    else: cat[f"another genus by 2% or more, shared {r['shared_rank']}"] += 1
for k, v in cat.most_common(): print(f"   {k:<48}{v:>6}  {v/len(multi):.1%}")
print("\nUnits nearest to another genus by 2% or more: identity to that genus, by shared rank")
for rank in ("family", "order", "class", "phylum", "none"):
    sel = [r for r in multi if f(r["best_other"]) - f(r["best_same_genus"]) >= 0.02 and r["shared_rank"] == rank]
    if sel:
        ids = [f(r["best_other"]) for r in sel]
        print(f"   {rank:<8}{len(sel):>5}  identity median {st.median(ids):.3f}; at 97% or more: {sum(i >= .97 for i in ids)}; "
              f"no congener at 80%: {sum(not r['best_same_genus'] for r in sel)}; single-record BINs: {sum(r['n_records'] == 1 for r in sel)}")
print("\nNo congener within 80% identity (isolated within the genus):", sum(not r["best_same_genus"] for r in multi),
      "; within 85%:", sum(f(r["best_same_genus"]) < .85 for r in multi), "; within 90%:", sum(f(r["best_same_genus"]) < .90 for r in multi))
print("single-record BINs overall:", f"{sum(r['n_records']==1 for r in rows)/len(rows):.0%}")
# clear contaminants / gross misidentifications: nearly identical to a unit of another order or beyond
gross = [r for r in rows if f(r["best_other"]) >= 0.95 and r["shared_rank"] in ("class", "phylum", "none") and f(r["best_other"]) - f(r["best_same_genus"]) >= 0.02]
print(f"\n95% or more identical to a unit of another order, class or phylum: {len(gross)}")
for r in sorted(gross, key=lambda r: -f(r["best_other"]))[:12]:
    print("   ", r["unit_id"], r["genus"], f"({tax.get(r['genus'],{}).get('class','')})", "->", r["other_genus"], f"({tax.get(r['other_genus'],{}).get('class','')})", r["best_other"], "same genus", r["best_same_genus"] or "-", "records", r["n_records"])
g = Counter(r["genus"] for r in multi if f(r["best_other"]) - f(r["best_same_genus"]) >= 0.02)
gs = [(n / size[x], n, x) for x, n in g.items()]
print("\ngenera with any such unit:", len(g), "of", sum(1 for x in size.values() if x >= 2), "; with a quarter or more of their units:", sum(a >= .25 for a, _, _ in gs))
print("most affected large genera:", [(x, n, size[x]) for a, n, x in sorted(gs, key=lambda t: -t[1])[:12]])
