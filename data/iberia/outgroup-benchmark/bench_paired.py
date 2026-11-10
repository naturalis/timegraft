import csv, sys, math
from collections import defaultdict
rows=list(csv.DictReader(open(sys.argv[1]),delimiter='\t')); base=sys.argv[2] if len(sys.argv)>2 else 'ref'
by=defaultdict(dict)
for r in rows:
    if r.get('same_root') in ('True','False'): by[r['genus']][r['method']]=r['same_root']=='True'
print(f"paired against {base}: genera where both give a root; b = only method right, c = only {base} right; exact two-sided sign test")
for m in ['sister_random','far_data','family_random','family_near','tolerant10','mad','midpoint']:
    pairs=[(d[m],d[base]) for d in by.values() if m in d and base in d]
    if not pairs: continue
    b=sum(x and not y for x,y in pairs); c=sum(y and not x for x,y in pairs); n=b+c
    p=min(1,2*sum(math.comb(n,k) for k in range(0,min(b,c)+1))/2**n) if n else 1
    print(f"{m:14s} n={len(pairs):4d}  method {sum(x for x,_ in pairs)/len(pairs):.0%} vs {base} {sum(y for _,y in pairs)/len(pairs):.0%}  b={b} c={c}  p={p:.3f}")
