# Outgroup benchmark (#19)

Scripts and per-genus results of the comparison of outgroup choices. `outgroup_benchmark.py` documents the methods; `bench_summary.py` and `bench_paired.py` produce the tables. The simulation is `simulate_run300.py` (300 genera, roots known). Run from the repository root with MAFFT on the path, on the output folders unpacked from `data/iberia/` (`--run results/iberia`), after running `timegraft.rooting_diagnostics` on them.

```

## simulated (300 genera)

| method | usable | passed test | same root as reference | crown off >25% | PD off >25% | median dist. to genus |
| --- | --- | --- | --- | --- | --- | --- |
| ref | 300 | 293/300 | 149/292 (51%) | 0/291 | 0/291 | 0.255 |
| sister_random | 287 | 282/287 | 145/281 (52%) | 4/274 | 1/274 | 0.244 |
| far_data | 263 | 255/263 | 149/254 (59%) | 10/247 | 1/247 | 0.280 |
| family_random | 300 | 283/300 | 126/282 (45%) | 60/274 | 27/274 | 0.600 |
| family_near | 300 | 286/300 | 107/285 (38%) | 78/277 | 30/277 | 0.462 |
| tolerant10 | 300 | 300/300 | 160/299 (54%) | 18/291 | 8/291 | 0.417 |
| mad | 300 | 300/300 | 143/299 (48%) | 50/291 | 29/291 | |
| midpoint | 300 | 300/300 | 137/299 (46%) | 50/291 | 29/291 | |

paired against ref: genera where both give a root; b = only method right, c = only ref right; exact two-sided sign test
sister_random  n= 275  method 52% vs ref 52%  b=20 c=20  p=1.000
far_data       n= 248  method 60% vs ref 59%  b=21 c=19  p=0.875
family_random  n= 275  method 45% vs ref 54%  b=23 c=47  p=0.006
family_near    n= 278  method 38% vs ref 53%  b=17 c=57  p=0.000
tolerant10     n= 292  method 54% vs ref 51%  b=24 c=14  p=0.143
mad            n= 292  method 49% vs ref 51%  b=27 c=34  p=0.443
midpoint       n= 292  method 47% vs ref 51%  b=29 c=42  p=0.154


## failed (184 genera)

| method | usable | passed test | same root as reference | crown off >25% | PD off >25% | median dist. to genus |
| --- | --- | --- | --- | --- | --- | --- |
| ref | 184 | 2/184 | 2/2 (100%) | 0/2 | 0/2 | 0.191 |
| sister_random | 131 | 57/131 | 0/1 (0%) | 0/1 | 0/1 | 0.202 |
| far_data | 113 | 101/113 | n/a | 0/0 | 0/0 | 0.274 |
| family_random | 163 | 93/163 | n/a | 0/0 | 0/0 | 0.240 |
| family_near | 163 | 101/163 | n/a | 0/0 | 0/0 | 0.184 |
| tolerant10 | 184 | 184/184 | 1/2 (50%) | 0/2 | 0/2 | 0.217 |
| mad | 184 | 184/184 | 0/2 (0%) | 1/2 | 0/2 | |
| midpoint | 184 | 184/184 | 0/2 (0%) | 1/2 | 0/2 | |

## stable (408 genera)

| method | usable | passed test | same root as reference | crown off >25% | PD off >25% | median dist. to genus |
| --- | --- | --- | --- | --- | --- | --- |
| ref | 408 | 402/408 | 402/402 (100%) | 0/402 | 0/402 | 0.214 |
| sister_random | 350 | 347/350 | 252/342 (74%) | 0/342 | 0/342 | 0.220 |
| far_data | 363 | 338/363 | 206/336 (61%) | 5/336 | 4/336 | 0.220 |
| family_random | 331 | 294/331 | 185/292 (63%) | 8/292 | 1/292 | 0.236 |
| family_near | 331 | 286/331 | 181/284 (64%) | 6/284 | 4/284 | 0.191 |
| tolerant10 | 408 | 408/408 | 278/402 (69%) | 9/402 | 6/402 | 0.228 |
| mad | 408 | 408/408 | 219/402 (54%) | 18/402 | 12/402 | |
| midpoint | 408 | 408/408 | 205/402 (51%) | 23/402 | 13/402 | |

paired against sister_random: genera where both give a root; b = only method right, c = only sister_random right; exact two-sided sign test
sister_random  n= 342  method 74% vs sister_random 74%  b=0 c=0  p=1.000
far_data       n= 286  method 61% vs sister_random 75%  b=23 c=63  p=0.000
family_random  n= 248  method 64% vs sister_random 75%  b=16 c=42  p=0.001
family_near    n= 244  method 64% vs sister_random 77%  b=18 c=50  p=0.000
tolerant10     n= 342  method 72% vs sister_random 74%  b=36 c=42  p=0.572
mad            n= 342  method 55% vs sister_random 74%  b=21 c=85  p=0.000
midpoint       n= 342  method 51% vs sister_random 74%  b=16 c=93  p=0.000
```
