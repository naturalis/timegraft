"""Perturb every TimeTree age in a copy of the cache (#56).

    python perturb_cache.py CACHE SEED OUT [SIGMA]

Each cached age (and everything stored with it: adjusted age, interval, study estimates) is multiplied
by exp(SIGMA z), with z standard normal and drawn per query from the seed and the query's key, so that
the factor of a query does not depend on the order or the content of the rest of the cache.
"""

import json
import math
import random
import sys

cache, seed, out = sys.argv[1:4]
sigma = float(sys.argv[4]) if len(sys.argv) > 4 else 0.3
data = json.load(open(cache))
n = 0
for key, entry in data.items():
    if not key.startswith("mrca:") or not isinstance(entry, dict) or entry.get("age") is None:
        continue
    f = math.exp(sigma * random.Random(f"{seed}:{key}").gauss(0.0, 1.0))
    for field in ("age", "adjusted_age", "ci_low", "ci_high", "fresh_age"):
        if entry.get(field) is not None:
            entry[field] = entry[field] * f
    if entry.get("estimates"):
        entry["estimates"] = [e * f for e in entry["estimates"]]
    for study in entry.get("studies") or []:
        if study.get("age") is not None:
            study["age"] = study["age"] * f
    n += 1
json.dump(data, open(out, "w"))
print(f"seed {seed}: {n} ages perturbed, sigma {sigma}")
