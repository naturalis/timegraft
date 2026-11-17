"""Merge API caches of pipeline runs and benchmark builds into one (#87).

    python merge_caches.py OUT CACHE [CACHE ...]

Each CACHE is a cache as the pipeline writes it (``resources/api_cache.json``, ``exp/*/cache.json``), plain or
gzipped (``.gz``); OUT is written gzipped if its name ends in ``.gz``. The caches are merged key by key. The
merge stops with an error:

- if two caches give different answers to the same question (the benchmark would no longer be reproducible
  from the merged cache);
- if a cache is a perturbed copy (``exp/noise/``): its ages are not TimeTree's;
- if a TimeTree answer is not in cache format 2 (#47), which the current pipeline writes: its version marker
  ``v`` must be at least 2, as ``timegraft.timetree.is_legacy()`` reads it.

Used to build ``timetree_cache.json.gz``, the cache every build of ``run_experiments.sh`` starts from, and to add
the answers of a later run to it.
"""

import gzip
import json
import sys
from collections import Counter
from pathlib import Path


CACHE_VERSION = 2  # timegraft.timetree.CACHE_VERSION


def load(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return json.load(handle)


def merge(paths):
    merged, origin = {}, {}
    for path in paths:
        if "/noise/" in str(path):
            sys.exit(f"{path} is a perturbed copy of the cache; leave it out")
        for key, value in load(path).items():
            if key.startswith("mrca:") and not (isinstance(value, dict) and value.get("v", 1) >= CACHE_VERSION):
                sys.exit(f"{path}: {key[:80]} is not in cache format 2 (#47)")
            if key in merged and merged[key] != value:
                sys.exit(f"{path} and {origin[key]} answer {key[:80]} differently")
            merged.setdefault(key, value)
            origin.setdefault(key, path)
    return merged


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out, paths = Path(sys.argv[1]), [Path(p) for p in sys.argv[2:]]
    merged = merge(paths)
    text = json.dumps(merged, sort_keys=True, separators=(",", ":"))
    if out.suffix == ".gz":
        # no file name and mtime 0 in the header, so that the same cache gives the same bytes wherever it is written
        with open(out, "wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", compresslevel=9,
                                                    mtime=0) as handle:
            handle.write(text.encode())
    else:
        out.write_text(text)
    kinds = Counter(k.split(":")[0] for k in merged)
    print(f"{len(paths)} caches, {len(merged)} answers ({', '.join(f'{n} {k}' for k, n in sorted(kinds.items()))}); "
          f"{out}: {out.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
