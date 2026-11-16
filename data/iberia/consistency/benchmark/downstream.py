"""Run the genus stage, grafting and placement on a backbone outside Snakemake (#81).

    python downstream.py --repo REPO --backbone DIR --units DIR --out DIR
                         [--overlay FILE ...] [--set key.path=value ...] [--genera FILE]
                         [--memo DIR] [--jobs N] [--threads N]

Runs the unchanged workflow scripts ``shard_genera.py``, ``build_shard.py`` (one per shard, ``--jobs`` at a
time), ``graft.py`` and ``place_units.py``, as the rules of ``workflow/rules/genus.smk``, ``graft.smk`` and
``placement.smk`` do, with the settings of ``config/config.yaml`` and ``config/iberia.yaml`` of the repository,
the overlays on top (e.g. ``legacy.yaml``) and then the ``--set`` values.

``--backbone`` is a folder as ``backbone.py`` leaves it (``backbone.tre`` and ``crown_ages.tsv``). ``--units``
holds the units of a run: ``screened_representatives.tsv``, ``screened_unassigned.tsv``, ``filtered.tsv`` and
``representatives.tsv`` (``results/<run>/units`` of a full run, or ``exp/legacy_units``). ``--genera`` (one genus
per line, as ``subsample.py`` writes ``genera.txt``) keeps only those genera: their units, their records, and the
unassigned units whose genera are all among them.

``--memo`` keeps the result of every ``build_genus()`` call in DIR, keyed on everything it reads: the genus, its
ingroup and outgroup sequences, its stem age, the marker and rooting settings, its own entries of the crown
tables, its rate, and the source of ``src/timegraft``. A later call with the same inputs takes the result from
DIR instead of computing it again. ``build_genus()`` is deterministic (``evaluate.py`` checks that the harness
reproduces the full runs), so this changes nothing but the time; it saves most of the genus stage where an
experiment leaves most genera's inputs as they were, as an injected misplacement does. The search of
``place_units.py`` (``vsearch_hits()``, the unplaced units against the records of the units in the tree) is kept
in the same way, keyed on the query and reference sequences and the search settings: it is the same in every run
with the same units in its tree. Both functions are wrapped from the outside; the workflow scripts are not
changed.

OUT gets ``provenance.tsv``, ``grafted.tre``, ``timegraft.tre``, ``units_placed.tsv``, ``subtrees/``, ``genera/``,
``units/representatives.tsv`` (the units of the run, for ``evaluate.py``), ``logs/``, ``steps.tsv`` (seconds per
step and the cache hits) and, last, ``seconds.txt``. The run is made in ``OUT.partial`` and moved to OUT when it
has finished, so OUT exists only for a complete run.
"""

import argparse
import csv
import hashlib
import json
import os
import pickle
import runpy
import shutil
import sys
import time
import types
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from backbone import NS, settings_for

csv.field_size_limit(sys.maxsize)


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows, columns):
    with open(path, "w", newline="") as handle:
        w = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def subset_units(units, genera, out):
    """The units of ``units`` restricted to ``genera`` (None: all), written to ``out``; returns the paths."""
    out.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name in ("screened_representatives.tsv", "screened_unassigned.tsv", "filtered.tsv", "representatives.tsv"):
        source = units / name
        if genera is None:
            os.symlink(source.resolve(), out / name)
        else:
            with open(source, newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                rows = list(reader)
                columns = reader.fieldnames
            if name == "screened_unassigned.tsv":
                # "genera" lists the genera of the unit's records as genus:count; keep a unit only if all are kept
                def keep(r):
                    names = {g.split(":")[0] for g in (r.get("genera") or "").split(";") if g}
                    return bool(names) and names <= genera
            else:
                def keep(r):
                    return r.get("genus") in genera
            write(out / name, [r for r in rows if keep(r)], columns)
        paths[name] = str(out / name)
    return paths


def code_version(repo):
    digest = hashlib.sha256()
    for path in sorted(Path(repo, "src/timegraft").glob("*.py")):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def install_memo(repo, memo):
    """Wrap timegraft.build.build_genus so that its results are kept in ``memo`` (see the module docstring)."""
    import timegraft.build as build

    original = build.build_genus
    version = code_version(repo)
    Path(memo).mkdir(parents=True, exist_ok=True)
    stats = {"hits": 0, "misses": 0}

    def memoised(genus, ingroup, outgroup, stem_age, marker_cfg, known_crowns, timetree_crowns=None, rooting=None,
                 rate=None, crown_basis=None):
        key = json.dumps([version, genus, ingroup, outgroup, repr(float(stem_age)), marker_cfg,
                          repr((known_crowns or {}).get(genus)), repr((timetree_crowns or {}).get(genus)), rooting,
                          repr(rate), (crown_basis or {}).get(genus)], sort_keys=True, default=repr)
        path = Path(memo) / f"{hashlib.sha256(key.encode()).hexdigest()}.pickle"
        if path.exists():
            stats["hits"] += 1
            with open(path, "rb") as handle:
                return pickle.load(handle)
        stats["misses"] += 1
        result = original(genus, ingroup, outgroup, stem_age, marker_cfg, known_crowns, timetree_crowns, rooting,
                          rate, crown_basis)
        partial = path.with_suffix(f".{os.getpid()}.tmp")
        with open(partial, "wb") as handle:
            pickle.dump(result, handle)
        os.replace(partial, path)
        return result

    build.build_genus = memoised
    return stats


def install_search_memo(repo, memo):
    """Wrap timegraft.placement.vsearch_hits in the same way: its result depends only on the query and reference
    sequences and the search settings, which are the same in every run that has the same units in its tree (the
    clean run, the injected errors and the perturbed ages of one setting)."""
    import timegraft.placement as placement

    original = placement.vsearch_hits
    version = code_version(repo)
    Path(memo).mkdir(parents=True, exist_ok=True)
    stats = {"hits": 0, "misses": 0}

    def memoised(queries, refs, min_identity, query_cov=0.9, maxaccepts=20, threads=1):
        digest = hashlib.sha256()
        for part in (version, json.dumps([list(map(list, queries)), repr(min_identity), repr(query_cov),
                                          repr(maxaccepts)])):
            digest.update(part.encode())
        for name, seq in refs:
            digest.update(f"{name}\t{seq}\n".encode())
        path = Path(memo) / f"search_{digest.hexdigest()}.pickle"
        if path.exists():
            stats["hits"] += 1
            with open(path, "rb") as handle:
                return pickle.load(handle)
        stats["misses"] += 1
        result = original(queries, refs, min_identity, query_cov, maxaccepts, threads)
        partial = path.with_suffix(f".{os.getpid()}.tmp")
        with open(partial, "wb") as handle:
            pickle.dump(dict(result), handle)
        os.replace(partial, path)
        return result

    placement.vsearch_hits = memoised
    return stats


def run_shard(job):
    """One build_shard.py, in a process of its own (so its logging goes to its own log file)."""
    repo, sm, memo = job
    sys.path.insert(0, f"{repo}/src")
    stats = install_memo(repo, memo) if memo else {"hits": 0, "misses": 0}
    start = time.time()
    runpy.run_path(f"{repo}/workflow/scripts/build_shard.py", init_globals={"snakemake": sm}, run_name="__main__")
    return time.time() - start, stats


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", required=True)
    ap.add_argument("--backbone", required=True)
    ap.add_argument("--units", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--overlay", action="append", default=[])
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--genera")
    ap.add_argument("--memo")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()
    repo = str(Path(args.repo).resolve())
    final = Path(args.out)
    if (final / "seconds.txt").exists():
        return 0
    out = final.with_name(final.name + ".partial")
    shutil.rmtree(out, ignore_errors=True)
    for sub in ("genera", "subtrees", "logs", "units"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    if args.memo:
        Path(args.memo).mkdir(parents=True, exist_ok=True)
    cfg = settings_for(repo, args.overlay, args.set)
    marker = cfg["markers"][cfg["marker"]]
    if marker.get("dating", {}).get("method", "mpl") != "mpl":
        sys.exit("downstream.py runs the default dating (mpl) only; lsd2 and chronos need their tools")
    rooting = dict(cfg.get("rooting") or {})
    shards = list(range(int(cfg["n_shards"])))
    genera = set(Path(args.genera).read_text().split()) if args.genera else None
    units = subset_units(Path(args.units), genera, out / "units")
    backbone = Path(args.backbone)
    src = f"{repo}/src"
    steps = []

    def script(name, sm):
        start = time.time()
        runpy.run_path(f"{repo}/workflow/scripts/{name}.py", init_globals={"snakemake": sm}, run_name="__main__")
        steps.append({"step": name, "seconds": f"{time.time() - start:.1f}"})

    sys.path.insert(0, src)
    start = time.time()
    script("shard_genera", types.SimpleNamespace(
        input=NS(reps=units["screened_representatives.tsv"], records=units["filtered.tsv"],
                 backbone=str(backbone / "backbone.tre")),
        output=NS(shards=[str(out / f"genera/shard_{i}.tsv") for i in shards], singles=str(out / "genera/singles.tsv"),
                  missing=str(out / "genera/missing_from_backbone.tsv"), taxonomy=str(out / "genera/taxonomy.tsv")),
        params=NS(src=src, n_shards=len(shards), rooting=rooting), log=[str(out / "logs/shard_genera.log")]))

    crown_table = f"{repo}/{cfg['crown_ages']}" if cfg.get("crown_ages") else []
    rate_table = marker["crown_age"].get("rate_table") or []
    jobs = [(repo, types.SimpleNamespace(
        input=NS(shard=str(out / f"genera/shard_{i}.tsv"), crown_ages=crown_table,
                 timetree_crowns=str(backbone / "crown_ages.tsv"), taxonomy=str(out / "genera/taxonomy.tsv"),
                 rate_table=f"{repo}/{rate_table}" if rate_table else [], lsd2=[]),
        output=NS(main=str(out / f"subtrees/shard_{i}.tsv"), set_aside=str(out / f"subtrees/set_aside_{i}.tsv")),
        params=NS(src=src, marker_cfg=marker, rooting=rooting), log=[str(out / f"logs/build_shard_{i}.log")]),
        args.memo) for i in shards]
    begin = time.time()
    hits = misses = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for seconds, stats in pool.map(run_shard, jobs):
            hits, misses = hits + stats["hits"], misses + stats["misses"]
    steps.append({"step": "build_shard", "seconds": f"{time.time() - begin:.1f}",
                  "memo_hits": hits, "memo_misses": misses})

    subtrees = [str(out / f"subtrees/shard_{i}.tsv") for i in shards]
    script("graft", types.SimpleNamespace(
        input=NS(backbone=str(backbone / "backbone.tre"), subtrees=subtrees, singles=str(out / "genera/singles.tsv")),
        output=NS(tree=str(out / "grafted.tre"), provenance=str(out / "provenance.tsv")),
        params=NS(src=src, prune_empty=cfg["prune_empty"]), log=[str(out / "logs/graft.log")]))
    search = install_search_memo(repo, args.memo) if args.memo else {"hits": 0, "misses": 0}
    script("place_units", types.SimpleNamespace(
        input=NS(tree=str(out / "grafted.tre"), unassigned=units["screened_unassigned.tsv"],
                 records=units["filtered.tsv"], reps=units["screened_representatives.tsv"],
                 set_aside=[str(out / f"subtrees/set_aside_{i}.tsv") for i in shards]),
        output=NS(tree=str(out / "timegraft.tre"), placements=str(out / "units_placed.tsv")),
        params=NS(src=src, enabled=cfg.get("place_unassigned_units", True), settings=marker["placement"],
                  rate=marker["crown_age"].get("rate")),
        threads=args.threads, log=[str(out / "logs/place_units.log")]))
    steps[-1].update(memo_hits=search["hits"], memo_misses=search["misses"])
    write(out / "steps.tsv", steps, ["step", "seconds", "memo_hits", "memo_misses"])
    (out / "seconds.txt").write_text(f"{time.time() - start:.1f}\n")
    shutil.rmtree(final, ignore_errors=True)
    os.replace(out, final)
    return 0


if __name__ == "__main__":
    sys.exit(main())
