"""Screening diagnostics on the Iberian units (#55).

Runs the screening of issue #55 on the representatives of the archive and
measures, against known errors put into the data:

- ``injection.tsv``: recall of mislabelled units (a unit given the name of
  another genus of its family, or of another order) for a grid of ``margin``
  and ``min_identity``, with the collateral flags that the error causes, and
  the share of units set aside without any injected error, with both units
  of a mutual pair set aside or only the further one (``mutual_pairs``);
- ``injection_classes.tsv``: the same at the defaults and at the thresholds
  suggested per class by the calibration, per class;
- ``numt_simulation.tsv``: recall of copies with neutral (unconstrained)
  substitutions, as a nuclear copy accumulates them, for a grid of the codon
  thresholds, against the flag rate on real pairs of nearest relatives;
- ``external.tsv``: the wider reference, by splitting the genera of each
  family between the run and an external reference, against the decisions
  with all genera in the run;
- ``runtime.tsv``: time per step;
- ``screened.tsv.gz`` and ``calibration.tsv``: the screening at the defaults,
  as ``screen_units`` writes them.

Usage (VSEARCH on the PATH):

    python benchmark.py REPRESENTATIVES FILTERED_RECORDS OUTDIR [THREADS]

with ``units/representatives.tsv`` of ``iberia-run-tolerant.tar.gz`` and the
records of ``records.tsv.gz`` without those in its ``units/rejected.tsv``.
"""

import csv
import gzip
import random
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

from timegraft.geneticcode import GeneticCodes, has_open_frame
from timegraft.screening import (
    CALIBRATION_COLUMNS, COLUMNS, SET_ASIDE, assess_translation, calibrate, compare_codons, judge_codons, screen,
    search_units,
)

reps_path, records_path, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
threads = int(sys.argv[4]) if len(sys.argv) > 4 else 2
out.mkdir(parents=True, exist_ok=True)


def read(path):
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(name, rows, columns=None):
    columns = columns or list(rows[0])
    opener = gzip.open if name.endswith(".gz") else open
    with opener(out / name, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


runtime = []


def timed(step, f, *args, **kwargs):
    start = time.time()
    value = f(*args, **kwargs)
    runtime.append({"step": step, "seconds": f"{time.time() - start:.1f}", "threads": threads})
    return value


reps = read(reps_path)
records = {r["record_id"]: r for r in read(records_path)}
taxon = {r["unit_id"]: records[r["record_id"]] for r in reps}
genus_of = {r["unit_id"]: r["genus"] for r in reps}
seqs = {r["unit_id"]: r["sequence"] for r in reps}
class_of = {u: t.get("class", "") for u, t in taxon.items()}
codes = GeneticCodes(5)
code_of = {u: codes.resolve(t)[0] for u, t in taxon.items()}
pairs = [(r["unit_id"], r["sequence"]) for r in reps]

# 1. The screening at the defaults --------------------------------------------
hits = timed("search, all against all", search_units, pairs, 0.80, 0.7, 10, threads)
hits, alignments = timed("search, all against all, with alignments and a second pass for one-sided hits",
                         search_units, pairs, 0.80, 0.7, 10, threads,
                         genus_of=genus_of)
rows = timed("decisions", screen, genus_of, hits)
timed("codon comparison", assess_translation, rows, alignments, seqs, code_of)
n_records = {r["unit_id"]: r["n_records"] for r in reps}
for row in rows:
    row["n_records"] = n_records[row["unit_id"]]
calibration = timed("calibration", calibrate, rows, class_of)
write("screened.tsv.gz", rows, COLUMNS)
write("calibration.tsv", calibration, CALIBRATION_COLUMNS)
suggested = {c["class"]: {"min_identity": float(c["suggested_min_identity"])}
             for c in calibration if c["suggested_min_identity"]}
base_decision = {r["unit_id"]: r["decision"] for r in rows}
print(Counter(base_decision.values()), Counter(r["translation"] for r in rows))

# 2. Mislabelled units ----------------------------------------------------------
members = defaultdict(list)
for unit, genus in genus_of.items():
    members[genus].append(unit)
genus_taxon = {genus_of[u]: taxon[u] for u in genus_of}
big = {g for g, us in members.items() if len(us) >= 3}
by_family, by_phylum = defaultdict(set), defaultdict(set)
for g in big:
    by_family[genus_taxon[g].get("family", "")].add(g)
    by_phylum[genus_taxon[g].get("phylum", "")].add(g)
eligible = [u for u in genus_of if genus_of[u] in big and base_decision[u] == "keep"]


def inject(seed, kind, share=0.02):
    """Give ``share`` of the eligible units the name of another genus: in the family, or in another order."""
    rng = random.Random(seed)
    labels, injected = dict(genus_of), set()
    for unit in rng.sample(eligible, int(share * len(eligible))):
        t = genus_taxon[genus_of[unit]]
        if kind == "family":
            options = sorted(by_family[t.get("family", "")] - {genus_of[unit]})
        else:  # another order of the same phylum
            options = sorted(g for g in by_phylum[t.get("phylum", "")] if genus_taxon[g].get("order") != t.get("order"))
        if options:
            labels[unit] = rng.choice(options)
            injected.add(unit)
    return labels, injected


def flagged(labels, **kwargs):
    return {r["unit_id"] for r in screen(labels, hits, **kwargs) if r["decision"] in SET_ASIDE}


grid = [(m, t, p) for m in (0.0, 0.01, 0.02, 0.03, 0.05) for t in (0.80, 0.82, 0.85, 0.88) for p in ("both", "further")]
baseline = {(m, t, p): flagged(genus_of, margin=m, min_identity=t, mutual_pairs=p) for m, t, p in grid}
start = time.time()
results, per_class = [], defaultdict(Counter)
for kind in ("family", "order"):
    for seed in range(5):
        labels, injected = inject(seed, kind)
        for m, t, p in grid:
            found = flagged(labels, margin=m, min_identity=t, mutual_pairs=p)
            tp = len(found & injected)
            collateral = len(found - injected - baseline[(m, t, p)])
            results.append({"kind": kind, "seed": seed, "margin": m, "min_identity": t, "mutual_pairs": p,
                            "injected": len(injected),
                            "recall": f"{tp / len(injected):.4f}",
                            "precision_of_changes": f"{tp / (tp + collateral):.4f}" if tp + collateral else "",
                            "collateral": collateral,
                            "baseline_set_aside_share": f"{len(baseline[(m, t, p)]) / len(genus_of):.4f}"})
        for setting, kwargs in (("default", {}), ("calibrated", {"class_of": class_of, "thresholds": suggested})):
            found = flagged(labels, **kwargs)
            base = flagged(genus_of, **kwargs)
            for unit in injected:
                per_class[(kind, setting, class_of[unit])]["injected"] += 1
                per_class[(kind, setting, class_of[unit])]["found"] += unit in found
            for unit in base:
                per_class[(kind, setting, class_of[unit])]["baseline_set_aside"] += 1
runtime.append({"step": f"injection grid ({len(results)} screenings of injected labels)", "seconds": f"{time.time() - start:.1f}",
                "threads": 1})

summary = defaultdict(list)
for r in results:
    summary[(r["kind"], r["margin"], r["min_identity"], r["mutual_pairs"])].append(r)
write("injection.tsv", [
    {"kind": k, "margin": m, "min_identity": t, "mutual_pairs": p, "replicates": len(v), "injected": v[0]["injected"],
     "recall": f"{sum(float(x['recall']) for x in v) / len(v):.4f}",
     "precision_of_changes": f"{sum(float(x['precision_of_changes'] or 0) for x in v) / len(v):.4f}",
     "collateral": f"{sum(x['collateral'] for x in v) / len(v):.1f}",
     "baseline_set_aside_share": v[0]["baseline_set_aside_share"]}
    for (k, m, t, p), v in sorted(summary.items())])
units_per_class = Counter(class_of.values())
write("injection_classes.tsv", [
    {"kind": k, "setting": s, "class": c, "units": units_per_class[c], "injected": n["injected"],
     "recall": f"{n['found'] / n['injected']:.4f}" if n["injected"] else "",
     "baseline_set_aside_share": f"{n['baseline_set_aside'] / 5 / units_per_class[c]:.4f}",
     "min_identity": suggested.get(c, {}).get("min_identity", 0.80) if s == "calibrated" else 0.80}
    for (k, s, c), n in sorted(per_class.items()) if units_per_class[c] >= 200])

# 3. Nuclear copies -------------------------------------------------------------
rng = random.Random(11)
candidates = [u for u in genus_of if has_open_frame(seqs[u], code_of[u]) and set(seqs[u].upper()) <= set("ACGT-")]
sample = rng.sample(candidates, 2000)
stats_by_rate = {}
for rate in (0.01, 0.02, 0.04, 0.08):
    stats = []
    for unit in sample:
        source = seqs[unit].replace("-", "")
        copy = list(source)
        for i in rng.sample(range(len(copy)), max(1, round(rate * len(copy)))):
            copy[i] = rng.choice([b for b in "ACGT" if b != copy[i]])
        copy = "".join(copy)
        s = compare_codons(1, copy, source, copy, code_of[unit])
        s["passes_filter"] = has_open_frame(copy, code_of[unit])
        stats.append(s)
    stats_by_rate[rate] = stats
real = [{"codons_changed": int(r["codons_changed"]), "aa_changes": int(r["aa_changes"]), "frameshift_stretches": 0}
        for r in rows if r.get("codons_changed")]
numt = []
for min_changed in (4, 8, 12):
    for share in (0.3, 0.4, 0.5, 0.6):
        settings = {"min_changed_codons": min_changed, "max_nonsynonymous_share": share}
        real_rate = sum(judge_codons(s, settings)[0] == "numt_like" for s in real) / len(real)
        for rate, stats in stats_by_rate.items():
            passing = [s for s in stats if s["passes_filter"]]
            numt.append({"min_changed_codons": min_changed, "max_nonsynonymous_share": share, "substitution_rate": rate,
                         "copies": len(stats), "pass_reading_frame_filter": len(passing),
                         "recall_all": f"{sum(judge_codons(s, settings)[0] == 'numt_like' for s in stats) / len(stats):.4f}",
                         "recall_after_filter": f"{sum(judge_codons(s, settings)[0] == 'numt_like' for s in passing) / len(passing):.4f}"
                         if passing else "",
                         "flag_rate_real_pairs": f"{real_rate:.4f}"})
write("numt_simulation.tsv", numt)

# 4. A wider reference ----------------------------------------------------------
rng = random.Random(3)
run_genera = set()
family_genera = defaultdict(list)
for genus in sorted(members):
    family_genera[genus_taxon[genus].get("family", "")].append(genus)
for family, genera in sorted(family_genera.items()):
    rng.shuffle(genera)
    run_genera.update(genera[: (len(genera) + 1) // 2])
run_units = [u for u in genus_of if genus_of[u] in run_genera]
ext_units = [u for u in genus_of if genus_of[u] not in run_genera]
run_families = {taxon[u].get("family") for u in run_units}
ext_in_scope = [u for u in ext_units if taxon[u].get("family") in run_families]
run_pairs = [(u, seqs[u]) for u in run_units]
run_genus = {u: genus_of[u] for u in run_units}
ext_ids = {f"external:{u}": genus_of[u] for u in ext_in_scope}
hits_run = timed("search, run half only", search_units, run_pairs, 0.80, 0.7, 10, threads,
                 genus_of=run_genus, alignments=False)
hits_ext = timed("search, run half against run and external reference", search_units, run_pairs, 0.80, 0.7, 10,
                 threads, external=[(f"external:{u}", seqs[u]) for u in ext_in_scope],
                 genus_of={**run_genus, **ext_ids}, alignments=False)
alone = {r["unit_id"]: r for r in screen(run_genus, hits_run)}
wider = {r["unit_id"]: r for r in screen(run_genus, hits_ext, external=ext_ids, reference_scope="run+external:family")}
full = {u: base_decision[u] for u in run_units}
ext_rows = []
for name, decisions in (("run half only", alone), ("run half with external reference", wider)):
    aside = {u for u, r in decisions.items() if r["decision"] in SET_ASIDE}
    truth = {u for u in run_units if full[u] in SET_ASIDE}
    ext_rows.append({
        "screening": name, "units": len(run_units), "reference_sequences": len(ext_in_scope) if "external" in name else 0,
        "set_aside": len(aside), "set_aside_with_all_genera": len(truth), "agree": len(aside & truth),
        "missed": len(truth - aside), "extra": len(aside - truth),
        "other_genus_by_external_hit": sum(1 for r in decisions.values()
                                           if r["decision"] == "other_genus" and r.get("best_other_scope") == "external"),
    })
write("external.tsv", ext_rows)
write("runtime.tsv", runtime)
for row in ext_rows:
    print(row)
