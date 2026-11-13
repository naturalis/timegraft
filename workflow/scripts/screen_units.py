"""Set aside reference units that do not fit their genus (issues #38, #55)."""

import logging
import sys
import time
from collections import Counter

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

from timegraft.filters import USES_GENETIC_CODE
from timegraft.geneticcode import GeneticCodes
from timegraft.screening import (
    CALIBRATION_COLUMNS, COLUMNS, apply_translation, assess_translation, calibrate, in_scope, not_coding,
    read_threshold_table, screen, search_units, set_aside,
)
from timegraft.seqio import read_tsv, write_tsv

RANKS = ("genus", "family", "order", "class", "phylum")

reps = read_tsv(snakemake.input.reps)
unassigned = read_tsv(snakemake.input.unassigned)
records = {r["record_id"]: r for r in read_tsv(snakemake.input.records)}
cfg = dict(snakemake.params.settings or {})
marker = snakemake.params.marker
reps_columns = list(reps[0]) if reps else ["unit_id", "genus", "n_records", "genus_share", "record_id", "sequence"]
unassigned_columns = list(unassigned[0]) if unassigned else ["unit_id", "reason", "genera", "n_records", "record_id", "sequence"]

# taxonomy of each unit, from the record that represents it
taxon = {r["unit_id"]: records.get(r["record_id"], {}) for r in reps}
class_rank = cfg.get("classes", "class") or "class"
class_of = {u: (t.get(class_rank) or "") for u, t in taxon.items()}
thresholds = dict(read_threshold_table(cfg["threshold_table"])) if cfg.get("threshold_table") else {}
thresholds.update(cfg.get("thresholds") or {})
translation = {**(cfg.get("translation") or {})}

if cfg.get("enabled", False) and reps:
    genus_of = {r["unit_id"]: r["genus"] for r in reps}
    sequences = {r["unit_id"]: r["sequence"] for r in reps}
    # a wider reference: sequences from outside the run, within the scope of the run's taxa
    reference = cfg.get("reference") or {}
    external, scope = {}, "run"
    if reference.get("path"):
        rank = reference.get("scope", "family")
        if rank != "all" and rank not in RANKS:
            raise ValueError(f"screening.reference.scope must be one of {RANKS} or all, not {rank!r}")
        run_taxa = {r: {t.get(r) for t in taxon.values() if t.get(r)} for r in RANKS}
        rows_ref = read_tsv(reference["path"])
        kept_ref = in_scope(rows_ref, run_taxa, rank, genus_of)
        for i, row in enumerate(kept_ref):
            ident = f"external:{row.get('unit_id') or row.get('record_id') or i}"
            external[ident] = row["genus"]
            sequences[ident] = row["sequence"]
        scope = f"run+external:{rank}"
        logging.info("external reference %s: %d of %d sequences within scope %s", reference["path"],
                     len(kept_ref), len(rows_ref), rank)
    coding = marker in USES_GENETIC_CODE and translation.get("enabled", True)
    started = time.time()
    found = search_units([(r["unit_id"], r["sequence"]) for r in reps], min_identity=cfg.get("min_identity", 0.80),
                         query_cov=cfg.get("query_cov", 0.7), maxaccepts=cfg.get("maxaccepts", 10),
                         threads=snakemake.threads, external=[(k, sequences[k]) for k in external],
                         genus_of={**external, **genus_of} if coding else None)
    hits, alignments = found if coding else (found, {})
    logging.info("search: %.1f s on %d threads", time.time() - started, snakemake.threads)
    floor = min([cfg.get("min_identity", 0.80)] + [t.get("min_identity", 1.0) for t in thresholds.values()])
    if floor < cfg.get("min_identity", 0.80):
        logging.warning("per-class min_identity below the search cut-off (%s): hits below it are never seen",
                        cfg.get("min_identity", 0.80))
    rows = screen(genus_of, hits, margin=cfg.get("margin", 0.02), min_identity=cfg.get("min_identity", 0.80),
                  external=external, class_of=class_of, thresholds=thresholds, reference_scope=scope,
                  mutual_pairs=cfg.get("mutual_pairs", "both"))
    if coding:
        settings = snakemake.params.filter or {}
        codes = GeneticCodes(settings.get("genetic_code", 5), settings.get("genetic_codes") or {},
                             settings.get("genetic_code_table") or None)
        code_of = {u: codes.resolve(t)[0] for u, t in taxon.items()}
        assess_translation(rows, alignments, sequences, code_of, translation)
        apply_translation(rows, translation.get("action", "flag"))
    else:
        not_coding(rows)
else:
    rows = [{"unit_id": r["unit_id"], "genus": r["genus"], "decision": "not_screened",
             "reason": "screening disabled for this marker", "threshold_class": class_of.get(r["unit_id"], "")}
            for r in reps]

n_records = {r["unit_id"]: r["n_records"] for r in reps}
for row in rows:
    row["n_records"] = n_records[row["unit_id"]]
aside = set_aside(rows)
decision = {r["unit_id"]: r["decision"] for r in rows}
kept = [r for r in reps if r["unit_id"] not in aside]
# set-aside units join the units without a usable genus, in the same format, to be placed
moved = [{"unit_id": r["unit_id"], "reason": f"screen_{decision[r['unit_id']]}", "genera": f"{r['genus']}:{r['n_records']}",
          "n_records": r["n_records"], "record_id": r["record_id"], "sequence": r["sequence"]}
         for r in reps if r["unit_id"] in aside]
write_tsv(snakemake.output.screened, rows, COLUMNS)
write_tsv(snakemake.output.reps, kept, reps_columns)
write_tsv(snakemake.output.unassigned, unassigned + moved, unassigned_columns)
# observed identities per class, to set per-class thresholds from (#55)
write_tsv(snakemake.output.calibration, calibrate(rows, class_of, floor=cfg.get("min_identity", 0.80))
          if cfg.get("enabled", False) else [], CALIBRATION_COLUMNS)
counts = Counter(decision.values())
logging.info("%d units screened: %s", len(rows), ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
logging.info("%d set aside, in %d genera", len(aside), len({r["genus"] for r in reps if r["unit_id"] in aside}))
flags = Counter(r.get("translation", "") for r in rows)
logging.info("codon comparison: %s", ", ".join(f"{k} {v}" for k, v in sorted(flags.items()) if k))
