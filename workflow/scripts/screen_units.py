"""Set aside reference units that do not fit their genus (issue #38)."""

import logging
import sys
from collections import Counter

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

from timegraft.screening import COLUMNS, screen, search_units, set_aside
from timegraft.seqio import read_tsv, write_tsv

reps = read_tsv(snakemake.input.reps)
unassigned = read_tsv(snakemake.input.unassigned)
cfg = dict(snakemake.params.settings or {})
reps_columns = list(reps[0]) if reps else ["unit_id", "genus", "n_records", "genus_share", "record_id", "sequence"]
unassigned_columns = list(unassigned[0]) if unassigned else ["unit_id", "reason", "genera", "n_records", "record_id", "sequence"]

if cfg.get("enabled", False) and reps:
    hits = search_units([(r["unit_id"], r["sequence"]) for r in reps], min_identity=cfg.get("min_identity", 0.80),
                        query_cov=cfg.get("query_cov", 0.7), maxaccepts=cfg.get("maxaccepts", 10),
                        threads=snakemake.threads)
    rows = screen({r["unit_id"]: r["genus"] for r in reps}, hits,
                  margin=cfg.get("margin", 0.02), min_identity=cfg.get("min_identity", 0.80))
else:
    rows = [{"unit_id": r["unit_id"], "genus": r["genus"], "decision": "not_screened"} for r in reps]

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
counts = Counter(decision.values())
logging.info("%d units screened: %s", len(rows), ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
logging.info("%d set aside, in %d genera", len(aside), len({r["genus"] for r in reps if r["unit_id"] in aside}))
