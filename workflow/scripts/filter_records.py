"""Apply the marker's filter plug-in to the record table."""

import logging
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

from timegraft.filters import get_check
from timegraft.seqio import read_records, write_tsv

records = read_records(snakemake.input[0])
check = get_check(snakemake.params.marker)
kept, rejected = [], []
for record in records:
    reason = check(record["sequence"], **snakemake.params.options)
    if reason is None:
        kept.append(record)
    else:
        rejected.append({**record, "reason": reason})

columns = list(records[0].keys()) if records else ["record_id", "unit_id", "genus", "sequence"]
write_tsv(snakemake.output.kept, kept, columns)
write_tsv(snakemake.output.rejected, rejected, ["record_id", "unit_id", "genus", "reason"])
logging.info("kept %d of %d records; rejected %d", len(kept), len(records), len(rejected))
