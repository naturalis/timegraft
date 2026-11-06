"""Choose one representative sequence per unit and assign its genus."""

import logging
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

from timegraft.seqio import read_records, write_tsv
from timegraft.units import REPRESENTATIVE_COLUMNS, UNASSIGNED_COLUMNS, representatives

records = read_records(snakemake.input[0])
reps, unassigned = representatives(
    records,
    max_records=snakemake.params.max_records,
    genus_agreement=snakemake.params.genus_agreement,
)
write_tsv(snakemake.output.reps, reps, REPRESENTATIVE_COLUMNS)
write_tsv(snakemake.output.unassigned, unassigned, UNASSIGNED_COLUMNS)
logging.info("%d units assigned to a genus, %d unassigned", len(reps), len(unassigned))
