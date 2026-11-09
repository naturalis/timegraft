"""Apply the marker's filter plug-in to the record table.

For markers with a reading-frame check (COI), the genetic code is resolved per
record from its taxonomy, and written to the output with where it came from.
"""

import logging
import sys
from collections import Counter

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

from timegraft.filters import USES_GENETIC_CODE, get_check
from timegraft.geneticcode import GeneticCodes
from timegraft.seqio import read_records, write_tsv

marker = snakemake.params.marker
options = dict(snakemake.params.options)
overrides = options.pop("genetic_codes", None) or {}
table = options.pop("genetic_code_table", None) or None
default_code = options.pop("genetic_code", 5)

records = read_records(snakemake.input[0])
check = get_check(marker)
codes = GeneticCodes(default_code, overrides, table) if marker in USES_GENETIC_CODE else None

kept, rejected = [], []
tally = Counter()  # (code, source, outcome)
for record in records:
    extra = {}
    if codes is not None:
        code, source, _ = codes.resolve(record)
        record["genetic_code"], record["genetic_code_source"] = str(code), source
        extra["genetic_code"] = code
    reason = check(record["sequence"], **options, **extra)
    if reason is None:
        kept.append(record)
    else:
        rejected.append({**record, "reason": reason})
    if codes is not None:
        tally[(record["genetic_code"], record["genetic_code_source"], reason or "kept")] += 1

columns = list(records[0].keys()) if records else ["record_id", "unit_id", "genus", "sequence"]
reject_columns = ["record_id", "unit_id", "genus", "reason"]
if codes is not None:
    reject_columns += ["genetic_code", "genetic_code_source"]
write_tsv(snakemake.output.kept, kept, columns)
write_tsv(snakemake.output.rejected, rejected, reject_columns)
logging.info("kept %d of %d records; rejected %d", len(kept), len(records), len(rejected))

if codes is not None:
    # A code whose stop-codon rejection rate stands out points to a wrong mapping.
    logging.info("genetic code\tsource\trecords\tstop_codon\trate")
    for code, source in sorted({(c, s) for c, s, _ in tally}, key=lambda k: (int(k[0]), k[1])):
        n = sum(v for (c, s, _), v in tally.items() if (c, s) == (code, source))
        stops = tally[(code, source, "stop_codon")]
        logging.info("%s\t%s\t%d\t%d\t%.3f", code, source, n, stops, stops / n)
    n_default = sum(v for (_, s, _), v in tally.items() if s == "default")
    if n_default:
        logging.warning("%d records matched no taxon and got the default genetic code %d", n_default, default_code)
