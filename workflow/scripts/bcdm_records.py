"""Filter a BOLD data package (BCDM TSV) into a record table."""

import logging
import sys

sys.path.insert(0, snakemake.params.src)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")

from timegraft.bcdm import convert

taxa = {rank: set(names if isinstance(names, list) else [names]) for rank, names in snakemake.params.taxa.items()}
counts = convert(snakemake.input[0], snakemake.output[0], marker=snakemake.params.marker, taxa=taxa,
                 countries=set(snakemake.params.countries) or None)
logging.info("%s", counts)
