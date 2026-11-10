# Iberian run, November 2026

Inputs and outputs of the timegraft runs on genera of western Iberia, described in `docs/science/results/2026-11-iberia.md` and used in the manuscript. Kept in the repository so that the published numbers can be checked and rebuilt.

| File | Contents |
| --- | --- |
| `records.tsv.gz` | Records cut from the BOLD public data package of 27 March 2026 with `timegraft.bcdm` (commands in `docs/technical/running.md`, "The Iberian run"): 102,794 records, 25,166 BINs, 1,232 genera |
| `timetree_cache.json` | TimeTree answers from the runs; copy to `resources/api_cache.json` to get the same backbone ages without querying TimeTree |
| `iberia-run-tolerant.tar.gz` | Output folder of the run with tolerant rooting, the current default and the source of the manuscript's numbers |
| `iberia-run.tar.gz` | Output folder of the run with the earlier strict rooting |
| `rooting-diagnostics.tar.gz` | `diagnostics/rooting/` of the strict run: `per_genus.tsv`, `random_pairs.tsv` (20 random outgroup pairs per genus) and `report.md` (#19) |
| `outgroup-benchmark/` | Scripts and per-genus results of the comparison of outgroup choices and rooting rules (#19) |

The output folders leave out `units/records.tsv` and `units/filtered.tsv`, which the workflow rebuilds from the records.

## How the runs were made

Both runs used `config/iberia.yaml` on top of `config/config.yaml`, on 10 November 2026, with 2 cores, Snakemake 9.27.1, Python 3.13, DendroPy 5.1.3, MAFFT 7.525, VSEARCH 2.30.0, OToL synthesis 16.1 and OTT 3.7.3.

| Run | timegraft commit | Rooting | Time |
| --- | --- | --- | --- |
| tolerant | `f0f55bd` | as in `config/iberia.yaml` | 36 minutes |
| strict | `96375a9` | `rule: strict, n_outgroup: 2, outgroup_choice: sister, fallback: midpoint` | 32 minutes |

To rebuild the tolerant run from the repository root:

```bash
mkdir -p resources/iberia
gunzip -c data/iberia/records.tsv.gz > resources/iberia/records.tsv
cp data/iberia/timetree_cache.json resources/api_cache.json
snakemake --cores 2 --use-conda --configfile config/config.yaml config/iberia.yaml
```

To look at the archived outputs instead:

```bash
mkdir -p results/iberia && tar xzf data/iberia/iberia-run-tolerant.tar.gz -C results/iberia
```

## Main numbers

Both runs: 1,046 genus subtrees plus 137 single-unit genera; 24,882 units with a genus, 215 without (212 placed); 727 records rejected; 24,894 tips in `timegraft.tre`, root at 750 Ma.

- Tolerant: all 1,046 subtrees rooted on their outgroup, 512 cleanly and 534 with units on the wrong side (median 2); 579 units of 351 genera among the outgroup.
- Strict: 852 rooted on the outgroup, 194 fell back to midpoint.
