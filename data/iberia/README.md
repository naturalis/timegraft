# Iberian run, November 2026

Inputs and outputs of the timegraft runs on genera of western Iberia, described in `docs/science/results/2026-11-iberia.md` and used in the manuscript. Kept in the repository so that the published numbers can be checked and rebuilt.

| File | Contents |
| --- | --- |
| `records.tsv.gz` | Records cut from the BOLD public data package of 27 March 2026 with `timegraft.bcdm` (commands in `docs/technical/running.md`, "The Iberian run"): 102,794 records, 25,166 BINs, 1,232 genera |
| `timetree_cache.json` | TimeTree answers from the runs; copy to `resources/api_cache.json` to get the same backbone ages without querying TimeTree |
| `iberia-run-tolerant.tar.gz` | Output folder of the run with tolerant rooting, bounded crown ages and calibrated rates: the current settings of `config/iberia.yaml` and the source of the manuscript's numbers |
| `calibration/` | The substitution rates per class used by the run, and how they were derived (#18) |
| `dating-benchmark/` | Tests of penalised likelihood, rates and bounds on the run as it was before crown ages changed (#36) |
| `iberia-run.tar.gz` | Output folder of the run with the earlier strict rooting (and the earlier crown age rule) |
| `rooting-diagnostics.tar.gz` | `diagnostics/rooting/` of the strict run: `per_genus.tsv`, `random_pairs.tsv` (20 random outgroup pairs per genus) and `report.md` (#19) |
| `outgroup-benchmark/` | Scripts and per-genus results of the comparison of outgroup choices and rooting rules (#19) |

The output folders leave out `units/records.tsv` and `units/filtered.tsv`, which the workflow rebuilds from the records.

## How the runs were made

Both runs used `config/iberia.yaml` on top of `config/config.yaml`, on 10 November 2026, with 2 cores, Snakemake 9.27.1, Python 3.13, DendroPy 5.1.3, MAFFT 7.525, VSEARCH 2.30.0, OToL synthesis 16.1 and OTT 3.7.3.

| Run | timegraft commit | Rooting | Time |
| --- | --- | --- | --- |
| tolerant | `f0f55bd`, rerun with `df16e3d`; genus subtrees onwards rerun on 11 November with the crown age changes of #36 | as in `config/iberia.yaml` | 36 minutes; rerun 28; 36 without the backbone |
| strict | `96375a9` | `rule: strict, n_outgroup: 2, outgroup_choice: sister, fallback: midpoint` | 32 minutes |

The tolerant run was rerun the same evening with `df16e3d`, which adds the root uncertainty columns to `provenance.tsv` ([#34](https://github.com/naturalis/timegraft/issues/34)), and gave otherwise identical outputs.

On 11 November 2026 it was rerun from the record filter onwards with the crown age changes of [#36](https://github.com/naturalis/timegraft/issues/36): bounds in place of the TimeTree threshold, rates per class (`calibration/`) and the literature table (`data/crown_ages/literature.tsv`). The backbone of the earlier run was kept as it was (the OToL files were not at hand), and the rates were calibrated on the earlier run's genus trees, which this run reproduces: calibrating on the archived run gives the same table. Tools: Snakemake 9.27.1, Python 3.13, DendroPy 5.1.3, MAFFT 7.505 and VSEARCH 2.27.0 on the path. Units, backbone and the rooting of all 1,046 genus trees are the same as before; crown ages, and with them the trees, are not. The archive now holds this run. The run before it is in the repository history (commit `ea8edb3`), and its `provenance.tsv` is kept as `dating-benchmark/provenance_before.tsv`.

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

- Crown ages in the tolerant run: 780 from the barcodes within all bounds, 134 lowered to the TimeTree age, 129 lowered to 0.9 of the stem age, 3 from the literature table; median 24.8 Ma. Total branch length of the tree 552,322 Ma (312,011 Ma under the earlier crown age rule).

- Tolerant: all 1,046 subtrees rooted on their outgroup, 512 cleanly and 534 with units on the wrong side (median 2); 579 units of 351 genera among the outgroup.
- Strict: 852 rooted on the outgroup, 194 fell back to midpoint.
