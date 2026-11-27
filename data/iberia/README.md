# Iberian run, November 2026

Inputs and outputs of the timegraft runs on genera of western Iberia, described in `docs/science/results/2026-11-iberia.md` and used in the manuscript. Kept in the repository so that the published numbers can be checked and rebuilt.

| File | Contents |
| --- | --- |
| `records.tsv.gz` | Records cut from the BOLD public data package of 27 March 2026 with `timegraft.bcdm` (commands in `docs/technical/running.md`, "The Iberian run"): 102,794 records, 25,166 BINs, 1,232 genera |
| `timetree_cache.json` | The TimeTree and ChecklistBank answers of the reference run (23,563: 11,119 TimeTree, 12,444 ChecklistBank for COL26.10 XR; cache format version 2, #47): everything the run asks, including the NCBI lineages and probe queries of the node check (#108), so that copying it to `resources/api_cache.json` gives the reference run without querying either service. Until #77 it held only the 9,820 answers of the runs before #101 and #108, all of which it still holds unchanged |
| `iberia-reference.tar.gz` | **The reference run**: output folder of the full run of the current pipeline at its default settings (#77; see "The reference run" below) |
| `iberia-run-tolerant.tar.gz` | Output folder of the run with tolerant rooting, unit screening, a pure-birth backbone and crown ages from TimeTree or the midpoint of their bounds, made on 12 November 2026 **before the consistency layer** (#46): the source of the manuscript's numbers, but superseded by `iberia-reference.tar.gz` for the pipeline as it now is. `units/screened_representatives.tsv` is left out (it is `units/representatives.tsv` without the units set aside in `units/screened.tsv`) |
| `screening/` | The analysis of units that do not fit their genus, and scripts comparing runs with and without screening (#38) |
| `calibration/` | Substitution rates per class fitted to these genera; tried, rejected and not used (#18) |
| `dating-benchmark/` | Tests of penalised likelihood, rates and bounds on the run as it was before crown ages changed (#36) |
| `iberia-run.tar.gz` | Output folder of the run with the earlier strict rooting (and the earlier crown age rule); superseded like the tolerant run |
| `rooting-diagnostics.tar.gz` | `diagnostics/rooting/` of the strict run: `per_genus.tsv`, `random_pairs.tsv` (20 random outgroup pairs per genus) and `report.md` (#19) |
| `consistency/` | Evidence for the consistency layer (#46) and its follow-ups (#89): a section per issue, each measured on the state of its own stage, and `benchmark/`, the release gate with the full runs of legacy and of the current pipeline and their comparison |
| `outgroup-benchmark/` | Scripts and per-genus results of the comparison of outgroup choices and rooting rules (#19) |
| `allgenus_nocf_13nov26.txt` | List of genera from Jorge Moutinho's checklist. Currently unused but may be input for later Iberian runs after refactoring |

The output folders leave out `units/records.tsv` and `units/filtered.tsv`, which the workflow rebuilds from the records.

## The reference run

`iberia-reference.tar.gz` is the full run of the pipeline as it now is, at its default settings: the full default run of the release gate of #108 (`consistency/benchmark/results_108/gate/`, where all 41 must metrics pass). It replaces the archives above as the reference for the current pipeline; they stay as the record of the analyses made on them.

- **Commit**: `bf58ca6` (branch of #124: `main` at `53a0630` with `backbone_build.ages.node_check.mode: report`). Between `53a0630` and this commit no file of `src/` or `workflow/` changed, and of `config/` only that default; `environment.yml` was added and the package description in `pyproject.toml` reworded.
- **Configuration**: `config/config.yaml` with `config/iberia.yaml`, unchanged.
- **Inputs**: `records.tsv.gz`; OToL synthesis 16.1 (`labelled_supertree.tre`, `annotations.json`) and OTT 3.7.3; the TimeTree and ChecklistBank answers in `timetree_cache.json`, which the run used as its cache and to which it added none.
- **Made** on 26 November 2026, on Linux with 2 cores: Snakemake 9.27.1, Python 3.13.18, DendroPy 5.1.3, MAFFT 7.505, VSEARCH 2.32.0. 47 minutes, peak memory 1,072 MB (`consistency/benchmark/results_108/gate/runtime.tsv`).
- **Left out of the archive**: `units/filtered.tsv`, which the workflow rebuilds from the records, as in the other archives.

Its main numbers (`consistency/benchmark/results_108/gate/end_to_end.tsv`, which also compares them with the legacy pipeline):

- 102,794 records, 727 rejected; 1,186 genera in the backbone (784 nodes), root at 750 Ma;
- 1,039 genus subtrees and 147 single-unit genera, with 23,003 units;
- 1,136 units set aside by screening and 432 rate outliers set aside in the genus stage (#73), with 215 units without a usable genus: 1,253 of these placed next to their nearest unit, 521 on the stem of their genus, 9 not placed;
- 24,924 tips in `timegraft.tre`, total branch length 643,604 Ma;
- crown ages: 262 from TimeTree (16 more demoted to a bound), 748 at the midpoint of their bounds, 10 at the barcode cap (#79), median 24.96 Ma;
- 346 genera with units on the wrong side of the root, 253 units of a genus among its outgroup.

To look at it:

```bash
mkdir -p results/iberia && tar xzf data/iberia/iberia-reference.tar.gz -C results/iberia
```

To make it again, run the commands in "How the runs were made" below on `main`; with the cache copied, no question goes to TimeTree or ChecklistBank.

## How the runs were made

This section describes the two earlier archives, `iberia-run-tolerant.tar.gz` and `iberia-run.tar.gz`; the reference run is described above.


Both runs used `config/iberia.yaml` on top of `config/config.yaml`, on 10 November 2026, with 2 cores, Snakemake 9.27.1, Python 3.13, DendroPy 5.1.3, MAFFT 7.525, VSEARCH 2.30.0, OToL synthesis 16.1 and OTT 3.7.3.

| Run | timegraft commit | Rooting | Time |
| --- | --- | --- | --- |
| tolerant | `f0f55bd`, rerun with `df16e3d`; genus subtrees onwards rerun on 11 November with the crown age changes of #36 | as in `config/iberia.yaml` | 36 minutes; rerun 28; 27 for genus subtrees onwards |
| strict | `96375a9` | `rule: strict, n_outgroup: 2, outgroup_choice: sister, fallback: midpoint` | 32 minutes |

The tolerant run was rerun the same evening with `df16e3d`, which adds the root uncertainty columns to `provenance.tsv` ([#34](https://github.com/naturalis/timegraft/issues/34)), and gave otherwise identical outputs.

On 11 November 2026 it was rerun with the crown age changes of [#36](https://github.com/naturalis/timegraft/issues/36): crown ages within recorded bounds in place of the TimeTree threshold, and the literature table (`data/crown_ages/literature.tsv`). The record filter and representatives were rerun and gave identical units; the genus subtrees, grafting and placement were then rebuilt (27 minutes). The backbone of the earlier run was kept as it was (the OToL files were not at hand). Tools: Snakemake 9.27.1, Python 3.13, DendroPy 5.1.3, MAFFT 7.505 and VSEARCH 2.27.0 on the path. Units, backbone and the rooting of all 1,046 genus trees are the same as before. The archive now holds this run. The run before it is in the repository history (commit `ea8edb3`), and its `provenance.tsv` is kept as `dating-benchmark/provenance_before.tsv`. A run in between, with rates fitted per class, was rejected and is not archived (`calibration/README.md`).

Later the same day it was rerun from screening onwards with the unit screening of [#38](https://github.com/naturalis/timegraft/issues/38) (31 minutes): 1,294 units were set aside and placed separately. The archive now holds this run; the run before it is at the commit before the screening change.

On 12 November 2026 it was rebuilt from the backbone onwards (45 minutes), with the pure-birth interpolation of undated backbone nodes ([#16](https://github.com/naturalis/timegraft/issues/16)) from the Open Tree release files, and crown ages from TimeTree where available and the midpoint of their bounds otherwise ([#42](https://github.com/naturalis/timegraft/issues/42)). TimeTree answers came from `timetree_cache.json`; all 314 TimeTree ages are as before. The archive now holds this run.

To run the pipeline on these records from the repository root (on `main` this gives the reference run, not the earlier archives, which were made with the commits listed above):

```bash
mkdir -p resources/iberia
gunzip -c data/iberia/records.tsv.gz > resources/iberia/records.tsv
cp data/iberia/timetree_cache.json resources/api_cache.json
snakemake --cores 2 --use-conda --configfile config/config.yaml config/iberia.yaml
```

To look at the earlier tolerant run instead:

```bash
mkdir -p results/iberia && tar xzf data/iberia/iberia-run-tolerant.tar.gz -C results/iberia
```

## Main numbers of the earlier archived runs

These describe the runs of 10 to 12 November, before the consistency layer; for the current pipeline, see "The reference run" above.


Both runs as first made, before screening: 1,046 genus subtrees plus 137 single-unit genera; 24,882 units with a genus, 215 without (212 placed); 727 records rejected; 24,894 tips in `timegraft.tre`, root at 750 Ma.
- Tolerant run, as archived: 24,909 tips; 1,036 genus subtrees and 147 single-unit genera; 1,294 units set aside by screening (819 placed next to their nearest unit, 474 on the stem of their genus, 1 unplaced); 215 units without a usable genus (210 placed). Crown ages: 261 from TimeTree, 772 at the midpoint of their bounds, 1 at 0.9 of the stem age, 2 from the literature table (cut to 0.9 of the stem age); median 24.8 Ma. Total branch length 638,002 Ma (360,668 Ma with every crown at its lower bound, 1,424,818 Ma at its upper bound; 312,011 Ma under the first crown age rule).

- Tolerant, before screening: all 1,046 subtrees rooted on their outgroup, 512 cleanly and 534 with units on the wrong side (median 2); 579 units of 351 genera among the outgroup. With screening: 1,036 subtrees, 670 cleanly and 366 with units on the wrong side; 276 units of 159 genera among the outgroup.
- Strict: 852 rooted on the outgroup, 194 fell back to midpoint.
