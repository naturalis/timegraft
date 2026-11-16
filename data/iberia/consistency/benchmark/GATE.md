# Release gate for the consistency layer: metrics and thresholds

Declared on 13 November 2026 for [#56](https://github.com/naturalis/timegraft/issues/56), before the benchmark was run. This file is committed on its own, ahead of the harness and the results, so that the history shows the order. The results and the status of every metric are in `REPORT.md`.

## What is compared

- **legacy**: the pipeline as on `main` (commit `e20fdf6`). Run once end to end from `main` itself. For the experiments that need many backbone builds, `consistency-layer` is run with the settings that switch the layer off (`legacy.yaml`); the gate first checks that this reproduces the backbone, screening and crown table of the run from `main` (metric 0).
- **default**: `consistency-layer` at its defaults.

Both on the Iberian records (`data/iberia/records.tsv.gz`), `config/iberia.yaml`, OTT 3.7.3, synthesis 16.1, COL release 316691, and the TimeTree cache of each branch. Every random choice has a fixed seed in the harness.

## What was already known when this was written

The issues of the epic each reported their own stage on the Iberian run, so some numbers were known before these thresholds were set: the conflict metrics of the age fit (#52), its hold-out error, the crown tails (#54), and the screening injections (#55). For those the thresholds below cannot count as blind; they are marked **known**. Not known in advance: anything end to end (no full run of the branch existed), the placement injections, both perturbation experiments, the diversity comparisons, runtime, memory and reproducibility. Those are marked **blind**.

## Decision rule

- **Ship as default**: every *must* metric passes.
- **Ship as an option** (merge, but with the legacy settings as defaults): no *must* metric of families 5 to 7 fails (nothing gets less robust, less stable or operationally worse), but one or more *must* metrics of families 1 to 4 fail (the layer does not deliver what it was built for).
- **Hold**: a *must* metric of families 5 to 7 fails, or metric 0 fails (the comparison itself is not sound).

*Report* metrics have no threshold; they are there to explain the others. A threshold of the form "not worse" allows the stated tolerance, because both arms carry noise.

## Metrics

| # | Family | Metric | Threshold | Kind | Known or blind |
| --- | --- | --- | --- | --- | --- |
| 0 | Soundness | `consistency-layer` with `legacy.yaml` gives the same `backbone.tre`, backbone node ages, set-aside units and TimeTree crown ages as the run from `main` | identical | must | blind |
| 1a | Known errors | *Onchidella* (a gastropod that OTT merges with a sea urchin homonym) in the final tree: the other genera in the smallest backbone clade that contains it | all Mollusca | must | known at backbone level, blind end to end |
| 1b | Known errors | Injected misplacements (see below), genera with at least two checks that can be evaluated: share flagged | at least 0.90 | must | blind |
| 1c | Known errors | Injected misplacements, all: share flagged | at least 0.50 | must | blind |
| 1d | Known errors | Flagged injected genera: share placed back in their own class | at least 0.90 | must | blind |
| 1e | Known errors | Genera without an injected error that are flagged, per replicate (not counting *Onchidella*) | at most 1 on average | must | blind |
| 1f | Known errors | Corrected injected genera: stem age within a factor 1.5 of their stem age without the error | at least 0.80 | must | blind |
| 1g | Known errors | Genera without an injected error whose stem age moves by more than 10% because of the errors: default against legacy | at most half | must | blind |
| 2a | Backbone conflicts | TimeTree calibrations off from the final age by more than a factor 1.5 | at most half of legacy | must | known |
| 2b | Backbone conflicts | Root mean square of calibration residuals over their uncertainty | lower than legacy | must | known |
| 2c | Backbone conflicts | Nodes older than their parent in the final backbone | 0 | must | blind |
| 3a | Crown tails | Genera whose crown age is more than 20 times the age from the barcodes, in the final provenance | 0 | must | known for the genus stage of #54, blind end to end |
| 3b | Crown tails | The same above 10 times | fewer than legacy | must | as 3a |
| 3c | Crown tails | Genera with a TimeTree crown used that is less than a third of the barcode age | 0 | must | as 3a |
| 3d | Crown tails | Genera with a mean root-to-tip distance below 0.02 and a crown age above 20 Ma | fewer than legacy | report | blind |
| 4a | Hold-out | Median error factor on held-out TimeTree node ages (50 random halves), default against legacy | not worse (within 2%) | must | known |
| 4b | Hold-out | The same at the 90% quantile | not worse (within 2%) | must | known |
| 5a | Perturbation | TimeTree ages perturbed (each cached age times a log-normal factor, sigma 0.3; 5 replicates): 90% quantile of the factor by which backbone node ages move | not worse than legacy (within 5%) | must | blind |
| 5b | Perturbation | The same: share of node ages that move by more than a factor 1.5 | not worse than legacy (within 1 percentage point) | must | blind |
| 5c | Perturbation | Half the genera dropped (5 replicates): 90% quantile of the factor by which the ages of the remaining clades move | not worse than legacy (within 5%) | must | blind |
| 5d | Perturbation | The same for genus stem ages | not worse than legacy (within 5%) | must | blind |
| 6a | Diversity | Faith's PD of the same random samples on the final trees of legacy and default (150 units; 100 genera, one unit; 100 genera, five units; 200 samples each): rank correlation | at least 0.90 in each design | must | blind |
| 6b | Diversity | Genera whose stem age differs by more than 10% between legacy and default: share with a recorded cause (new anchor, re-anchoring, another identifier sent to TimeTree, or a fitted age that departs from the precedence age on its path to the root) | at least 0.90 | must | blind |
| 6c | Diversity | PD of genus samples (100 genera, 200 samples) on the backbone: rank correlation between the unperturbed backbone and those of 5a, default against legacy | not lower (within 0.01) | must | blind |
| 6d | Diversity | PD ratio default over legacy per design, and the change in total tree length | none | report | blind |
| 7a | Operations | Wall time of the full run, two cores, warm cache | at most 1.5 times legacy | must | blind |
| 7b | Operations | Peak memory of the full run | at most 1.5 times legacy, and below 7 GB | must | blind |
| 7c | Operations | Backbone and screening rebuilt from the same inputs and cache | byte-identical outputs | must | blind |
| 7d | Operations | Backbone rebuilt with the network blocked, from the cache of the run | succeeds, identical outputs | must | blind |
| 7e | Operations | Test suite | passes | must | known |

## The injected misplacements (metrics 1b to 1g)

The one real misplacement in the run, *Onchidella*, came about because OTT attaches the genus name to another animal. The injection reproduces exactly that in the input, and leaves the pipeline untouched: in a copy of the OTT taxonomy, a genus of the run is merged into a genus-rank taxon elsewhere in the animals that is a node of the synthesis and not in the run. That taxon takes the genus's name, its NCBI ids (after its own, as for *Onchidella*, whose OTT entry lists the sea urchin's id first) and its species as children; the genus's own taxon is left empty under another name. The pipeline then matches the genus to the wrong taxon, as it did for *Onchidella*. (The NCBI ids and the species were added to this description after the first commit of this file and before any injection was run, to make the injected error the same kind of merge.)

Five replicates. In each, 24 genera drawn at random from the genera with a backbone position in the unperturbed default run: 12 moved to a taxon of another phylum, 12 to a taxon of another class of their own phylum. The backbone is then built with the legacy and the default settings. "Its own class" is that of the genus's records. A check "can be evaluated" as `concordance.tsv` reports it.

Misplacements within a class are outside the default `flag_level: class` and are not injected.

## The perturbations (metrics 5a to 5d, 6c)

- **Ages.** In a copy of the TimeTree cache, every cached age is multiplied by exp(0.3 z), z standard normal, seeded per replicate. The backbone is rebuilt offline with each setting, and every node age is compared with the unperturbed backbone of the same setting.
- **Genera.** Half the genera of the run are dropped at random (their records and units), seeded per replicate; the backbone is rebuilt with each setting (TimeTree is asked for the new picks), and each clade of the smaller backbone is compared with the most recent common ancestor of the same genera in the full backbone of the same setting.

Both compare each setting with itself, so they measure how much the result depends on the noise, not which result is right.

## Version 2

Declared on 16 November 2026 for [#74](https://github.com/naturalis/timegraft/issues/74), before any run that it judges, and committed on its own like version 1. Version 1 above stays as it was, so that the report of #56 and the reruns of #89 and #73 can be reproduced (`evaluate.py --version 1`). Version 2 is the gate for later changes to the layer (`evaluate.py --version 2`, the default).

### Why a second version

Running version 1 showed where its metrics were weak (`REPORT.md`, "What to do", point 6):

1. **6b was too easy.** It counted a changed stem as explained when *any* node on its path to the root had a fitted age that departed from the precedence age, which held for all 175 changed stems. [#84](https://github.com/naturalis/timegraft/issues/84) records where each node's age comes from, and `compare_backbones.py` follows each changed stem to the nodes that moved it.
2. **7b as a ratio meant little at this scale.** It failed at 1.32 GB against an allowed 1.23 GB while its absolute limit was 7 GB.
3. **Rooting was not measured.** The screening fix of #55 put more units on the wrong side of the root ([#75](https://github.com/naturalis/timegraft/issues/75)), and no metric saw it.
4. **The two injection variants were not metrics.** `inject_own` and `inject_name` explained 1f, 1g and the detection limit, but were added after the first results.
5. **Everything stopped at the backbone.** The injection and perturbation experiments did not reach crown ages or PD ([#81](https://github.com/naturalis/timegraft/issues/81)).

### What is compared

As in version 1: **legacy** is the pipeline without the layer (`main` at `e20fdf6` for the full run; `legacy.yaml` on the branch under test for the experiments, checked by metric 0), and **default** is the branch under test at its defaults. The full run of legacy can be reused from an earlier gate as long as its inputs are the same (records, configuration, OTT, synthesis, COL release and cache), as the reruns of #89 and #73 did. Metric 6b compares the two clean backbones of the harness (`legacy_clean` and `default_clean`), not the backbones of the full runs: the run from `e20fdf6` predates the columns of #84, and metric 0 and 7c check that the clean backbones equal those of the full runs.

### Decision rule

As in version 1, with two changes:

- **Families 5 to 8** guard robustness, stability, operations and rooting: a failing *must* metric of these families, or of metric 0, means **hold**. A failing *must* metric of families 1 to 4 alone means **ship as an option**. All *must* metrics passing means **ship as default**.
- **Metrics that were not run.** The downstream metrics (1n, 1o, 5e, 5f, 6e, 6f) need the harness of #81. Where their outputs do not exist, `evaluate.py` reports them as *not run*, and the decision is reached on the other metrics and says that it was. Once the harness exists, every run of the gate includes them, and a *must* metric that was not run counts as failed.

### Metrics

Metrics not listed here are as in version 1, with the same threshold and kind. For every metric of version 1, the result for the current default is known from the rerun of #73 (`results_73/`): of the thresholds below, only those of the downstream metrics are blind.

| # | Family | Metric | Threshold | Kind | Known or blind |
| --- | --- | --- | --- | --- | --- |
| 1a | Known errors | As in version 1, stated as evaluated since #73: the units of *Onchidella*'s subtree, leaving out units set aside by screening or in the genus stage (`subtrees/set_aside_*.tsv`), which are placed by nearest neighbour and not taken to belong to the genus | all Mollusca | must | known |
| 1h | Known errors | `inject_own` (the merged taxon holds only the genus's own species): share flagged | at least 0.90 | must | known (1.00) |
| 1i | Known errors | `inject_own`: flagged genera back in their own class | at least 0.90 | must | known (0.99) |
| 1j | Known errors | `inject_own`: corrected genera with a stem age within a factor 1.5 of their stem age without the error | at least 0.80 | must | known (0.82) |
| 1k | Known errors | `inject_name` (only the name is merged: no NCBI ids, no species): share flagged | at least 0.50 | must | known (1.00) |
| 1l | Known errors | `inject_name`: flagged genera back in their own class | at least 0.90 | must | known (1.00) |
| 1m | Known errors | `inject_name`: corrected genera with a stem age within a factor 1.5 of their stem age without the error | none | report | known (0.72) |
| 1n | Known errors | Injection, carried through the genus stage: corrected injected genera with a crown age within a factor 1.5 of their crown age in the clean run of the same setting | at least 0.80 | must | blind |
| 1o | Known errors | Injection, carried through: PD of the designs of 6a on the final tree with the errors against the clean final tree of the same setting, rank correlation, default against legacy | not lower (within 0.01) | must | blind |
| 5e | Perturbation | Ages perturbed, carried through: 90% quantile of the factor by which crown ages move against the clean run of the same setting | not worse than legacy (within 5%) | must | blind |
| 5f | Perturbation | Half the genera dropped, carried through: the same for the crown ages of the remaining genera | not worse than legacy (within 5%) | must | blind |
| 6b | Diversity | Genus stems whose age differs by more than 10% between the clean legacy and default backbones: share traced to a cause by `compare_backbones.py` (a cause other than a calibration gained or lost: see its documentation) | at least 0.95 | must | known (1.00) |
| 6e | Diversity | Ages perturbed, carried through: PD of the designs of 6a on the final trees, rank correlation with the clean final tree of the same setting, default against legacy | not lower (within 0.01) | must | blind |
| 6f | Diversity | Half the genera dropped, carried through: the same, for samples drawn from the remaining genera | not lower (within 0.01) | must | blind |
| 7b | Operations | Peak memory of the full run (the process tree, as `monitor.py` records it) | below 2,000 MB | must | known (1,045 MB) |
| 7f | Operations | Peak memory of the full run, default over legacy | none | report | known (1.26) |
| 8a | Rooting | Genera with units on the wrong side of the root (`violations` above 0 in `provenance.tsv`) | not more than legacy (within 5%) | must | known (346 against 366) |
| 8b | Rooting | Units of a genus that the root puts among its outgroup (`misplaced` in `provenance.tsv`) | not more than legacy (within 5%) | must | known (248 against 276) |
| 8c | Rooting | Crowns that depend strongly on the root: `crown_age_max` at least 1.25 times `crown_age_min` in `provenance.tsv` | none | report | known (11 against 20) |

Why these thresholds:

- **1h to 1l** carry over the thresholds of the injection of the gate (1b or 1c, 1d, 1f) unchanged. 1h takes the 0.90 of 1b because nearly every `inject_own` genus has two checks that can be evaluated (118 of 120), as in the injection of the gate; 1k takes the 0.50 of 1c because only 69 of 120 `inject_name` genera do.
- **1m is a report.** Its value (0.72 in the rerun of #73) is below the 0.80 of 1f, and a failing *must* metric of family 1 would make the verdict "ship as an option" for the whole layer on the weakest of three variants. It is reported so that [#70](https://github.com/naturalis/timegraft/issues/70) can show what it changes.
- **6b** at 0.95 rather than 1.00 leaves room for chains that `compare_backbones.py` cannot follow (it stops after twelve steps), without letting an untraced share of more than one stem in twenty pass.
- **7b** at 2,000 MB is about twice the peak of the rerun of #73 and well within what the pipeline is meant to run on; the ratio is kept as a report (7f).
- **8a and 8b** allow 5% because both counts move with the screening and the placement of units, which both settings share only in part.
- **The downstream metrics** take the thresholds of their backbone counterparts: 1n of 1f, 5e and 5f of 5a and 5c, 1o, 6e and 6f of 6c.

### The downstream experiments (1n, 1o, 5e, 5f, 6e, 6f)

The harness of #81 runs the genus stage, grafting and placement on the backbones of two replicates (seeds 1 and 2) of the injection, the age perturbation and the subsampling, for both settings, and once on the clean backbones. `evaluate.py` reads its outputs from `exp/downstream/`:

```
exp/downstream/clean/{legacy,default}/           provenance.tsv, timegraft.tre, units/representatives.tsv
exp/downstream/{inject,noise,sub}/s{1,2}/{legacy,default}/   the same
```

Each experiment is compared with the clean run of the same setting, as the backbone experiments are. Crown ages are compared for genera with a subtree in both; for the injection, over the corrected injected genera (flagged and back in their own class on the backbone). PD uses the three designs of 6a (seed 6), drawn from the units in both final trees, and for 6f from the units of the remaining genera. The correlations of 1o, 6e and 6f are the means over the replicates.
