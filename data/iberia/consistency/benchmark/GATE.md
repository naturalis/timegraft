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

The one real misplacement in the run, *Onchidella*, came about because OTT attaches the genus name to another animal. The injection reproduces exactly that in the input, and leaves the pipeline untouched: in a copy of the OTT taxonomy, the name of a genus of the run is moved to a genus-rank taxon elsewhere in the animals that is a node of the synthesis and not in the run, and the genus's own taxon gets another name. The pipeline then matches the genus to the wrong taxon, as it did for *Onchidella*.

Five replicates. In each, 24 genera drawn at random from the genera with a backbone position in the unperturbed default run: 12 moved to a taxon of another phylum, 12 to a taxon of another class of their own phylum. The backbone is then built with the legacy and the default settings. "Its own class" is that of the genus's records. A check "can be evaluated" as `concordance.tsv` reports it.

Misplacements within a class are outside the default `flag_level: class` and are not injected.

## The perturbations (metrics 5a to 5d, 6c)

- **Ages.** In a copy of the TimeTree cache, every cached age is multiplied by exp(0.3 z), z standard normal, seeded per replicate. The backbone is rebuilt offline with each setting, and every node age is compared with the unperturbed backbone of the same setting.
- **Genera.** Half the genera of the run are dropped at random (their records and units), seeded per replicate; the backbone is rebuilt with each setting (TimeTree is asked for the new picks), and each clade of the smaller backbone is compared with the most recent common ancestor of the same genera in the full backbone of the same setting.

Both compare each setting with itself, so they measure how much the result depends on the noise, not which result is right.
