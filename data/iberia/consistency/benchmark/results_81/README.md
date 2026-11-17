# The experiments carried through the genus stage (#81)

The injection and perturbation experiments of the gate of #56 rebuilt only the backbone. [#81](https://github.com/naturalis/timegraft/issues/81) carries them through the genus stage, grafting and placement, so that their effect on **crown ages** and on the **final tree** (and so on PD) can be measured: metrics 1n, 1o, 5e, 5f, 6e and 6f of gate version 2 (`../GATE.md`, declared in #74 before these runs, all blind).

## Method

`../downstream.py` runs the unchanged workflow scripts `shard_genera.py`, `build_shard.py`, `graft.py` and `place_units.py` on a given backbone, with the settings of the configuration (and `../legacy.yaml` for legacy). `../run_downstream.sh` makes the 14 runs of gate version 2: the clean backbones, and replicates 1 and 2 of the injection, the age perturbation and the subsampling, for both settings. Each is compared with the clean run of the same setting.

The backbones and units are those of the injection rerun of #70 (`../results_70/`, current `main` for the default arm, with the clean default backbone byte-identical to that of #73), and of the #73 benchmark for the age perturbation, the subsampling and the legacy arm, none of which #70, #74 or #84 changed. With `B=data/iberia/consistency/benchmark` and the working folder of `../results_70/`:

```
PATH=/path/to/vsearch:$PATH $B/run_downstream.sh $BENCH70 . 2
PYTHONPATH=src python $B/evaluate.py $BENCH70 out                       # the gate metrics, downstream.tsv
PYTHONPATH=$B python $B/analyse_downstream.py $BENCH70/exp out/injection_genera.tsv out
python $B/results_81/pd_seed_variation.py $BENCH70/exp/downstream pd_seed_variation.tsv
```

**The harness reproduces the full runs.** The clean default run is byte-identical to the default full run of #73 in every shard, `provenance.tsv`, `grafted.tre`, `timegraft.tre` and `units_placed.tsv`. The clean legacy run gives the same 1,183 genera as the full run of `main` at `e20fdf6`, with every column the two `provenance.tsv` share equal, and byte-identical `grafted.tre`, `timegraft.tre` and placements. For that, `legacy.yaml` gains `rate_outlier_action: flag`: `main` did not set long branches aside, and the setting is read in the genus stage only, so the backbone and the screening of the legacy arm are unaffected.

**Time.** A clean run takes 40 to 47 minutes on two cores: 22 to 29 in the genus stage and about 16 in `place_units`. `downstream.py` keeps every `build_genus()` result, keyed on everything the function reads, and reuses it when an experiment calls it with the same inputs. The injection runs then rebuild only the genera whose inputs the errors changed: 2 and 98 of about 1,040 for the default (replicates 1 and 2), 265 for legacy in replicate 1; the default's genus stage of replicate 1 takes 3 seconds instead of 29 minutes. The age perturbation changes every stem age and reuses nothing. Placement is not cached: its time goes into attaching units to the tree. The 14 runs took about 6.3 hours of computing.

## Result

`downstream_metrics.tsv` (the gate metrics) and `downstream.tsv` (per run):

| # | Metric | Legacy | Default | Threshold | Status |
| --- | --- | --- | --- | --- | --- |
| 1n | Injection: corrected genera with a crown age within 1.5 of the clean run | (21 of 42) | 42 of 42 | at least 0.80 | pass |
| 1o | Injection: PD rank correlation with the clean final tree (three designs) | 0.969; 0.952; 0.962 | 1.000; 1.000; 1.000 | not lower than legacy (0.01) | pass |
| 5e | Ages perturbed: 90% quantile of the factor by which crown ages move | 1.479 | 1.364 | at most legacy × 1.05 | pass |
| 5f | Half the genera dropped: the same | 1.677 | 1.653 | at most legacy × 1.05 | pass |
| 6e | Ages perturbed: PD rank correlation with the clean final tree | 0.871; 0.851; 0.874 | 0.934; 0.941; 0.955 | not lower than legacy (0.01) | pass |
| 6f | Half the genera dropped: the same | 0.972; 0.958; 0.935 | 0.961; 0.952; 0.926 | not lower than legacy (0.01) | **fail** (150 units at random: 0.011 below; sampling noise, see below) |

**Injected errors hardly reach the final tree.** Of the 48 genera injected in replicates 1 and 2, 42 have a crown age in both runs (the others are single-unit genera). With the default all 42 are within a factor 1.5 of their crown age in the clean run: 40 unchanged, *Idas* 45.8 against 44.8 Ma (factor 1.02) and *Sphincterochila* 24.4 against 20.1 Ma (1.21). With legacy 21 are within 1.5 (12 unchanged), and the median factor is 1.50 (`downstream_injected.tsv`). The errors also move the crowns of other genera with legacy, 26 and 42 by more than 10%, nearly all through their stems, and with the default 0 and 1 (`downstream_collateral.tsv`). The default's final tree with the errors gives nearly the same PD ranking as without them (1o: rank correlation 0.9997 to 1.0000 per replicate and design).

**Noisy TimeTree ages move crowns less.** With every cached age multiplied by a random factor, the 90% quantile of the crown age changes is 1.34 and 1.40 for the default against 1.47 and 1.48 for legacy, and a crown moves by more than 10% in 46% and 43% of genera against 71% and 49% (`downstream_runs.tsv`). PD rankings hold up better (6e: 0.93 to 0.95 against 0.85 to 0.87). This carries the gain on the backbone (5a, 5b) through to what the pipeline is for.

**Half the genera dropped: crowns as robust; PD rankings no different, though 6f fails.** Crown ages move about as much with either setting (5f). With the gate's sample seed the PD rank correlation is lower for the default in all three designs, by 0.011, 0.006 and 0.009, so 6f fails on the first. That difference is sampling noise (next section).

## Why 6f fails: sampling noise

`evaluate.py` draws the samples of 6f for each setting separately, from the units in that setting's own trees, so legacy and default are judged on different samples. `pd_seed_variation.py` recomputes 6f exactly as `evaluate.py` does for the sample seeds 6 to 25 (`pd_seed_variation.tsv`; the gate uses 6):

| Design | Mean difference (default − legacy) | Standard deviation across seeds | Seeds where the design fails |
| --- | --- | --- | --- |
| 150 units at random | −0.0004 | 0.0069 | 2 of 20 |
| 100 genera, one unit each | +0.0024 | 0.0095 | 2 of 20 |
| 100 genera, five units each | −0.0028 | 0.0107 | 3 of 20 |

Averaged over seeds, the default is as robust to subsampling as legacy in every design. The difference varies with the seed by about as much as the threshold allows, and 6f as declared, with all three designs required to pass, fails for 7 of the 20 seeds, seed 6 among them. 6f as declared cannot tell the two settings apart, and its failure says nothing about the default. 1o and 6e are measured in the same way but pass by margins (0.03 to 0.08) well beyond this noise.

[#104](https://github.com/naturalis/timegraft/issues/104) proposes a gate version 3 that measures 1o, 6e and 6f on paired samples (the same samples for both settings, drawn from units in all the trees compared) with more draws, leaving their thresholds as they are.

**Corrected.** An earlier version of this section explained the failure by the placement of the long branches set aside since #73, from a comparison of PD with and without placed units. That comparison drew different samples too, and its differences are within the same noise; it is withdrawn, and [#103](https://github.com/naturalis/timegraft/issues/103), which built on it, is closed. Most set-aside units are placed within their own genus already (308 of the 313 whose nearest unit is in their genus, in the clean default run).

## What this is not

This is not a run of gate version 2. The downstream metrics above are final for these backbones, but the gate also needs fresh full runs of the current `main` (metrics 0, 2c, 3a to 3d, 6a, 7a to 7e, 8a to 8c) and its experiments built by `run_experiments.sh` from one checkout. Here 7c and 7d cannot be evaluated: they compare the clean backbones of #84 and #70 with runs of #73 that predate their columns.
