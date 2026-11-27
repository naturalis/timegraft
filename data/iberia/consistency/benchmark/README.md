# Benchmark and release gate of the consistency layer (#56)

Whether `consistency-layer` should become the default: the metrics and thresholds, the harness that measures them on the Iberian run, and the outcome.

| File | Contents |
| --- | --- |
| `GATE.md` | The metrics, their thresholds and the decision rule, committed before anything was run: version 1 for #56, version 2 (#74), version 3 for later changes (#104; see [below](#three-versions-of-the-gate)) |
| `REPORT.md` | The outcome: status of every metric, what the numbers mean, and the recommendation |
| `results/metrics.tsv` | The gate as a table: metric, legacy, default, threshold, status |
| `results/*.tsv` | One table per experiment (see below) |
| `legacy.yaml` | Settings that switch the layer off: `consistency-layer` run as `main` ran |
| `backbone.py`, `screen.py` | Run the unchanged workflow scripts for the backbone and the screening outside Snakemake, with settings from the configuration files plus overlays |
| `downstream.py`, `run_downstream.sh` | The genus stage, grafting and placement on a backbone outside Snakemake, with the unchanged workflow scripts, and the runs of the downstream metrics of version 2 (#81; see [below](#the-downstream-runs)) |
| `inject.py` | Merges genera of the run into taxa elsewhere in a copy of the OTT taxonomy, as OTT did with *Onchidella* |
| `perturb_cache.py` | Multiplies every TimeTree age in a copy of the cache by a random factor |
| `timetree_cache.json.gz`, `merge_caches.py` | The benchmark cache: every TimeTree and ChecklistBank answer of the benchmark builds and full runs (#87; see [below](#seeds-and-inputs)), and the script that merges caches into it |
| `subsample.py` | Keeps a random half of the genera |
| `run_experiments.sh` | All backbone builds of the gate: clean, injected (with both variants since version 2), perturbed, subsampled, ablations, reproducibility. Exits non-zero if a build fails; rerunning it retries only the failed builds |
| `run_posthoc.sh` | Two variants of the injection, added after the first results of #56 to explain them: not part of version 1, part of version 2, where `run_experiments.sh` builds them |
| `run_injection.sh`, `evaluate_injection.py` | The injection part of the gate again, for [#69](https://github.com/naturalis/timegraft/issues/69): legacy, the default, and a control arm with `reanchor.check_ids: false`, each against its own clean build; metrics 1d, 1f and 1g, and why each other stem moves (`results/injection_69*.tsv`) |
| `monitor.py` | Runs a command and records its wall time and peak memory: of the largest single process over its whole life, from the kernel, and of the whole process tree, sampled every 0.2 s |
| `compare_backbones.py` | For [#84](https://github.com/naturalis/timegraft/issues/84): traces every stem age that differs by more than a factor between two backbones (format 5 `node_ages.tsv`) to its causes |
| `results_84/` | The 175 changed stems of #56 traced with `compare_backbones.py`, and how they were obtained |
| `results_70/` | The injection metrics (1a to 1m of version 2) after re-anchoring within the target family, [#70](https://github.com/naturalis/timegraft/issues/70), against #69 alone, and how they were obtained |
| `results_101/` | The injection metrics (1a to 1m of version 2) after finding a moved genus's species by name inside its target, [#101](https://github.com/naturalis/timegraft/issues/101), against #70, and the time and memory of the Iberian build with and without it |
| `results_79/` | The crown-tail metrics (3a to 3d) after bounding midpoint crowns by the barcodes, [#79](https://github.com/naturalis/timegraft/issues/79): the crowns that change, the genera exempt because their units lump species, and the check against TimeTree that rejected centring the midpoint on the cap |
| `results_108/` | [#108](https://github.com/naturalis/timegraft/issues/108): TimeTree's raw against its adjusted times on the backbone metrics (step 1), whether TimeTree's tree refines the NCBI taxonomy, and which node each answer dates, with the probe of unnamed answers and its validation (step 2) |
| `results_81/` | The injection, age perturbation and subsampling carried through the genus stage with `run_downstream.sh` ([#81](https://github.com/naturalis/timegraft/issues/81)): the downstream metrics of version 2, details per run and per injected genus, and why 6f fails (sampling noise; #104) |
| `evaluate.py`, `benchlib.py` | Compute every metric from the runs and write `results/`: `--version 3` (the default), `--version 2` or `--version 1` |
| `run_lineage_alone.sh`, `evaluate_lineage_alone.py` | The detection metrics again, for [#71](https://github.com/naturalis/timegraft/issues/71): the default (`concordance.lineage_alone_level: class`) and a control with it `off`, on the clean run and on both injection variants (merged ids and name only); flagged by rule, back in their own class, false flags (`results/lineage_alone*.tsv`) |

The pipeline is never patched for an experiment: errors and noise go into copies of its inputs (the taxonomy, the cache, the units), and the workflow scripts run as they are.

## Three versions of the gate

**Version 1** was declared for #56 and judged the layer as it was merged; the reruns after #89 (`results_89/`) and #73 (`results_73/`) used it too. Running it showed where its metrics were weak, so [#74](https://github.com/naturalis/timegraft/issues/74) declared **version 2** (section "Version 2" of `GATE.md`), committed on its own before any run that it judges. Version 2 is the gate for later changes to the layer. What it changes:

| | Version 1 | Version 2 |
| --- | --- | --- |
| 6b, changed stems explained | any fitted age on the path to the root counts as a cause (175 of 175 in #56; 159 on a stricter reading) | each changed stem followed to the nodes that moved it with `compare_backbones.py` (#84), on the two clean backbones; at least 0.95 traced |
| 7b, peak memory | at most 1.5 times legacy and below 7 GB | below 2,000 MB; the ratio is a report (7f) |
| Rooting | not measured | 8a genera with units on the wrong side of the root and 8b units among the outgroup, not more than legacy (within 5%); 8c crowns that depend strongly on the root, a report |
| Injection variants | `run_posthoc.sh`, not metrics | 1h to 1m: flagged, back in own class, stem within 1.5, for `inject_own` and `inject_name` (the stem share of `inject_name` is a report) |
| Downstream | the experiments stop at the backbone | 1n, 1o, 5e, 5f, 6e, 6f: crown ages and PD of the injection, age noise and subsampling carried through the genus stage ([#81](https://github.com/naturalis/timegraft/issues/81)) |
| Decision | must metrics of families 5 to 7 or metric 0 failing: hold | the same with family 8; must metrics not run are listed, and the decision says it was reached without them |

`evaluate.py --version 1` reproduces the tables of `results_73/` (values and statuses; the order of rows in `holdout.tsv` and of the phyla in the note of 1a vary between runs). Version 2 needs `node_ages.tsv` format 5 in the two clean backbones (built by any checkout from #84 on) and the two variants in `exp/`, and writes `stem_causes.tsv` and `decision.txt` besides the tables of version 1.

**Version 3** ([#104](https://github.com/naturalis/timegraft/issues/104)) changes only how the PD metrics of the downstream experiments, 1o, 6e and 6f, are measured. In version 2 each setting drew its own samples, so the two were judged on different samples, and 6f varied with the sample seed by as much as its threshold: it failed for 7 of 20 seeds although the settings did not differ on average (`results_81/pd_seed_variation.tsv`). Version 3 draws each replicate's samples once, from the units in all four trees compared, uses them for both settings, and takes 1,000 samples per design instead of 200 (`benchlib.paired_pd`). The noise of 6f falls to a standard deviation of 0.001 to 0.002 (`results_81/pd_paired_noise_1000.tsv`). The thresholds and every other metric are those of version 2; `GATE.md` gives the reason for changing a metric after a result, and the values already known on the runs of #81. `evaluate.py --version 2` gives what it gave before.

### The downstream runs

[#81](https://github.com/naturalis/timegraft/issues/81) carries the injection, the age perturbation and the subsampling through the genus stage, so that version 2 can measure what they do to crown ages and PD (1n, 1o, 5e, 5f, 6e, 6f). `run_downstream.sh` runs `downstream.py` on the clean backbones and on replicates 1 and 2 of the three experiments, for both settings: 14 runs, in `exp/downstream/` as `GATE.md` lays out (`clean/{legacy,default}/` and `{inject,noise,sub}/s{1,2}/{legacy,default}/`, each with `provenance.tsv`, `timegraft.tre` and `units/representatives.tsv`).

`downstream.py` runs `shard_genera.py`, `build_shard.py` (one process per shard), `graft.py` and `place_units.py` as the workflow does, on the units of the full run of the same setting (`exp/legacy_units` for legacy). For the subsampling it keeps only the units of the genera of `genera.txt` and the unassigned units whose genera are all kept, and then the records of the kept units. It keeps every `build_genus()` result in `exp/downstream/memo`, keyed on everything the function reads, and takes it from there when a later run calls it with the same inputs: the clean runs go first, and the experiments then rebuild only the genera whose inputs they changed. `build_genus()` is deterministic, so the cache changes the time and nothing else. A run is made in a `.partial` folder and moved into place when it has finished.

`legacy.yaml` gains `rate_outlier_action: flag` for the legacy arm: `main` at `e20fdf6` did not set long branches aside from their subtrees. The setting is read in the genus stage only, so the backbone and the screening of the legacy arm (metric 0) are unaffected.

Since the harness exists, `evaluate.py` (version 2) counts a downstream metric whose runs are missing as failed, as `GATE.md` requires; with a partly filled `exp/downstream/` it refuses to evaluate.

## Rerunning it

Needed: the OTT 3.7.3 taxonomy, the synthesis 16.1 tree and its `annotations.json` (as the workflow downloads them to `resources/otol`), MAFFT, VSEARCH and Snakemake (all in the `timegraft` conda environment: `conda activate timegraft`, see [installation](../../../../docs/technical/installation.md#with-conda-recommended)), and network access to TimeTree for the injected and subsampled builds (their picks are not in the cache). About three hours on two cores.

```bash
BENCH=/path/to/work                       # scratch folder, about 5 GB
git worktree add --detach $BENCH/main e20fdf6     # legacy: the pipeline before the consistency layer
git worktree add --detach $BENCH/branch main      # the code under test
for w in main branch; do (
  cd $BENCH/$w
  mkdir -p resources/iberia
  ln -s /path/to/otol resources/otol      # opentree16.1/{labelled_supertree.tre,annotations.json}, ott3.7.3/taxonomy.tsv
  gunzip -c data/iberia/records.tsv.gz > resources/iberia/records.tsv
  if [ $w = main ]; then                  # legacy: its own cache (format 1, #47)
    cp data/iberia/timetree_cache.json resources/api_cache.json
  else                                    # default: the Iberian cache together with the benchmark cache
    python data/iberia/consistency/benchmark/merge_caches.py resources/api_cache.json \
        data/iberia/timetree_cache.json data/iberia/consistency/benchmark/timetree_cache.json.gz
  fi
  # the full run, with wall time and peak memory (monitor.py is in this folder)
  python /path/to/this/folder/monitor.py $BENCH/${w}_run.json \
      snakemake --cores 2 --configfile config/config.yaml config/iberia.yaml
); done
(cd $BENCH/branch && python -m pytest -q tests | tail -1 > $BENCH/tests.txt)
./run_experiments.sh $BENCH $BENCH/branch $BENCH/branch/resources/otol 2    # add "offline" to forbid any API request
./run_posthoc.sh $BENCH $BENCH/branch $BENCH/branch/resources/otol 2    # version 1 only: run_experiments.sh builds them now
./run_downstream.sh $BENCH $BENCH/branch 2                                # version 2: the downstream runs (#81)
PYTHONPATH=$BENCH/branch/src python evaluate.py $BENCH results               # version 3; --version 2 for #74, --version 1 for #56
```

The two full runs must not overlap with each other or with anything else, or the time comparison (metric 7a) is void.

## Seeds and inputs

- Records: `data/iberia/records.tsv.gz`. Configuration: `config/config.yaml` with `config/iberia.yaml`.
- TimeTree and ChecklistBank answers: for the full runs, `data/iberia/timetree_cache.json` of `e20fdf6` for the legacy run, and that of the code under test merged with the benchmark cache for the default run; and `timetree_cache.json.gz` in this folder for the builds of `run_experiments.sh` ([#87](https://github.com/naturalis/timegraft/issues/87)). The benchmark cache holds every answer that the builds of the gate and the full default run asked for: 28,644 (15,913 TimeTree, 12,731 ChecklistBank), 1.01 MB gzipped. It started at 19,084; #101 added the 401 answers its injection builds asked for (`results_101/`), #108 the 4,078 NCBI lineages and probe queries of the node check (`results_108/`), and the gate run of #108 with `node_check.mode: report` the 5,081 answers of its builds (mostly probe queries for the re-query answers of the subsampled backbones, which take TimeTree about 7 s each: those builds took about an hour each online). It was merged by `merge_caches.py` from 101 caches: those of all builds of the reruns of #73 and #70, the clean builds of #84, and the full default run; not the perturbed copies of the noise runs, and not the run of `main` at `e20fdf6`, whose cache predates cache format 2 (#47) and which is only ever reused. No two caches answered a question differently. The default full run starts from the benchmark cache as well. Until #77, `data/iberia/timetree_cache.json` (then 9,820 answers) lacked the answers that #101 and #108 added; it now holds the 23,563 answers of the reference run, but the benchmark cache also holds those that only the experiments ask, so merging the two keeps the default run independent of what the cache of the next reference run happens to contain. Started from that cache alone, the default run of 25 November 2026 asked about 13,700 questions online: its backbone step took 2 h 22 min instead of minutes, which voids the time comparison (7a), and one live TimeTree answer differed from the one in the benchmark cache, so `run_experiments.sh` stopped when it merged the caches. `run_experiments.sh` starts every build from it together with the cache of the full run (`exp/start_cache.json`). With `offline` as its fifth argument, any request to TimeTree or ChecklistBank fails the build: on 17 November 2026 all 62 builds ran that way, with the network also blocked, and every backbone (tree, crown ages and every column of `node_ages.tsv`) and the default screening came out identical to the builds made online; `evaluate.py` gave the same value for every metric. A change to the pipeline that asks new questions needs the network once; add its answers afterwards with `python merge_caches.py timetree_cache.json.gz timetree_cache.json.gz $BENCH/exp/*/cache.json ...` (leave out `exp/noise/`, which the script refuses).
- Replicates 1 to 5 everywhere; the seed of a replicate is its number. Injection: `inject.py` (which genera, which hosts). Noise: `perturb_cache.py` (factor per query from the seed and the query key). Subsampling: `subsample.py`. Sample draws for PD: seeds 6 and 56 in `evaluate.py`. Hold-out: seed 52 in `../agefit/holdout.py`.

## Tables in `results/`

| Table | Contents |
| --- | --- |
| `metrics.tsv` | The gate |
| `injection_genera.tsv` | Every injected genus: where it was merged to, which checks could be made and failed, whether it was flagged and moved back, and its stem age against the run without the error, for both settings |
| `injection_variants.tsv`, `inject_own_genera.tsv`, `inject_name_genera.tsv` | The injection of the gate next to the two variants of `run_posthoc.sh`: a merge whose taxon holds only the genus's own species (the OTT entry of *Onchidella* holds those of both genera, as the injection of the gate does), and a name moved without NCBI ids or species |
| `injection_replicates.tsv` | Per replicate: genera flagged without an injected error, and other genera whose stem age moved, for both settings |
| `noise.tsv` | Per replicate and setting: how far node and stem ages move when the TimeTree ages are perturbed, and the rank correlation of backbone PD |
| `subsample.tsv` | Per replicate and setting: how far clade and stem ages move when half the genera are dropped |
| `holdout.tsv.gz` | Held-out TimeTree ages against the ages each method gives without them |
| `pd.tsv` | Faith's PD of the same random samples on the final trees of the two full runs |
| `stem_changes.tsv` | Genera whose stem age differs by more than 10% between the two full runs, with the causes the outputs record |
| `ablation.tsv` | The default backbone with one feature switched back at a time |
| `runtime.tsv` | Seconds per workflow rule, wall time and peak memory of the two full runs |
| `memory.tsv` | Peak memory of the backbone and screening steps on their own, measured with `monitor.py` around `backbone.py` and `screen.py`; made by hand, not by `evaluate.py` |
| `stem_causes.tsv` | Version 2: each stem that differs by more than 10% between the clean backbones, its causes and the chain behind them (`compare_backbones.py`) |
| `downstream.tsv` | Version 2: per downstream run, how far crown ages move and the PD rank correlation per design, against the clean run of the same setting |
| `decision.txt` | Version 2: the decision of the rule in `GATE.md`, and what it was reached without |
| `end_to_end.tsv` | The two full runs side by side: genera, units set aside, rooting, crown ages, placement, tree length |
