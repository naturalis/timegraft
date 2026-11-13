# Benchmark and release gate of the consistency layer (#56)

Whether `consistency-layer` should become the default: the metrics and thresholds, the harness that measures them on the Iberian run, and the outcome.

| File | Contents |
| --- | --- |
| `GATE.md` | The metrics, their thresholds and the decision rule, committed before anything was run |
| `REPORT.md` | The outcome: status of every metric, what the numbers mean, and the recommendation |
| `results/metrics.tsv` | The gate as a table: metric, legacy, default, threshold, status |
| `results/*.tsv` | One table per experiment (see below) |
| `legacy.yaml` | Settings that switch the layer off: `consistency-layer` run as `main` ran |
| `backbone.py`, `screen.py` | Run the unchanged workflow scripts for the backbone and the screening outside Snakemake, with settings from the configuration files plus overlays |
| `inject.py` | Merges genera of the run into taxa elsewhere in a copy of the OTT taxonomy, as OTT did with *Onchidella* |
| `perturb_cache.py` | Multiplies every TimeTree age in a copy of the cache by a random factor |
| `subsample.py` | Keeps a random half of the genera |
| `run_experiments.sh` | All backbone builds of the gate: clean, injected, perturbed, subsampled, ablations, reproducibility |
| `run_posthoc.sh` | Two variants of the injection, added after the first results to explain them; not part of the gate |
| `monitor.py` | Runs a command and records its wall time and the peak memory of its process tree |
| `evaluate.py`, `benchlib.py` | Compute every metric from the runs and write `results/` |

The pipeline is never patched for an experiment: errors and noise go into copies of its inputs (the taxonomy, the cache, the units), and the workflow scripts run as they are.

## Rerunning it

Needed: the OTT 3.7.3 taxonomy, the synthesis 16.1 tree and its `annotations.json` (as the workflow downloads them to `resources/otol`), MAFFT, VSEARCH, Snakemake, and network access to TimeTree for the injected and subsampled builds (their picks are not in the cache). About three hours on two cores.

```bash
BENCH=/path/to/work                       # scratch folder, about 5 GB
git worktree add $BENCH/main main         # e20fdf6
git worktree add $BENCH/branch consistency-layer
for w in main branch; do (
  cd $BENCH/$w
  mkdir -p resources/iberia
  ln -s /path/to/otol resources/otol      # opentree16.1/{labelled_supertree.tre,annotations.json}, ott3.7.3/taxonomy.tsv
  gunzip -c data/iberia/records.tsv.gz > resources/iberia/records.tsv
  cp data/iberia/timetree_cache.json resources/api_cache.json
  # the full run, with wall time and peak memory (monitor.py is in this folder)
  python /path/to/this/folder/monitor.py $BENCH/${w}_run.json \
      snakemake --cores 2 --configfile config/config.yaml config/iberia.yaml
); done
(cd $BENCH/branch && python -m pytest -q tests | tail -1 > $BENCH/tests.txt)
./run_experiments.sh $BENCH $BENCH/branch $BENCH/branch/resources/otol 2
./run_posthoc.sh $BENCH $BENCH/branch $BENCH/branch/resources/otol 2    # optional: the two variants
PYTHONPATH=$BENCH/branch/src python evaluate.py $BENCH results
```

The two full runs must not overlap with each other or with anything else, or the time comparison (metric 7a) is void.

## Seeds and inputs

- Records: `data/iberia/records.tsv.gz`. Configuration: `config/config.yaml` with `config/iberia.yaml`.
- TimeTree and ChecklistBank answers: `data/iberia/timetree_cache.json` of each branch. The cache after the benchmark, with the queries of the injected and subsampled builds, is not kept: those builds ask TimeTree again.
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
| `end_to_end.tsv` | The two full runs side by side: genera, units set aside, rooting, crown ages, placement, tree length |
