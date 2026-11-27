# Which part of the consistency layer gives its robustness (#80)

The gate shows that the consistency layer makes the backbone ages more robust than the pipeline before it (legacy):
- against noise in TimeTree's ages (metrics 5a, 5b);
- against dropping half the genera (5c, 5d).

Until now only the whole layer was compared with legacy. Here each part of the layer is switched back on its own, the same perturbation experiments are run for it, and each setting is compared with its own clean backbone, as the gate does.

## How it was made

On the gate run of #108 (`../results_108/gate/`, commit `bf58ca6`, whose pipeline is that of `main`), with its legacy and default builds as they are. Six ablations:

| Ablation | Setting | What it switches back |
|---|---|---|
| `precedence_ages` | `ages.method=precedence` | the weighted age fit of #52, back to the precedence rule |
| `mrca_anchoring` | `anchoring.rule=mrca`, `anchoring.prefer_sampled=false` | anchoring of #51 |
| `identifiers_from_ott` | `col.join=ott` | identifiers from COL XR (#48), back to OTT |
| `no_reanchoring` | `reanchor.mode=report` | re-anchoring of #50 (with #69, #70, #101) |
| `global_crowns` | `crown_query.basis=global` | the crown queries of #54 |
| `raw_timetree_time` | `ages.timetree_time.use=raw` | TimeTree's adjusted times of #108 |

The first five are the gate's ablation builds. The sixth came after #80 was written. For each ablation the script builds:

- its clean backbone without re-queries, the reference for noise (the gate's `default_norequery`);
- five noise replicates;
- five subsampled replicates against the gate's clean ablation build.

The seeds and the subsampled genus sets are those of the gate. The noise replicates run offline from caches perturbed as in the gate. The factor of a query depends only on the seed and the query's key, so the questions that legacy and the default ask have the factors they had in the gate.

```bash
./run_ablations.sh BENCH REPO OTOL 2      # BENCH: a gate folder after run_experiments.sh; 66 builds
python attribute.py BENCH .               # noise.tsv, subsample.tsv
```

The builds started from the benchmark cache together with the caches of the gate's ablation builds. The subsampled builds still asked TimeTree 476 new questions, and the longest of them took about five minutes. These answers, and those of the gate's ablation builds, are now in the benchmark cache, so a rerun asks nothing. `attribute.py` reproduces the gate's 5a to 5d exactly for legacy and the default.

## Results

Pooled over the five seeds (`noise.tsv`, `subsample.tsv`; one row per seed there too). The factor is that by which an age moves against the arm's own clean backbone.

**TimeTree's ages times a random factor (sigma 0.3)**

| Arm | 90% quantile, node ages (5a) | Share moving by more than 1.5 (5b) | 90% quantile, genus stems |
|---|---|---|---|
| Legacy | 1.769 | 0.210 | 1.750 |
| Default | 1.492 | 0.092 | 1.492 |
| `precedence_ages` | **1.775** | **0.218** | **1.732** |
| `mrca_anchoring` | 1.492 | 0.094 | 1.492 |
| `identifiers_from_ott` | 1.482 | 0.091 | 1.483 |
| `no_reanchoring` | 1.492 | 0.095 | 1.493 |
| `global_crowns` | 1.492 | 0.092 | 1.492 |
| `raw_timetree_time` | 1.487 | 0.093 | 1.489 |

**Half the genera dropped**

| Arm | 90% quantile, clade ages (5c) | Share moving by more than 1.5 | 90% quantile, genus stems (5d) |
|---|---|---|---|
| Legacy | 1.429 | 0.092 | 1.866 |
| Default | 1.255 | 0.056 | 1.495 |
| `precedence_ages` | **1.329** | 0.073 | **1.651** |
| `mrca_anchoring` | 1.279 | 0.066 | **1.649** |
| `identifiers_from_ott` | 1.238 | 0.055 | 1.494 |
| `no_reanchoring` | 1.270 | 0.057 | 1.506 |
| `global_crowns` | 1.255 | 0.056 | 1.495 |
| `raw_timetree_time` | 1.271 | 0.062 | **1.604** |

## What it means

**Robustness to noisy ages comes from the weighted fit alone.** Without it the backbone is exactly as sensitive as legacy (5a 1.775 against 1.769, 5b 0.218 against 0.210). With it, the worst of the five seeds (1.505) is better than the best seed without it (1.642). Every other ablation stays within 0.01 of the default. This is what the fit is for: precedence lets one answer replace another, so a noisy answer moves its node by its full error. The fit weighs every calibration of a node and its neighbours together, so one noisy answer moves a node by only part of its error.

**Robustness to sampling has several sources, and the fit is the largest.** The default gains 0.174 on 5c over legacy.

- Without the fit, 0.074 of that is lost, about 40%, and the loss shows in every seed (0.04 to 0.19).
- Anchoring (0.024), TimeTree's adjusted times (0.016) and re-anchoring (0.015) each account for less.
- The rest of the gap is not attributed to one part. Legacy also differs in its units, from the screening of #55, and the parts may not add up when switched back one at a time.

For stems (5d), the fit and anchoring each account for about 0.15 of the 0.37 gained, and the adjusted times for about 0.11. But the stem quantile varies between seeds as much as between some arms (1.36 to 1.73 for the default), so these shares are rough.

**Two ablations change nothing, as expected, or slightly the other way.**

- `global_crowns` gives identical backbones, because crown queries date genus crowns, not backbone nodes.
- `identifiers_from_ott` is as robust as the default or slightly more (5c 1.238), so the COL identifiers of #48 are not a source of robustness. They earn their place for what they anchor (`../../backbone-col/`), not for this.

**The fit is worth its complexity on this evidence.** On held-out TimeTree ages it is no better than precedence (gate metric 4a, 1.309 against 1.308). But it accounts for all of the robustness to noisy ages, and for the largest share of the robustness to sampling. So the issue that #80 asked for "if the fit does not account for the gain" is not needed.
