# TimeTree's times and the nodes it dated (#108)

Step 1 chooses which of TimeTree's two times to calibrate with; step 2 checks which node each answer dates.

# Step 1: which of TimeTree's two times to calibrate with

Every TimeTree answer has two times:
- **`age`:** the summary over studies, a robust average (the Huber estimator, by TimeTree's FAQ).
- **`adjusted_age`:** the time after TimeTree smooths its tree to make it ultrametric, where conflicting studies left an ancestor younger than a descendant.

The backbone was calibrated with `age` only. Step 1 of [#108](https://github.com/naturalis/timegraft/issues/108) records both and compares four ways of using them (`ages.timetree_time.use`):

| Variant | Calibration |
| --- | --- |
| `raw` | `age` (before #108) |
| `downweight` | `age`, with its sigma multiplied by 3 where `age` and `adjusted_age` differ by more than a factor 1.5 |
| `replace` | `age`, but `adjusted_age` where the two differ by more than a factor 1.5 |
| `adjusted` | `adjusted_age` throughout |

## Method

Every build is offline from the caches of the benchmark of #73, on the Iberian units of #73 and the full OTT 3.7.3 and synthesis 16.1, with `ages.requery.mode: off` (re-querying changes no age):
- **Clean:** one build per variant, from `exp/default_clean/cache.json` of #73.
- **Noise:** five builds per variant, from `exp/noise/s1` to `s5/cache.json` of #73. These are the union cache with every TimeTree answer multiplied by a random factor (sigma 0.3); `perturb_cache.py` multiplies the adjusted time by the same factor as the summary.

The metrics are those of the backbone part of the gate, computed as `evaluate.py` does but between the variants rather than against legacy:
- 2a and 2b: `metadata.json`, `ages.final`.
- 2c: nodes older than their parent.
- 4a and 4b: `agefit/holdout.py` with 50 repeats.
- 5a and 5b: node ages of each noise build against the clean build of the same variant.

They also include how far each variant moves the backbone from `raw` (node ages and stems by more than 10% or 1.5), and the stem of *Eusergestes* (#78). `compare_variants.py` computes the table from the builds; its paths are those of the machine it ran on. `raw` gives a backbone byte-identical to `main`.

## Result

`variants.tsv`:

| | `raw` | `downweight` | `replace` | `adjusted` |
| --- | --- | --- | --- | --- |
| Calibrations | 387 | 387 | 387 | 330 |
| 2a calibrations off by more than 1.5 | 5 | 5 | 4 | 2 |
| 2b weighted RMS of residuals | 0.384 | 0.352 | 0.354 | 0.218 |
| 2c nodes older than their parent | 0 | 0 | 0 | 0 |
| 4a hold-out, median error factor | 1.318 | 1.319 | 1.309 | 1.272 |
| 4b hold-out, 90% quantile | 2.485 | 2.485 | 2.324 | 2.271 |
| 5a noise, 90% quantile of node moves | 1.487 | 1.483 | 1.492 | 1.535 |
| 5b noise, share of nodes moving more than 1.5 | 0.093 | 0.091 | 0.092 | 0.115 |
| Noise, 90% quantile of stem moves | 1.489 | 1.484 | 1.492 | 1.561 |
| Nodes moved from `raw` by more than 10% / 1.5 | 0 / 0 | 18 / 0 | 53 / 22 | 220 / 40 |
| Stems moved from `raw` by more than 10% | 0 | 27 | 88 | 356 |
| Stem of *Eusergestes* and *Deosergestes* (Ma) | 7.7 | 7.7 | 44.6 | 44.6 |

**Metrics 2a, 2b, 4a and 4b cannot decide this alone.** Adjusted times are consistent by construction, so they fit and predict each other better, and these metrics favour any variant that uses more of them. Robustness (5a, 5b) compares like with like, because each variant is measured against itself.

**`downweight` does nothing where it is aimed.** The answers whose two times differ most conflict with no neighbour in our backbone, so the fit keeps them at their age whatever their weight. The shrimp stem stays at 7.7 Ma.

**`adjusted` changes the backbone broadly.** It leaves 57 calibrations out, all `interpolated_unresolved_in_timetree`. Where TimeTree's raw time for a child is older than its parent's, its adjusted time is the parent's: TimeTree pools the two into one time, as an isotonic fit does. Our backbone then treats the child as not resolved by TimeTree and interpolates it. It is also the least robust variant: 5a goes from 1.487 to 1.535 and 5b from 9.3% to 11.5%. That is still well within the gate's thresholds against legacy, which were 1.77 and 21% in #56.

**`replace` takes the adjusted time only where the two differ markedly.** It changes 17 calibrations (`replaced.tsv`); 11 of them become older and 6 younger, 8 rest on a single study, and 9 are unnamed TimeTree nodes. Robustness is unchanged (1.492 and 9.2%). It moves 53 node ages and 88 stems by more than 10%. Among the 17:

| Genera | Raw (Ma) | Adjusted | Age before | Age after |
| --- | --- | --- | --- | --- |
| *Deosergestes*, *Eusergestes*, *Sergia* | 7.7 | 44.6 | 7.7 | 44.6 |
| *Ampelisca*, *Apolochus*, *Dexamine* | 22.9 | 102.0 | 22.9 | 102.0 |
| *Caryophyllia*, *Dendrophyllia*, *Oculina* | 49.7 | 188.2 | 49.7 | 188.2 |
| *Anodonta*, *Potomida*, *Unio* | 61.3 | 134.5 | 61.3 | 134.5 |
| *Aurelia*, *Cassiopea*, *Pelagia* | 205.2 | 385.2 | 205.2 | 385.2 |
| *Ameiurus*, *Ictalurus*, *Pangasianodon* | 119.9 | 64.9 | 119.9 | 64.9 |

The catfish node shows the change can go either way: there TimeTree's adjusted time is younger than its summary.

**The factor.** `factors.tsv` gives `replace` with factors 1.25, 1.5 and 2.0:

| Factor | Nodes moved by more than 10% | Stems | 2a | 5a | 5b |
| --- | --- | --- | --- | --- | --- |
| 1.25 | 100 | 164 | 1 | 1.492 | 0.092 |
| 1.5 | 53 | 88 | 4 | 1.492 | 0.092 |
| 2.0 | 9 | 27 | 5 | 1.487 | 0.093 |

Robustness is the same for all three, and all three fix the shrimp stem. The factor sets how far TimeTree's smoothing is allowed to reach into the backbone.

## Decision

The default is `replace` with factor 1.5:
- **It fixes the clear cases.** These are the answers that TimeTree's own tree contradicts by more than a factor 1.5, the threshold the gate uses for a severe conflict.
- **It costs no robustness.**
- **It leaves the 250 smaller adjustments alone.** Those are where `adjusted` merges nodes and loses robustness.

`raw` remains available and gives the ages before #108.

## What this does not show

- **Whether the adjusted times are right.** The gate has no external truth for backbone ages. The adjusted times are TimeTree's own consistent estimate, and several are more plausible than the raw ones, for example a stony-coral split of 188 rather than 50 Ma. But that is TimeTree's judgement, not an independent check.
- **Whether TimeTree dated the intended node.** That is step 2 of #108, and it may change some of these calibrations again.
- **The downstream metrics.** 88 stems change, so genus crowns and PD move. They are measured in the full gate run after step 2.

# Which node TimeTree dated (#108, step 2)

## Step 2.1: TimeTree's tree refines the NCBI taxonomy

`ncbi_refinement.py` tests this on every named TimeTree answer:
- the 165 answers calibrating the backbone (`refinement_backbone.tsv`);
- the 2,195 in the benchmark cache (`refinement_cache.tsv`).

If TimeTree's tree only resolves polytomies of the NCBI taxonomy, a named node must be the NCBI common ancestor of the ids it recognised. NCBI lineages come from ChecklistBank's copy of the NCBI taxonomy (dataset 2169).

| | Same node | TimeTree's node deeper | Lineage unknown |
| --- | --- | --- | --- |
| Backbone | 163 | 0 | 2 |
| Whole cache | 2,158 | 6 | 31 |

The six deeper answers are NCBI release differences, not departures from NCBI:
- **Five name Pancrustacea** where ChecklistBank's NCBI has *Allotriocarida* or *Altocrustacea* (ids above 3,700,000, added recently) in between.
- **One names the genus *Halichondria*** where NCBI now has a subgenus of the same name in between.

So TimeTree's tree refines a somewhat older NCBI taxonomy. For every answer, named or not, the node TimeTree dated therefore lies within the NCBI common ancestor of the ids it recognised.

## Step 2.2: the check

`timegraft.ttcheck`, run by the backbone build with `ages.node_check`, compares the genera of the run within that common ancestor with those of the backbone node.

**Unnamed answers that come out `wider`.** These are clades TimeTree made when resolving an NCBI polytomy, and they may still match the node. So they are probed: TimeTree is asked again with one outside genus added, nearest first, up to three genera it knows. If TimeTree's node id (`mrca_ttid`) stays the same, the genus is inside its node.

**The probe was validated on the named `wider` answers** (`probe_validation_named.tsv`, made with `probe_check.py ... named`). For those, NCBI already says the outside genera are inside TimeTree's node. All 31 the probe could decide came back inside, none falsely outside; 10 could not be decided.

`node_check.tsv` gives the outcome per answer on the Iberian backbone (step 1's default `replace`, with `node_check.mode: reject`):

| Outcome | Named | Unnamed | Used as calibration |
| --- | --- | --- | --- |
| `match` | 118 | 112 | yes |
| `probed_outside` | | 97 | yes |
| `unchecked` | 1 | 16 | yes |
| `wider` | 41 | 62 | no |
| `narrower` | 5 | 2 | no |

110 of the 454 answers are rejected, and their nodes are interpolated. Before this check, the named comparison of the issue found 83 of 117 answers dating the intended node. With NCBI lineages instead of OTT's, 118 of 164 named answers match. After the check, every named calibration matches by construction. (`node_check_prototype.tsv` and `probe_check_prototype.tsv` are the same check made outside the pipeline with `node_check.py` and `probe_check.py`; they agree with it up to the choice of probes.)

`ages.node_check.mode: reject` is the default.

The other part of step 2 in the issue, querying with more ids per side, was not done. With the check, whether an answer is used no longer depends on which few ids TimeTree recognises, only on which node it dated. Re-asking every node with new picks would change all 454 answers at once and hide the effect of the check.

## Result

`step2.tsv`, measured as in step 1 (`check` is `replace` with the node check):

| | `raw` (before #108) | `replace` (step 1) | `replace` + check |
| --- | --- | --- | --- |
| Calibrations | 387 | 387 | 335 |
| 2a calibrations off by more than 1.5 | 5 | 4 | 3 |
| 2b weighted RMS of residuals | 0.384 | 0.354 | 0.238 |
| 2c nodes older than their parent | 0 | 0 | 0 |
| 4a / 4b hold-out | 1.318 / 2.485 | 1.309 / 2.324 | 1.306 / 2.371 |
| 5a / 5b noise | 1.487 / 0.093 | 1.492 / 0.092 | 1.488 / 0.091 |
| Nodes moved from `raw` by more than 10% / 1.5 | | 53 / 22 | 128 / 50 |
| Stems moved from `raw` by more than 10% | | 88 | 179 |

The noise builds of the check use the perturbed caches of #73, with the probe queries and lineages of the clean build added. Only TimeTree's node ids matter for the check, so the same answers are rejected with and without noise.

`stem_causes.tsv` (`compare_backbones.py` between `raw` and the check) traces all 179 stems that move by more than 10%:
- 128 to `node_check`;
- 35 to `timetree_time` (step 1);
- 20 to the fit around them.

## Test cases

| Node | Genera | TimeTree (Ma) | Check | Before (Ma) | After |
| --- | --- | --- | --- | --- | --- |
| `n1373` | *Galathowenia*, *Magelona*, *Owenia* | 394.5 | `wider` (41 genera outside) | 394.5 | 134.3 |
| `n1811` | *Deosergestes*, *Eusergestes*, *Sergia* (#78) | 7.7 (adjusted 44.6) | `probed_outside` | 7.7 | 44.6 |
| `n490` | *Somniosus*, *Squalus* | 138.6 | `wider` | 107.0 | 78.7 |
| `n360` to `n363` | grenadiers (11, 9, 8, 7 genera) | 48.1, 48.1, 64.8, 48.1 | `wider`, `wider`, `narrower`, `wider` | 54.3, 53.9, 53.8, 48.1 | 48.2, 48.0, 48.0, 47.6 |

**The grenadiers** no longer take TimeTree's answer for one 20-genus taxon. They are interpolated now, but their neighbours bound them closely (46.3 Ma below, about 48 above), so their ages stay close together.

**The *Lepadogaster* and barnacle nodes of the issue** do not change. Their TimeTree age equals their parent's, so they were already interpolated before #108; they are now recorded as `wider` as well.

## What this does not show

- **Whether the interpolated ages are better than the rejected answers.** The gate has no external truth for backbone ages. Each rejected answer is TimeTree's age for another node, so it cannot be the age of the node asked about. For example, the annelid node got the age of a taxon holding 44 genera of the run.
- **The 17 unchecked answers.** These are still used.
- **The downstream metrics and the rest of the gate.** They are left to the full gate run with both steps merged.
