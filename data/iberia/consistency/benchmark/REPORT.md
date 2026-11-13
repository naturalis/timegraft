# Release gate for the consistency layer: outcome

Run on 13 November 2026 for [#56](https://github.com/naturalis/timegraft/issues/56), on the Iberian data. `consistency-layer` at `b2e30a8` against `main` at `e20fdf6`. The metrics, thresholds and decision rule are in `GATE.md`, committed before the run; how to rerun it is in `README.md`; the numbers are in `results/`.

## Verdict

**25 of 28 must metrics pass. Three fail: 1f, 1g and 7b.**

By the rule declared in advance, that is a **hold**, because 7b (peak memory) is an operations metric. Without 7b the rule gives **ship as an option**, because 1f and 1g are known-error metrics.

Reading the three failures for what they show, the recommendation is different from the rule's, and that is for the maintainer to accept or not:

- **7b fails on a threshold that was set badly, not on a problem.** Peak memory rose from 0.76 to 1.28 GB, a factor 1.68 against an allowed 1.5. The absolute limit in the same metric was 7 GB. All of the excess over the allowed 1.14 GB is the backbone step reading `annotations.json` whole (0.23 GB), which is easy to avoid.
- **1f and 1g are real gaps, found by the benchmark, in how completely an injected error is undone.** They do not show the default doing worse than `main` on any real result: on the same injected errors `main` leaves 104 of 120 genera in the wrong class, the default 2. But two things the correction does not yet do are now measured, and both have a clear fix (below).
- **Everything else improved or stayed the same**: the known real error is fixed end to end, conflicts and implausible crowns fall, the backbone is less sensitive to noise in TimeTree and to which genera are sampled, sample rankings are stable, and the run is reproducible offline.

**Recommendation: merge `consistency-layer` into `main` with its defaults, and open the follow-ups listed under "What to do".** None of the gaps changes the Iberian result itself: the only flagged genus there is *Onchidella*, which is corrected, and none of its TimeTree queries picked up a species of the other animal (see below; that is luck, not design). Whoever prefers to honour the rule to the letter should do follow-ups 1 to 4 first and rerun the gate (three hours).

## The gate

Legacy is `main`; default is `consistency-layer`. "Known" and "blind" say whether the result could have been known when the threshold was set (`GATE.md`).

| # | Metric | Legacy | Default | Threshold | Status |
| --- | --- | --- | --- | --- | --- |
| 0 | `consistency-layer` with the layer switched off reproduces `main` | | `backbone.tre` byte-identical; 0 of 24,882 screening decisions and 0 of 312 TimeTree crowns differ | identical | pass (blind) |
| 1a | *Onchidella*'s nearest genus in the final tree | an echinoderm | a mollusc (in the backbone: among 325 molluscs) | all Mollusca | pass |
| 1b | Injected misplacements with two checks possible: flagged | no check exists | 118 of 118 | at least 0.90 | pass (blind) |
| 1c | Injected misplacements, all: flagged | no check exists | 118 of 120 | at least 0.50 | pass (blind) |
| 1d | Flagged genera placed back in their own class | 16 of 120 in their own class | 117 of 118 | at least 0.90 | pass (blind) |
| 1e | Genera without an injected error flagged, per replicate | | 0, 0, 0, 0, 0 | at most 1 | pass (blind) |
| 1f | Corrected genera whose stem age is within a factor 1.5 of the one without the error | 21 of 120 | 85 of 117 (0.73) | at least 0.80 | **fail** (blind) |
| 1g | Other genera whose stem age moves by more than 10% because of the errors | 248 | 297 | at most half of legacy | **fail** (blind) |
| 2a | TimeTree calibrations off from the final age by more than a factor 1.5 | 26 of 381 | 5 of 387 | at most half | pass (known) |
| 2b | RMS of calibration residuals over their uncertainty | 1.39 | 0.38 | lower | pass (known) |
| 2c | Nodes older than their parent | 0 | 0, in the backbone and the final tree | 0 | pass (blind) |
| 3a | Crown age more than 20 times the barcode age | 7 | 0 | 0 | pass |
| 3b | Crown age more than 10 times the barcode age | 40 | 25 | fewer | pass |
| 3c | TimeTree crown used that is less than a third of the barcode age | 5 | 0 | 0 | pass |
| 3d | Root-to-tip distance below 0.02 with a crown above 20 Ma | 4 | 1 | none | report |
| 4a | Held-out TimeTree ages: median error factor | 1.319 | 1.318 | not worse | pass (known) |
| 4b | Held-out TimeTree ages: error factor at the 90% quantile | 2.50 | 2.49 | not worse | pass (known) |
| 5a | TimeTree ages perturbed: factor by which node ages move, 90% quantile | 1.77 | 1.49 | not worse | pass (blind) |
| 5b | TimeTree ages perturbed: node ages that move by more than a factor 1.5 | 21.0% | 9.3% | not worse | pass (blind) |
| 5c | Half the genera dropped: factor by which clade ages move, 90% quantile | 1.43 | 1.27 | not worse | pass (blind) |
| 5d | Half the genera dropped: the same for genus stems | 1.87 | 1.60 | not worse | pass (blind) |
| 6a | PD of the same samples on the two final trees: rank correlation | | 0.989 (150 units), 0.990 (100 genera, one unit), 0.972 (100 genera, five units) | at least 0.90 each | pass (blind) |
| 6b | Stem ages that differ by more than 10%: share with a recorded cause | | 175 of 175 as declared; 159 of 175 on a stricter reading | at least 0.90 | pass (blind) |
| 6c | Backbone PD of genus samples: rank correlation with and without perturbed ages | 0.912 | 0.942 | not lower | pass (blind) |
| 6d | Total branch length of the final tree; PD ratio default over legacy | 638,002 Ma | 630,769 Ma (-1.1%); PD 1.011, 1.013, 0.996 | none | report |
| 7a | Wall time of the full run, two cores | 48.0 min | 52.6 min (1.09 times) | at most 1.5 times | pass (blind) |
| 7b | Peak memory of the full run | 758 MB | 1,277 MB (1.68 times) | at most 1.5 times, and below 7 GB | **fail** (blind) |
| 7c | Backbone and screening rebuilt | | byte-identical: backbone twice, harness against workflow, screening twice | identical | pass (blind) |
| 7d | Backbone rebuilt with the network blocked | | built; all files identical | identical | pass (blind) |
| 7e | Test suite | | 235 passed, 2 skipped (no LSD2, no R) | passes | pass |

## What the numbers mean

### The comparison is sound (metric 0)

With `legacy.yaml`, `consistency-layer` gives the backbone of `main` byte for byte, the same screening decision for every unit and the same TimeTree crowns. So every difference below comes from the layer's settings and from nothing else on the branch, and each part of the layer can still be switched off.

### Known errors (1a to 1g)

**The real case is fixed, end to end.** In the final tree of `main`, the nearest genus to the gastropod *Onchidella* is a sea urchin; in the default it is a mollusc.

**Injected merges are found and moved back.** 120 genera in five replicates were merged in a copy of OTT into a taxon of another phylum or class, which is what OTT did to *Onchidella*. The default flags 118 and places 117 of those back in their own class, with no genus flagged that had no error. The two it misses have no NCBI id of their own, so the merged-id check has nothing to compare; both are listed for review. `main` has no check: 16 of the 120 end up in their own class by accident of where the host sits.

**Two things the correction does not yet do** (1f and 1g). Of the 117 corrected genera, 32 have a stem age more than a factor 1.5 from the one they had without the error (`main`: 99 of 120), and the errors move the stems of more *other* genera in the default (297) than in `main` (248). `results/injection_genera.tsv` and a rerun with a variant of the injection (`results/injection_variants.tsv`, added afterwards and not part of the gate) separate two causes:

1. **The species of the wrong taxon still go to TimeTree.** When OTT merges two genera, the merged taxon holds the species of both. After re-anchoring, the genus is in the right place, but its queries to TimeTree still include species of the other animal, and the node above it gets an age of several hundred Ma (a clade of three genera at 721 Ma, in replicate 1). The fit down-weights some of these as outliers, but not all, and the ages around the genus rise. In the variant where the merged taxon holds only the genus's own species, 96 of 117 stems are within 1.5 (0.82) and 142 other stems move instead of 297, against 251 for `main`.
2. **The genus goes back to its family, not to its place within the family.** Re-anchoring targets the family of the reference lineage. Where the genus used to sit inside a resolved family, its stem becomes the family's crown, and its former sister loses its sister. This accounts for all 21 remaining misses in the variant, and for most of the 142 other stems.

Neither shows in the Iberian run, but the first could have. OTT's entry for *Onchidella* holds 32 names of the gastropod and 36 of the sand dollar *Peronella*, two of them with an NCBI id. No TimeTree query for a node of this run contains one of those ids, and the genus has no TimeTree crown. With other genera in the run, and so other picks, one could. The second does not show because no other Iberian genus is in its family, so the family is its place.

**The limit of detection** (also from a variant, not part of the gate). If only the name moves to the wrong taxon, with no second NCBI id to give the merge away, the default flags 29 of 120 and lists the other 91 for review: none goes unnoticed, but three quarters are not corrected. The reason is the rule that two checks must fail. The lineage check fails for all 120. In the unperturbed run it fails for one genus of 1,218 (*Onchidella*); the 68 other lineage differences are at family or order level and stay below the `flag_level`. The barcode-neighbour check compares phyla, so it cannot fail for a move within a phylum, and it needs five hits. On this evidence the lineage check alone, at class level or above, would have caught all 240 injected errors of both kinds with no false flag.

### Backbone conflicts (2a to 2c) and held-out ages (4a, 4b)

The fit leaves 5 calibrations more than a factor 1.5 from their final age, where the precedence rule left 26, and the largest factor falls from 6.5 to 3.1. That was known from #52 and is what the fit optimises. The independent test is the hold-out: where a node has no calibration of its own, the two methods predict TimeTree's age equally well (1.319 and 1.318). So the fit honours more of the evidence without predicting unseen ages any better. Its measurable gain is in robustness, next.

### Robustness (5a to 5d)

This is the strongest new result, and it was blind.

- **Noise in TimeTree.** With every TimeTree age multiplied by a random factor (sigma 0.3), 21% of backbone node ages move by more than a factor 1.5 under `main` and 9% under the default; the 90% quantile of the move falls from 1.77 to 1.49. The default is better in each of the five replicates.
- **Which genera are sampled.** With half the genera dropped, the ages of the remaining clades move less under the default (90% quantile 1.27 against 1.43; better in five of five replicates), and so do genus stems (1.60 against 1.87; four of five).

Which part of the layer produces this was not tested separately. The fit is the likely source for the noise result, since it weighs all calibrations together where the precedence rule lets one age remove another.

### Crown ages (3a to 3d)

End to end, no genus has a crown age more than 20 times what its barcodes give (7 before), 25 are above 10 times (40 before), and no TimeTree crown is used that is less than a third of the barcode age (5 before). 15 TimeTree crowns are demoted to an upper bound. The sum of crown ages falls by 3.4%.

### Diversity (6a to 6d)

The layer corrects parts of the tree; it does not reshuffle it. Faith's PD of the same random samples ranks them almost identically on the two final trees (rank correlation 0.97 to 0.99), and PD itself changes by about 1%. Under noise in TimeTree, the ranking of genus samples on the backbone is more stable with the default (0.94 against 0.91; four of five replicates).

175 genus stems differ by more than 10% between the two runs, 136 of them older, 76 by more than a factor 1.5. The largest single change is the stem of the shrimps *Eusergestes* and *Deosergestes*, from 180 to 8 Ma: with COL's identifier for *Eusergestes*, TimeTree now dates their split, where the age was interpolated before. Whether an age that young is right is worth a look by someone who knows the group. As declared, metric 6b counts a stem as explained when any node on its path to the root has a fitted age that departs from the precedence age, and that holds for all 175, which makes the metric too easy. On a stricter reading (a new anchor, another identifier, or a fitted age at the nearest dated ancestor) 159 of 175 are explained. The other 16 are not explained by this reading. Some share their stem with a genus that is, such as *Deosergestes* with *Eusergestes*; the outputs record the TimeTree queries of every node, but the metric does not follow them.

### Operations (7a to 7e)

- **Time**: 52.6 against 48.0 minutes. Screening takes 297 s instead of 102 s (the second search of #55); the backbone 25 s instead of 17 s. The genus subtrees (51 to 53 minutes of CPU) and the placement of units (16 minutes) dominate both runs and did not change.
- **Memory**: the peak of both runs is the backbone step. On `consistency-layer` it takes 908 MB with the layer off and 1,143 MB with it on; without reading `annotations.json` it takes 910 MB (`results/memory.tsv`).
- **Reproducibility**: the default backbone is byte-identical when built twice, when built by the harness and by the workflow, and when built with the network blocked. The screening is byte-identical when run twice.

## What each part of the layer contributes

The default backbone with one part switched back at a time (`results/ablation.tsv`), and what the two full runs show.

| Part | Effect on the Iberian run | Verdict |
| --- | --- | --- |
| Evidence and provenance (#47) | No change to any result. Makes every difference in this report traceable | Keep |
| COL as identifier hub (#48) | Three genera join the backbone (1,183 to 1,186). Switched back: 10 stems move by more than 10% | Keep |
| Placement concordance and re-anchoring (#49, #50) | One genus flagged and moved, correctly. On injected merges: 118 of 120 found, none wrongly. Two gaps (follow-ups 1 and 2) and a detection limit (follow-up 3) | Keep; improve |
| Majority-clade anchoring (#51) | 62 genera anchored on their majority clade. Switched back: 23 stems move by more than 10%, 14 by more than a factor 1.5 | Keep |
| Weighted fit of ages (#52) | The largest effect: switched back, 136 stems move by more than 10% and 58 by more than 1.5; calibrations off by more than 1.5 go from 5 to 26. No better on held-out ages; most likely the source of the gain in robustness | Keep |
| Re-query of conflicting nodes (#53) | Changes no age at the default (`classify`): the default without it is identical. A diagnostic | Keep as it is; do not switch to `apply` on this evidence |
| Sampled crowns and plausibility (#54) | 83 TimeTree crowns differ from the global ones; the extreme crowns are gone | Keep |
| Screening (#55) | 1,139 units set aside instead of 1,294, because congeners are no longer crowded out of the hits. 265 units flagged as likely nuclear copies and 432 as long branches, not acted on | Keep; see the note on rooting |

**A side effect to know about.** The corrected screening keeps 155 more units in their genus, and some of them sit on the wrong side of the root: 386 genera have such units instead of 366, and 314 units fall among the outgroup instead of 276. The number of genera whose crown age depends strongly on the root is unchanged (19 against 20). This was not a gate metric.

The two full runs side by side are in `results/end_to_end.tsv`.

## What to do

1. **After re-anchoring, send only the genus's own species to TimeTree.** For a flagged genus, keep the species whose NCBI lineage agrees with the reference lineage (the lookup that the merged-id check already makes), and drop the others from the node queries and the global crown. Expected from the variant: stems within 1.5 from 0.73 to about 0.82, and other stems moved from 297 to about 142. *Onchidella* keeps its gastropod species and loses the sand dollars.
2. **Re-anchor within the family where the genus's species say where.** After choosing the target taxon, anchor on the majority clade of the genus's species inside it (the rule of #51), and fall back to the taxon when none are there. This addresses the remaining 21 misses and most of the remaining collateral. For *Onchidella* nothing changes: the synthesis has its species under the merged genus, with the sea urchins, so none lies inside the family.
3. **Let the lineage check flag on its own at class level or above.** Now 91 of 120 name-only errors are listed for review and left in place. On this run the lineage check at that level has no false positive in 1,218 genera. Keep the two-check rule for differences below class level.
4. **Read `annotations.json` for the nodes needed only.** Brings the peak from 1.28 to about 1.05 GB, under the declared limit. Low priority.
5. **Leave these switched as they are**: `requery.mode: classify`, `translation.action: flag`, `rate_outlier_action: flag`, `mutual_pairs: both`. Nothing here measured them end to end. `mutual_pairs: further` has the best case of the four (#55), and the rise in units on the wrong side of the root is a reason to test it with rooting as the metric.
6. **Tighten the gate before it is used again**: replace 6b with the stricter reading; state memory as an absolute limit; add rooting (units on the wrong side of the root) as a metric; add the two injection variants as declared metrics.

Follow-ups 1 to 3 change the result only for flagged genera, so they can go in after the merge without changing the Iberian tree.

## What this benchmark cannot tell

- **One data set, one real error.** Everything is the Iberian COI run. The injected errors are of one kind, chosen to match that one real case, and the variants show how much the kind matters.
- **TimeTree is the reference for ages.** The hold-out and the conflict metrics measure agreement with TimeTree, not with the truth. The robustness results do not have that weakness: they compare each setting with itself.
- **The injections and perturbations stop at the backbone.** They were not carried through the genus subtrees (26 minutes per run). Stem ages bound crown ages, so the direction is clear, but the size of the effect on the final tree was not measured.
- **Blind and known.** Nine of the thresholds were set knowing the stage-level result. The three failures and the robustness results were blind.
- **Not tested per part**: which part produces the robustness, and the options left off at the default.
