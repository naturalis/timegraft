# Screening reference units (November 2026)

Crown ages rest on the barcode distances within a genus ([crown ages](2026-11-dating.md)). [#38](https://github.com/naturalis/timegraft/issues/38) asked how many units do not belong in the genus they are assigned to, and what they do to those distances. This page reports what an all-against-all comparison of the Iberian units found, the screening step built on it, and what that step did to the Iberian tree. Scripts are in [`data/iberia/screening`](../../../data/iberia/screening/README.md).

In short: most units fit their genus, gross contaminants are rare, and about 5% of units are set aside. Screening makes rooting much more stable and removes inflated crown ages in about a hundred genera, but it does not change the overall picture of saturation or the agreement with TimeTree. The data are not bad in general; the age problem is in the marker, not in the identifications.

## What the comparison found

All 24,882 unit representatives were searched against each other with VSEARCH (global alignment, identity at least 80%, query coverage at least 70%; two minutes on two cores). For units of genera with at least two units:

| Nearest relative in the data | Units | Share |
| --- | --- | --- |
| A unit of the same genus | 22,510 | 91.0% |
| Another genus, by less than 2% | 602 | 2.4% |
| Another genus, by 2% or more | 1,043 | 4.2% |
| Nothing at 80% identity | 583 | 2.4% |

- **Gross contaminants are rare.** Twelve units are at least 95% identical to a unit of another order or class: a *Balanus* BIN at 95% to the octocoral *Leptogorgia*, *Radix* and *Gyraulus* at 99%, *Haliclona* and *Halichondria* at 98%, and a few pairs of sea slug and polychaete genera.
- **Most units closer to another genus share a family with it**, at a median identity of 86%. That points to genera that are not natural groups, or to misidentifications within a family, more than to contamination.
- **Single-record BINs are over-represented** among suspect units (57% of those set aside, 35% overall), as expected if some are poor sequences or one-off misidentifications.
- **The rooting step only caught part of this.** Of the 579 units it had flagged as falling among the outgroup, 398 are also suspect here.

## The screening step

`screen_units` runs between `representatives` and `shard_genera` ([configuration](../../technical/configuration.md#screening-units)). A unit is set aside when it is at least 2 percentage points more similar to a unit of another genus than to any unit of its own genus (`other_genus`), or when no unit of its own genus reaches 80% identity while other units of the genus do (`isolated`). A hit counts for both units, because heuristic searches do not always find a pair both ways. Genera with one unit are not screened, and a genus never loses all its units.

Set-aside units are not discarded. They are placed like units without a usable genus: next to their nearest unit at the age implied by their distance, or, if nothing reaches the placement threshold, on the stem of their own genus halfway between its crown and stem ages.

## The Iberian run with screening

| | Before | After |
| --- | --- | --- |
| Units set aside | | 1,294 in 508 genera: 850 `other_genus`, 444 `isolated` |
| Units kept because their whole genus would have gone | | 60 |
| Set-aside units placed next to their nearest unit | | 819 |
| Set-aside units on the stem of their own genus | | 474 |
| Genera with a subtree | 1,046 | 1,036 (10 now have a single unit) |
| Genera with units on the wrong side of the root | 534 | 366 |
| Units of a genus among its outgroup | 579 | 276 |
| Crown age varies by more than a quarter over plausible roots | 160 genera | 50 |
| ... by more than a factor of two | 37 | 17 |
| Median crown age | 8.9 Ma | 8.2 Ma |
| Total branch length | 329,193 Ma | 349,277 Ma |

- **Rooting is much cleaner.** Fewer genera have units on the wrong side of the root, and the crown age depends on the choice of root in a third as many genera.
- **About a hundred genera change.** In the 488 genera that lost units, mean root-to-tip distance fell by a median of 7%, by more than a quarter in 96 genera and by more than half in 32; their crown ages fell accordingly (*Balanus* from 10.5 to 7.9 Ma). In some genera the distance rose slightly, because the root moved once the outliers were gone.
- **Agreement with TimeTree does not improve.** For 260 genera with a TimeTree crown age below the stem age, the barcode minimum is at a median of 0.32 of TimeTree before and after, and the midpoint within a factor of two for 75% before and 73% after. Saturation is a property of the marker, not of a few bad units: removing them leaves the relation between distance and crown age as it was.
- **The total branch length rises by 6%**, from where the set-aside units are attached. The 474 on genus stems sit halfway between crown and stem ages, which gives each a long terminal branch.

## What it does to sample diversity

Faith's PD of the same random samples on the trees before and after screening (200 samples per design):

| Sample | PD after over before | Rank correlation of samples |
| --- | --- | --- |
| 150 units at random | 1.00 | 0.98 |
| 100 genera, one unit each | 1.00 | 0.98 |
| 100 genera, five units each | 1.08 | 0.90 |
| 100 genera, five units each, set-aside units excluded | 1.00 | |

The change for congener-rich samples comes entirely from the placement of set-aside units, not from the crown ages of the genera that kept their units. Where a set-aside unit occurs in a sample with its congeners, the halfway attachment on the stem adds a long branch. This is a convention, like the estimate within the crown age bounds: the unit diverged before the crown of the rest of its genus, so it belongs somewhere on the stem, but where on the stem is not known. `placement.stem_fraction` sets it.

## What this means

- The data are mostly sound. Bad identifications explain inflated crown ages in about a hundred genera, not the general mismatch with TimeTree.
- Screening is worth keeping for the stability of rooting and for those genera, at a cost of two minutes per run.
- The thresholds (2% and 80%) are judgements. They were not tuned to make the results look better; they were set before the rerun, from the comparison above.
- Screening cannot tell a misidentified unit from a genus that is not a natural group. Both are set aside in the same way, which is the right thing for crown ages, but the second calls for a taxonomic fix, not a data fix. `units/screened.tsv` lists the evidence for each decision for those who want to follow it up.
- It is off for ITS until the thresholds are tested on fungal data.

Follow-up in [#55](https://github.com/naturalis/timegraft/issues/55): every decision now carries its reason and thresholds, and an injection benchmark measures how many mislabelled units the thresholds catch, a codon comparison flags likely nuclear copies, a relative-rate check lists long branches within genera, and an external reference can be searched alongside the run. See the [record](../../../data/iberia/consistency/README.md#screening-55) and [configuration](../../technical/configuration.md#screening-units), which also sets out how to calibrate screening for ITS.
