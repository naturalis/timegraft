# The barcode cap on midpoint crowns (#79): the crown-tail metrics of the gate

[#79](https://github.com/naturalis/timegraft/issues/79) asked why *Marthasterias* has a crown of 20.6 Ma with a mean root-to-tip distance of 0.015, the last genus counted by metric 3d. The answer, in the issue, is that the crown is the geometric midpoint between the barcode age (1.29 Ma) and 0.9 times a stem of 366 Ma. That stem is interpolated in an unresolved node of starfish from two orders. Nothing in the data says the sampled units are old: both are *M. glacialis*, and TimeTree knows none of the genera around it.

This change bounds such a crown by the barcodes as well. Where the stem is the upper bound, the barcode age times the allowance of #54 (`rate_factor` × the saturation factor) bounds it too, and the midpoint is kept below that bound. That does not apply where the units lump species, or where TimeTree or the table bound the crown.

## Method

The crown tails (3a to 3d) are read from the provenance of the full default run, which is the output of its genus stage. The clean default run of the downstream harness reproduces that stage on the same backbone; metric 7c checks the two agree. So that run was repeated with this branch:

```
PYTHONPATH=$B python $B/downstream.py --repo . --backbone $BENCH70/exp/default_clean \
    --units $BENCH73/branch/results/iberia/units --out clean_default --memo memo --jobs 2 --threads 2
```

The run went through to placement: its 1,786 units are placed in the same genera with the same nearest units as before, only their attachment ages follow the crowns. The run before is `exp/downstream/clean/default` of the benchmark of #70, and legacy is the full run of `main` at `e20fdf6`, as in the gate.

The run was made twice. The second run, reported here, followed three corrections from the review of #114:
- the cap is not applied where the table gives the genus only a minimum;
- species are counted on the units of each tree a crown is chosen for, including the alternative root without misplaced units;
- species are reduced to genus and specific epithet, so subspecies and authors do not count twice.

Compared with the first run, 84 genera count fewer species. One, *Echinolittorina* (62 units, 68 species named before and 51 after), is no longer exempt; its upper bound is now the barcode cap, but its crown does not change. Nothing else changed apart from the wording of `crown_decision`.

## Result

`crown_tails.tsv`:

| Metric | Legacy | Default before | Default after | Threshold |
| --- | --- | --- | --- | --- |
| 3a crown more than 20 times the barcode age | 7 | 0 | 0 | 0 |
| 3b the same above 10 times | 40 | 26 | 26 | fewer than legacy |
| 3c TimeTree crown used below a third of the barcode age | 5 | 0 | 0 | 0 |
| 3d root-to-tip below 0.02 and crown above 20 Ma | 4 | 1 | **0** | report |

3b does not move. A crown at the barcode cap is 10 to 11 times its barcode age (`rate_factor` 10 times a saturation factor slightly above 1), so it is still counted above 10.

**What changes.** Nine crowns change (`changed_crowns.tsv`), all moved down to the barcode cap (`crown_source` `midpoint_barcode_max`):

| Genus | Units | Species | Barcode (Ma) | Stem (Ma) | Crown before | Crown after |
| --- | --- | --- | --- | --- | --- | --- |
| *Apseudes* | 2 | 1 | 1.80 | 244.4 | 19.9 | 19.3 |
| *Elona* | 2 | 1 | 0.74 | 286.3 | 13.9 | 7.7 |
| *Euconulus* | 10 | 5 | 2.73 | 466.4 | 33.9 | 30.7 |
| *Gibbulinella* | 2 | 1 | 2.26 | 466.4 | 30.8 | 24.8 |
| *Gittenbergeria* | 2 | 1 | 1.22 | 286.3 | 17.8 | 12.9 |
| *Heptranchias* | 4 | 1 | 1.19 | 256.3 | 16.6 | 12.5 |
| *Marthasterias* | 2 | 1 | 1.29 | 366.5 | 20.6 | 13.6 |
| *Sphaerocoryne* | 3 | 1 | 1.51 | 276.9 | 19.4 | 16.0 |
| *Xerotricha* | 2 | 1 | 0.88 | 194.4 | 12.4 | 9.1 |

The upper bound (`crown_upper`, `crown_upper_source` `barcodes`) changes for 228 more genera whose midpoint was already below the cap. That narrows the upper PD bound of the diversity step (`faith_pd_crown_max`), but no crown age. The crown range over alternative roots (`crown_age_min`, `crown_age_max`) changes for the nine and for *Myxilla*, one of whose other roots now meets the cap (minimum 16.8 to 9.8 Ma; its crown, 40.9 Ma, is unchanged). Every other column of `provenance.tsv` is identical in all 1,186 genera.

**Lumped units.** Nine genera have more species named among the records of their units than units, and a barcode cap below their stem bound (`exempt_lumped.tsv`). For seven of them the midpoint is below the cap anyway. The two it exempts are *Callogorgia* (an octocoral, 4 species in 3 units; crown 19.5 Ma, would have been 18.5) and *Paratimea* (a sponge, 5 species in 3 units; 14.3 Ma, would have been 9.9). In both groups mitochondrial COI evolves very slowly, so an old crown with little divergence can be genuine there. Species are counted as binomials only (genus and specific epithet, without `sp.`, `cf.` and the like). With that count *Apseudes* and *Xerotricha* have one species each, so they are capped, not exempt, unlike what the comment on the issue expected.

**Why the midpoint is not centred on the cap.** The first version replaced the stem bound by the barcode bound and took the midpoint against it. That changed 185 crowns in the six shards that had finished, each to about 3.2 times its barcode age (√10.5).

It was checked against the 244 genera whose TimeTree crown, of the species sampled, is plausible and below the stem bound (`validation.tsv`). The midpoints were computed without TimeTree, as if it had no crown for them. For the 67 of them where the cap binds:

| Midpoint | Median factor from TimeTree | Within a factor 2 | Median ratio to TimeTree |
| --- | --- | --- | --- |
| Against the stem bound (before) | 1.58 | 0.69 | 1.17 |
| Against the barcode bound (rejected) | 1.63 | 0.60 | 0.93 |
| Against the stem bound, kept below the barcode bound (this change) | 1.58 | 0.69 | 1.17 |

Centred on the cap, midpoints were further from TimeTree and too young. The barcode bound is a worst case (the slowest rate the allowance admits), not a typical age. Kept below the cap, the midpoint changes none of these 244 genera: none of them has a midpoint above its cap. So the change touches only crowns that the barcodes cannot account for.

## What this does not show

- **The stem itself.** *Marthasterias*'s stem stays at 366 Ma, against an Asteriidae crown younger than about 200 Ma. That is a gap in the TimeTree coverage of the backbone (#78, #108). The crowns here are bounded by the barcodes, not corrected through the stem.
- **The other gate metrics.** Nine crown ages change, so the PD metrics (1o, 5e, 6e, 6f) would move a little. They are left to the full gate run after #108.
