# Assumptions and limitations

What timegraft takes for granted, and where it is known to be weak. Evidence comes from the [first run on real data](results/2026-11-iberia.md). Each weakness links to the issue tracking it.

## Assumptions

- **The backbone is right.** Branching order above genus level comes from the Open Tree of Life synthesis and ages from TimeTree, and both are taken as given. Phylogenetic diversity is dominated by deep branches, so errors there matter more than errors within genera.
- **Genera are natural groups.** One tree is built per genus and grafted onto one tip. Genera that are not natural groups (not monophyletic) break this.
- **Barcodes are informative within a genus.** The genus is used as the default range in which barcode distances track evolutionary distance. The validation step can measure this per group.
- **A roughly constant rate of change within each genus.** Both dating methods (mean path length, LSD2) assume clock-like evolution within a genus.
- **A known COI rate** (0.0115 substitutions per site per million years per lineage, about 2.3% pairwise) where no crown age is known. It also sets how far back placed sample sequences attach, and the threshold for rejecting TimeTree crown ages.
- **One representative per BIN** stands in for the variation within it.
- **The genetic code follows from the taxonomy.** The reading-frame check on COI uses a mitochondrial code per record, looked up from its genus, family, order, class or phylum (see [configuration](../technical/configuration.md#genetic-codes)). Records whose names match nothing get a default code, and misidentified records get the code of the wrong group.

## Known weaknesses

| Weakness | Evidence | Issue |
| --- | --- | --- |
| TimeTree covers marine invertebrates poorly, so many backbone nodes are interpolated, and even spacing tends to make shallow undated nodes too old | 261 of 779 backbone nodes dated directly in the Iberian run; *Hediste*/*Nereis* at 291 Ma in a small test | [#16](https://github.com/naturalis/timegraft/issues/16), [#17](https://github.com/naturalis/timegraft/issues/17) |
| TimeTree crown ages are inflated for genera that are not natural groups; the check that rejects them depends on the assumed COI rate | 179 of 305 TimeTree crown ages rejected in the Iberian run (e.g. *Haliclona* at 444 Ma) | [#18](https://github.com/naturalis/timegraft/issues/18) |
| The Open Tree of Life contains misplacements, which nothing corrects automatically | the gastropod *Onchidella* among the sea urchins, because OTT merged it with an echinoderm homonym | [#20](https://github.com/naturalis/timegraft/issues/20) |
| A genus and its relatives often overlap in barcode distance, so some units end up on the wrong side of the root. Rooting now tolerates this, but the root of a genus remains uncertain where many units are misplaced | 534 of 1,046 genera have units on the wrong side (median 2); 579 units of 351 genera sit among the outgroup, some of them probably misidentified. The spread is now recorded per genus: 254 genera have more than one plausible root, and the crown age varies by more than 25% over those roots in 149 | [#19](https://github.com/naturalis/timegraft/issues/19), [#34](https://github.com/naturalis/timegraft/issues/34) |
| Placement by nearest neighbour is less accurate than likelihood-based placement, and its thresholds are not calibrated | at 80% identity the attachment age is already about 9 Ma | [#21](https://github.com/naturalis/timegraft/issues/21) |
| BINs and UNITE species hypotheses are delimited differently | PD across markers is comparable in time but not in the granularity of the tips | [#15](https://github.com/naturalis/timegraft/issues/15) |

## What timegraft is not

It is not a tree of life inferred from barcodes. The deep structure is borrowed from published syntheses; barcodes contribute only the shallow structure and the placement of samples. Results should be presented that way, with the provenance tables that record where every age came from.
