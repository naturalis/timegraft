# Roadmap

Where timegraft stands and what comes next. Each item links to its issue, where the details and discussion live.

## Done

- Reference units from record tables and BOLD data packages (TSV and Parquet), with marker filters and one representative per unit
- Dated genus-level backbone from Open Tree of Life release files and TimeTree, with checks on TimeTree crown ages and on conflicting node ages
- One dated subtree per genus (K2P or k-mer distances, BioNJ, outgroup rooting, mean path length or LSD2), grafted onto the backbone
- Placement of units without a genus and of metabarcoding ASVs
- Faith's PD, MPD, UniFrac and PhyloSor
- Validation of barcode distances against a reference tree
- A first run on real data: 25,000 BINs in 1,232 genera ([results](../science/results/2026-11-iberia.md))
- Crown ages as a range, from a barcode minimum to TimeTree and the stem age, with the bounds recorded and diversity reported at both ([#36](https://github.com/naturalis/timegraft/issues/36); [results](../science/results/2026-11-dating.md)); penalised likelihood and rates fitted to TimeTree tested and not adopted
- Screening of reference units that do not fit their genus, before genus subtrees are built ([#38](https://github.com/naturalis/timegraft/issues/38); [results](../science/results/2026-11-screening.md))

## Next: getting the ages right

The weakest parts of the first real run were the ages, so they come first.

| Issue | Topic |
| --- | --- |
| [#18](https://github.com/naturalis/timegraft/issues/18) | A rate cannot be fitted to these data (distances saturate within genera). Decide whether the estimate within the bounds becomes the midpoint; set rates per class from the literature where they are known to differ (elasmobranchs) |
| [#16](https://github.com/naturalis/timegraft/issues/16) | Better ages for backbone nodes that TimeTree does not date |
| [#17](https://github.com/naturalis/timegraft/issues/17) | Literature crown ages for groups that TimeTree covers poorly: the table and its format are in place (`data/crown_ages/literature.tsv`, three genera); it needs filling, with the PhD project |
| [#20](https://github.com/naturalis/timegraft/issues/20) | Errors in the OToL backbone, and a way to correct them |
| [#19](https://github.com/naturalis/timegraft/issues/19) | Rooting failures in one genus in five |

## Next: the coastal marine invertebrate application (CBMA PhD project)

| Issue | Topic |
| --- | --- |
| [#23](https://github.com/naturalis/timegraft/issues/23) | Build the tree from the genera found in the samples |
| [#24](https://github.com/naturalis/timegraft/issues/24) | Validation against independent reference trees: how deep does COI carry signal, per phylum |
| [#21](https://github.com/naturalis/timegraft/issues/21) | Calibrate placement thresholds; compare with EPA-ng |

## Then: fungal ITS (eDentity)

[#15](https://github.com/naturalis/timegraft/issues/15) collects what is needed: a UNITE reader with species hypotheses as units, ITSx extraction, a check of k-mer distances within genera, crown ages without a rate, a fungal backbone, coarse placement via 5.8S, and test data.

## Later

| Issue | Topic |
| --- | --- |
| [#25](https://github.com/naturalis/timegraft/issues/25) | Embeddings as a distance method and for nearest-neighbour search |
| [#26](https://github.com/naturalis/timegraft/issues/26) | Smaller follow-ups: speed of the representative step, genus conflicts, cache, TimeTree licence, treePL |

## Writing it up

A [preprint](../manuscript/manuscript.md) describing the method and the first application is being drafted alongside the code.
