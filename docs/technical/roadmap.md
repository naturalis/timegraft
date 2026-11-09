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

## Next: getting the ages right

The weakest parts of the first real run were the ages, so they come first.

| Issue | Topic |
| --- | --- |
| [#18](https://github.com/naturalis/timegraft/issues/18) | Calibrate the COI substitution rate, which sets most genus crown ages |
| [#16](https://github.com/naturalis/timegraft/issues/16) | Better ages for backbone nodes that TimeTree does not date |
| [#17](https://github.com/naturalis/timegraft/issues/17) | Literature crown ages for groups that TimeTree covers poorly |
| [#20](https://github.com/naturalis/timegraft/issues/20) | Errors in the OToL backbone, and a way to correct them |
| [#19](https://github.com/naturalis/timegraft/issues/19) | Rooting failures in one genus in five |

## Next: the coastal marine invertebrate application (CBMA PhD project)

| Issue | Topic |
| --- | --- |
| [#23](https://github.com/naturalis/timegraft/issues/23) | Build the tree from the genera found in the samples |
| [#24](https://github.com/naturalis/timegraft/issues/24) | Validation against independent reference trees: how deep does COI carry signal, per phylum |
| [#21](https://github.com/naturalis/timegraft/issues/21) | Calibrate placement thresholds; compare with EPA-ng |
| [#22](https://github.com/naturalis/timegraft/issues/22) | Genetic code per taxon |

## Then: fungal ITS (eDentity)

[#15](https://github.com/naturalis/timegraft/issues/15) collects what is needed: a UNITE reader with species hypotheses as units, ITSx extraction, a check of k-mer distances within genera, crown ages without a rate, a fungal backbone, coarse placement via 5.8S, and test data.

## Later

| Issue | Topic |
| --- | --- |
| [#25](https://github.com/naturalis/timegraft/issues/25) | Embeddings as a distance method and for nearest-neighbour search |
| [#26](https://github.com/naturalis/timegraft/issues/26) | Smaller follow-ups: speed of the representative step, genus conflicts, cache, TimeTree licence, treePL |

## Writing it up

A [preprint](../manuscript/manuscript.md) describing the method and the first application is being drafted alongside the code.
