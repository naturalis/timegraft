# Roadmap

Where timegraft stands and what comes next. Each item links to its issue, where the details and discussion live.

## Done

- Reference units from record tables and BOLD data packages (TSV and Parquet), with marker filters, a genetic code per record ([#22](https://github.com/naturalis/timegraft/issues/22)) and one representative per unit
- Dated genus-level backbone from Open Tree of Life release files and TimeTree; undated nodes placed by a pure-birth model on clade sizes in the synthesis ([#16](https://github.com/naturalis/timegraft/issues/16); [results](../science/results/2026-11-interpolation.md))
- One dated subtree per genus (K2P or k-mer distances, BioNJ, tolerant outgroup rooting with the root uncertainty recorded, mean path length or LSD2), grafted onto the backbone ([#19](https://github.com/naturalis/timegraft/issues/19), [#34](https://github.com/naturalis/timegraft/issues/34))
- Crown ages as a range, from a barcode minimum to TimeTree and the stem age, with the bounds recorded and diversity reported at both ([#36](https://github.com/naturalis/timegraft/issues/36), [#42](https://github.com/naturalis/timegraft/issues/42); [results](../science/results/2026-11-dating.md)). Penalised likelihood and rates fitted to TimeTree were tested and not adopted ([#18](https://github.com/naturalis/timegraft/issues/18)); a table for literature crown ages is in place, but filling it systematically was not pursued ([#17](https://github.com/naturalis/timegraft/issues/17))
- Screening of reference units that do not fit their genus, before genus subtrees are built ([#38](https://github.com/naturalis/timegraft/issues/38); [results](../science/results/2026-11-screening.md))
- Placement of units without a genus and of metabarcoding ASVs; Faith's PD, MPD, UniFrac and PhyloSor; validation of barcode distances against a reference tree
- **The consistency layer** ([epic #46](https://github.com/naturalis/timegraft/issues/46), merged in November 2026): checks and corrections between the inputs and the dated tree.
  - Evidence and provenance behind every backbone age ([#47](https://github.com/naturalis/timegraft/issues/47)).
  - The Catalogue of Life as identifier hub ([#48](https://github.com/naturalis/timegraft/issues/48)).
  - Placement concordance checks and re-anchoring of misplaced genera, which fixes *Onchidella* ([#49](https://github.com/naturalis/timegraft/issues/49), [#50](https://github.com/naturalis/timegraft/issues/50)).
  - Majority-clade anchoring ([#51](https://github.com/naturalis/timegraft/issues/51)).
  - A weighted fit of backbone ages ([#52](https://github.com/naturalis/timegraft/issues/52)) and re-queries of conflicting nodes ([#53](https://github.com/naturalis/timegraft/issues/53)).
  - TimeTree crowns of the sampled species with a plausibility test ([#54](https://github.com/naturalis/timegraft/issues/54)).
  - Wider screening with reasons, codon and relative-rate diagnostics ([#55](https://github.com/naturalis/timegraft/issues/55)).
  - A release gate on the Iberian run passed 25 of 28 metrics ([#56](https://github.com/naturalis/timegraft/issues/56); [report](../../data/iberia/consistency/benchmark/REPORT.md)).
- Decided against: embeddings as a distance method ([#25](https://github.com/naturalis/timegraft/issues/25)). They separate species but are no better than alignment-based or k-mer distances beyond short distances, and need a GPU.

## Next: finishing the consistency layer

[Epic #89](https://github.com/naturalis/timegraft/issues/89) collects the follow-ups of the benchmark. The first group closes the metrics that failed:

| Issue | Topic |
| --- | --- |
| [#69](https://github.com/naturalis/timegraft/issues/69) | After re-anchoring, send only the genus's own species to TimeTree |
| [#70](https://github.com/naturalis/timegraft/issues/70) | Re-anchor within the target family, not on it |
| [#71](https://github.com/naturalis/timegraft/issues/71) | Let the lineage check flag on its own at class level or above |
| [#72](https://github.com/naturalis/timegraft/issues/72) | Read the synthesis annotations only for the nodes needed (memory) |
| [#75](https://github.com/naturalis/timegraft/issues/75) | More units on the wrong side of the root after the screening fix |
| [#77](https://github.com/naturalis/timegraft/issues/77), [#88](https://github.com/naturalis/timegraft/issues/88) | A fresh Iberian reference run from `main`, and release notes |

The rest of #89 covers results to check, what the benchmark did not measure, provenance and operations.

## Next: the coastal marine invertebrate application (CBMA PhD project)

| Issue | Topic |
| --- | --- |
| [#23](https://github.com/naturalis/timegraft/issues/23) | Build the tree from the genera found in the samples |
| [#65](https://github.com/naturalis/timegraft/issues/65) | Diversity metrics for the project |
| [#24](https://github.com/naturalis/timegraft/issues/24) | Validation against independent reference trees: how deep does COI carry signal, per phylum |
| [#21](https://github.com/naturalis/timegraft/issues/21) | Calibrate placement thresholds; compare with EPA-ng |
| [#20](https://github.com/naturalis/timegraft/issues/20) | Errors in the OToL backbone: the mechanism is in place (#49, #50); remaining cases and curated overrides |

## Then: fungal ITS (eDentity)

[#15](https://github.com/naturalis/timegraft/issues/15) collects what is needed: a UNITE reader with species hypotheses as units, ITSx extraction, a check of k-mer distances within genera, crown ages without a rate, a fungal backbone, coarse placement via 5.8S, and test data. Screening of ITS units needs calibrating first ([configuration](configuration.md#screening-its)); [#82](https://github.com/naturalis/timegraft/issues/82) proposes doing that, and running the release gate, on eDentity's data.

## Later

| Issue | Topic |
| --- | --- |
| [#26](https://github.com/naturalis/timegraft/issues/26) | Smaller follow-ups: speed of the representative step, genus conflicts, cache, TimeTree licence, treePL |

## Writing it up

A [preprint](../manuscript/manuscript.md) describing the method and the first application is being drafted alongside the code. Its methods and numbers predate the consistency layer ([status](../manuscript/README.md#status)).
