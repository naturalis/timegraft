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
- **Follow-ups of the consistency layer** ([epic #89](https://github.com/naturalis/timegraft/issues/89)), done so far:
  - Correction of misplaced genera: only the genus's own species go to TimeTree after re-anchoring, a moved genus is attached within its target family, found by name where OTT gave its name to another taxon, and the lineage check flags on its own at class level or above ([#69](https://github.com/naturalis/timegraft/issues/69), [#70](https://github.com/naturalis/timegraft/issues/70), [#101](https://github.com/naturalis/timegraft/issues/101), [#71](https://github.com/naturalis/timegraft/issues/71)).
  - Memory: the synthesis annotations are read only for the nodes needed ([#72](https://github.com/naturalis/timegraft/issues/72)).
  - Rooting and screening: the rooting effect of the screening fix measured, the long-branch figures redone, and long branches set aside by default ([#75](https://github.com/naturalis/timegraft/issues/75), [#76](https://github.com/naturalis/timegraft/issues/76), [#73](https://github.com/naturalis/timegraft/issues/73)); the second screening search made cheaper ([#85](https://github.com/naturalis/timegraft/issues/85)).
  - Traceability: what changed each node age is recorded ([#84](https://github.com/naturalis/timegraft/issues/84)); the TimeTree cache of the benchmark is kept for offline reruns ([#87](https://github.com/naturalis/timegraft/issues/87)).
  - Re-querying of conflicting nodes kept as a diagnostic only ([#86](https://github.com/naturalis/timegraft/issues/86)); a midpoint crown age is bounded by the barcodes where only the stem bounds it ([#79](https://github.com/naturalis/timegraft/issues/79)).
  - Backbone calibrations use TimeTree's adjusted time where it contradicts the raw one, and a TimeTree answer only if it dates the intended node ([#108](https://github.com/naturalis/timegraft/issues/108), steps 1 and 2).
  - The release gate, tightened in version 2 and with the PD metrics on paired samples in version 3 ([#74](https://github.com/naturalis/timegraft/issues/74), [#104](https://github.com/naturalis/timegraft/issues/104)). Reruns after the follow-ups pass 27 of 28 *must* metrics; 1g (how far the stems of other genera move around a corrected error) still fails ([report](../../data/iberia/consistency/benchmark/REPORT.md)).
  - Representative selection sped up with a vectorised k-mer distance matrix, bit for bit the same ([#112](https://github.com/naturalis/timegraft/issues/112)).
- Decided against: embeddings as a distance method ([#25](https://github.com/naturalis/timegraft/issues/25)). They separate species but are no better than alignment-based or k-mer distances beyond short distances, and need a GPU.

## Next: finishing the consistency layer

What remains of [epic #89](https://github.com/naturalis/timegraft/issues/89):

| Issue | Topic |
| --- | --- |
| [#118](https://github.com/naturalis/timegraft/issues/118) | Date the backbone by matching TimeTree's tree (step 3 of [#108](https://github.com/naturalis/timegraft/issues/108)), once TimeTree's answers come; whether TimeTree caches may stay in this public repository |
| [#77](https://github.com/naturalis/timegraft/issues/77) | A fresh Iberian reference run from `main` |
| [#80](https://github.com/naturalis/timegraft/issues/80) | Which part of the layer gives the robustness gain |
| [#83](https://github.com/naturalis/timegraft/issues/83) | Revisit the judgement-set thresholds; needs the second data set of [#113](https://github.com/naturalis/timegraft/issues/113) |

A first release follows once there is a minimum viable product; separate release notes for the layer were not needed ([#88](https://github.com/naturalis/timegraft/issues/88)).

## Next: the coastal marine invertebrate application (CBMA PhD project)

| Issue | Topic |
| --- | --- |
| [#23](https://github.com/naturalis/timegraft/issues/23) | Build the tree from the genera found in the samples |
| [#65](https://github.com/naturalis/timegraft/issues/65) | Diversity metrics for the project |
| [#24](https://github.com/naturalis/timegraft/issues/24) | Validation against independent reference trees: how deep does COI carry signal, per phylum |
| [#21](https://github.com/naturalis/timegraft/issues/21) | Calibrate placement thresholds; compare with EPA-ng |
| [#20](https://github.com/naturalis/timegraft/issues/20) | Errors in the OToL backbone: the mechanism is in place (#49, #50); remaining cases and curated overrides |

## Then: fungal ITS (eDentity), and a second data set for the gate

[#113](https://github.com/naturalis/timegraft/issues/113) makes ITS a working, tested and benchmarked marker, and uses it to test whether the consistency layer generalises beyond the Iberian COI run. It merges the earlier issues on extending to ITS and on benchmarking a second data set ([#15](https://github.com/naturalis/timegraft/issues/15), [#82](https://github.com/naturalis/timegraft/issues/82)). In six parts:

1. Reference data: a UNITE reader with species hypotheses as units, ITSx extraction, a chimera check, and fungal test data with an end-to-end test.
2. Distances and rooting: the k-mer size chosen by validation against a reference tree, and a check of outgroups without signal.
3. Backbone, crown ages and placement: coverage of a fungal backbone, how often crown ages come from TimeTree, and coarse placement via 5.8S.
4. Calibration of ITS screening, which is off until then ([configuration](configuration.md#screening-its)).
5. The release gate on the ITS data and, if available, a second COI data set, with the harness freed of its Iberian assumptions.
6. Documentation, including what comparing PD across markers means.

## Writing it up

A [preprint](../manuscript/manuscript.md) describing the method and the first application is being drafted alongside the code. Its methods and numbers predate the consistency layer ([status](../manuscript/README.md#status)).
