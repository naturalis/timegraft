# Prior art

timegraft combines three ideas that each exist elsewhere: a fixed reference tree onto which single-gene data are attached, a dated synthesis of published phylogenies, and divide-and-conquer tree building in which clades are built separately and joined on a backbone. This page lists the closest related work and how timegraft differs. Earlier attempts at Naturalis itself (Bactria, the MDDB/ITS work, evospaice) are described under [background](background.md).

## Fixed backbones extended with single genes

- **Greengenes2** ([McDonald et al. 2023](https://www.nature.com/articles/s41587-023-01845-1)) is the closest precedent: a genome-based tree of bacteria and archaea kept fixed, with full-length 16S sequences inserted and over 21 million short 16S sequences placed onto it using DEPP. timegraft follows the same architecture for BOLD and UNITE, but with a dated backbone, so that branch lengths are in time rather than substitutions.
- **DEPP** ([Jiang et al.](https://doi.org/10.1101/2021.01.22.427808)) learns an embedding in which distances between sequences of one gene match (the square root of) distances on a reference species tree, and uses it to place new sequences. It needs a reference tree to learn from; in timegraft it would serve clade recognition for placement, not the estimation of deep distances.
- **Phylogenetic placement** tools (pplacer, EPA-ng, APPLES) attach query sequences to a fixed tree, by likelihood or by distance. On fishes, EPA-ng placed 70 to 78% of COI sequences correctly and the distance-based APPLES did worse ([Fernando et al. 2025](https://doi.org/10.1002/ece3.70817)), which sets expectations for the distance-based placement in timegraft.
- **T-BAS** ([Carbone et al. 2017](https://doi.org/10.1093/bioinformatics/btw808)) places fungal sequences on curated multilocus backbones. It relies on alignable loci (5.8S, LSU, SSU and protein-coding genes) and leaves the variable part of ITS out.

## Synthesis trees and dated trees

- **Open Tree of Life** synthesises published phylogenies and a reference taxonomy (OTT) into one tree of all named life. timegraft takes its genus-level topology from the synthesis release files.
- **TimeTree 5** ([Kumar et al. 2022](https://doi.org/10.1093/molbev/msac174)) synthesises published divergence times, with a REST API keyed by NCBI taxon ids. timegraft takes node ages and genus crown ages from it. Its coverage of marine invertebrates turned out to be patchy (see the [first results](results/2026-11-iberia.md)).
- **Chesters' insect synthesis trees** ([Chesters 2020](https://doi.org/10.1111/syen.12414); [insectphylo.org](https://fis.leuphana.de/en/publications/launching-insectphyloorg-a-new-hub-facilitating-construction-and-/)) combine phylogenomic backbones, mitogenomes and barcodes for tens of thousands of insect species, including a Lepidoptera synthesis applied to PD of caterpillar barcodes. They are alignment-based, species-level and limited to insects, but they are a candidate backbone for insects and the most direct comparison for COI.
- **Megatree tools for plants and other groups** take a dated backbone and add missing species at genus or family level, as polytomies or randomly within their clade. timegraft does the same at genus level, but resolves and dates the genus from barcode data rather than adding units without structure.

## Divide and conquer on a backbone

- **SUPERSMART** (Antonelli et al. 2017, *Systematic Biology*) estimates large dated phylogenies by decomposing them into clades that are inferred separately and grafted onto a backbone of exemplars. Bactria followed the same logic for barcodes. timegraft keeps the decomposition but stops estimating the backbone from the same data; the backbone comes from synthesis instead.
- **Branch lengths on a fixed shape** (Wei & Koslicki, [preprint](https://www.biorxiv.org/content/10.1101/2024.07.29.605688v2)): bottom-up assignment of branch lengths to a taxonomy from leaf distances, used in evospaice.

## How timegraft differs

| | Greengenes2 | Chesters' syntheses | Bactria | timegraft |
| --- | --- | --- | --- | --- |
| Group | bacteria, archaea | insects | animals (COI) | animals (COI), fungi (ITS) planned |
| Tips | sequences (ASVs) | species | BOLD records or BINs | BINs or UNITE SHs |
| Backbone | inferred from genomes | inferred from phylogenomics | inferred from COI exemplars | taken from OToL and TimeTree |
| Branch lengths | substitutions | substitutions | substitutions | millions of years |
| Barcodes used for | placement | species-level structure | everything | structure and lengths within genera, placement |
| Alignment | 16S alignment | per-locus alignments | profile HMM | within genus only, or none |

What timegraft adds: a dated tree of all reference units of a marker for the genera of interest; usable for markers that cannot be aligned across the tree; rebuildable from public releases; and explicit about where every age comes from.
