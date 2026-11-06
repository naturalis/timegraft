# Background

timegraft builds dated phylogenies of DNA barcode reference units (BINs for COI, UNITE species hypotheses for ITS) so that phylogenetic diversity (PD) can be computed for metabarcoding samples. It is the third attempt at this goal at Naturalis. This document records what the earlier attempts did, what they taught us, and how timegraft relates to work elsewhere.

## The goal

Metabarcoding gives lists of taxa per sample. Comparing samples by those lists ignores how related the taxa are: two samples can share no species and still be close relatives throughout, or share many species and differ deeply in the rest. Phylogenetic diversity measures (Faith's PD, mean pairwise distance, UniFrac, PhyloSor) use a tree instead, but only if its branch lengths mean something, and only if the tree covers the reference units that samples are identified to. No such tree exists for the barcode reference libraries we use (BOLD for animals, UNITE for fungi).

Users at Naturalis and partners:

- **eDentity**: soil metabarcoding, mainly fungal ITS, with COI for invertebrates.
- **A PhD project with CBMA (University of Minho)**: PD and genetic diversity metrics for coastal marine invertebrate communities of western Iberia, from COI metabarcoding (task 1 of the work plan builds on this pipeline).

## Earlier attempts

### Bactria (barcode-constrained-phylogeny)

[Bactria](https://github.com/naturalis/barcode-constrained-phylogeny) builds a COI tree by divide and conquer under an Open Tree of Life (OToL) constraint:

1. BOLD records are split into family-level chunks.
2. Each family is aligned against a COI profile HMM and its tree inferred with RAxML-NG, constrained by OToL.
3. Two exemplars per family are selected, merged into a backbone alignment, and the backbone inferred with RAxML-NG under the OToL constraint.
4. Family subtrees are grafted onto the backbone.

What it taught us:

- **The backbone does not scale.** Step 3 is a single, non-parallel maximum-likelihood search over all families, and it was the step that stalled.
- **Branch lengths of two analyses do not fit together.** Backbone branches came out much longer than the equivalent branches in the family subtrees. The Bactria documentation already attributed this to saturation and homoplasy in the broad taxon sample. Rescaling the two was never solved.
- **Exemplars and subtree rooting are fragile**, partly because of sparse taxon sampling within families.
- **It only works for an alignable marker** with a curated HMM (the one used was built for terrestrial insects).

### MDDB tree and ITS placement (LIACS)

For fungal ITS, an earlier project built a reference tree for the MDDB data, and a LIACS thesis (supervised by Rutger Vos) studied placement of ITS reads on it with a two-step approach: find the right subtree with BLAST, then place within it with pplacer. Lessons that carry over: alignment-free and distance-based approaches are tractable for ITS, but need taxonomic guidance; ITS cannot be aligned across fungi as a whole, so placement has to happen within a subtree.

### evospaice (Microsoft collaboration)

[evospaice](https://github.com/naturalis/evospaice) was a hackathon project with Microsoft (2026) to build a scaled reference tree from barcode embeddings: insect COI (the Leray fragment) from BOLD, embedded with a DNA language model (Omni-DNA-20M; the brief also mentions DNABERT-S, both fine-tuned with a contrastive objective for identification). The taxonomy gave the rough shape of the tree; neighbour joining on embedding distances resolved the bushy nodes and set branch lengths, bottom-up.

The brief was explicit about the main risk: agreement with BLAST shows that an embedding is a good identifier, not that its distances form a good metric over larger distances. The key experiment confirmed it. For 1,000 pairs of tips, tree distance in the embedding tree was compared with K2P distance from an HMM alignment:

- Up to a tree distance of roughly 0.15 to 0.2, both rise together (with considerable scatter).
- Beyond that, K2P plateaus around 0.13 while tree distance keeps growing (Pearson r = 0.45, Spearman ρ = 0.60).

Neither axis was ground truth, so the plot cannot say which one is off, but the conclusion holds either way: COI carries reliable distance information only up to roughly genus or young family level, so no distance derived from it (embedding, k-mer or alignment) can date or resolve deep nodes. On top of that, NJ over large distances, the choice of representatives, and the bottom-up branch length assignment can distort deep branches further.

The collaboration with Microsoft has ended. The repository is public domain and the embeddings computed for about 7 million insect records are on Naturalis's own Azure subscription, so code and data can be reused.

## What we concluded

1. **Do not estimate deep structure from barcodes.** Take topology and ages above the informative window from external sources.
2. **Use barcodes only where they carry signal**, within genera, and measure that window per marker and clade rather than assume it.
3. **Put everything in one unit, time**, so that PD is comparable across markers, clades and studies, and so that subtrees from different analyses fit together without rescaling.
4. **Keep it alignment-free where needed** (ITS) and parallel by construction (every genus independently).
5. **Embeddings are useful for search and coarse placement**; whether they beat k-mers or K2P for distances within a genus is an empirical question for the validation step.

This is the design timegraft implements; see [method.md](method.md).

## Related work

- **Greengenes2** ([McDonald et al. 2023](https://www.nature.com/articles/s41587-023-01845-1)) is the closest precedent: a genome-based tree kept fixed, with over 21 million 16S sequences placed onto it using DEPP. timegraft is in effect the same architecture for BOLD and UNITE, with a dated backbone.
- **DEPP** ([Jiang et al.](https://doi.org/10.1101/2021.01.22.427808)) learns an embedding in which distances match (the square root of) distances on a reference tree, and places single-gene sequences on a species tree. It needs a reference tree to learn from; its use here would be clade recognition for placement, not estimating deep distances.
- **Chesters' insect synthesis trees** ([Chesters 2020](https://doi.org/10.1111/syen.12414); [insectphylo.org](https://fis.leuphana.de/en/publications/launching-insectphyloorg-a-new-hub-facilitating-construction-and-/)) combine phylogenomic backbones with barcodes for tens of thousands of insect species, including a Lepidoptera synthesis applied to PD of caterpillar barcodes. Alignment-based, species-level and insects only, but a candidate backbone for insects and the most direct comparison.
- **T-BAS** ([Carbone et al. 2017](https://doi.org/10.1093/bioinformatics/btw808)) places fungal sequences on curated multilocus backbones; it relies on alignable loci and leaves the variable part of ITS out.
- **COI placement on backbones** ([Fernando et al. 2025](https://doi.org/10.1002/ece3.70817)): on fishes, 70 to 78% of COI sequences were placed correctly by EPA-ng, and the distance-based APPLES did worse. A benchmark for the cost of distance-based placement.
- **TimeTree 5** ([Kumar et al. 2022](https://doi.org/10.1093/molbev/msac174)) synthesises published divergence times and has a REST API keyed by NCBI taxon ids. Its coverage of marine invertebrates turned out to be patchy (see the [first results](results/2026-11-iberia.md)).
- **Wei & Koslicki** ([preprint](https://www.biorxiv.org/content/10.1101/2024.07.29.605688v2)): bottom-up branch length assignment on a taxonomy, used in evospaice.

## What timegraft adds

A dated tree, in millions of years, of all reference units of a marker for the genera of interest; alignment-free where needed and therefore usable for both COI and ITS; regularly rebuildable from public releases; and explicit about where each age comes from (TimeTree, enclosing taxa, interpolation, barcode rate, or a user's table). It is not a tree of life inferred from barcodes and should not be presented as one.
