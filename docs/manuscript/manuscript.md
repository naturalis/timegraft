---
title: "timegraft: dated phylogenies of DNA barcode reference libraries for phylogenetic diversity of metabarcoding samples"
author:
  - name: Rutger A. Vos
    affiliation: Naturalis Biodiversity Center, Leiden, The Netherlands
# TODO: agree the author list (CBMA co-authors, DCC contributors)
date: "Draft, November 2026"
bibliography: references.bib
link-citations: true
abstract: |
  Phylogenetic diversity measures summarise how much evolutionary history a community holds and how much two communities share, but applying them to metabarcoding data requires a tree that covers the reference units samples are identified to, with branch lengths that mean the same thing throughout. For the large barcode reference libraries of animals (BOLD) and fungi (UNITE), no such tree exists, and attempts to infer one from the barcodes themselves fail at the deep end: short markers saturate, cannot be aligned across distant groups, or both. We present timegraft, a workflow that takes the topology and ages of the tree down to genus level from published syntheses (Open Tree of Life and TimeTree) and uses barcodes only within genera, where they carry reliable signal. Each genus is resolved and dated independently from its reference units, and grafted onto the dated backbone, so that the final tree is ultrametric, in millions of years, with a BIN or species hypothesis at every tip. Metabarcoding sequences are then matched to reference units or placed near their closest one, and Faith's phylogenetic diversity, mean pairwise distance, UniFrac and PhyloSor are computed per sample. Every age in the tree is traceable to its source. Applied to COI reference data for 1,232 genera of western Iberian invertebrates and fishes (25,578 BINs), timegraft built a dated tree of 24,894 tips in under half an hour on two processor cores. TimeTree dated about a third of the backbone nodes directly; its genus crown ages were often inflated for genera that are not monophyletic, and conflicts between the two syntheses exposed errors in the Open Tree taxonomy. [TODO: summarise the validation and the CBMA application once available.] timegraft is open source and marker-agnostic by design; an extension to fungal ITS is under way.
---

<!--
Build: pandoc manuscript.md --citeproc -o manuscript.pdf   (see README.md in this folder)
Status: first draft. Numbers in Results come from the first run on real data
(docs/science/results/2026-11-iberia.md) and will change as the method is calibrated.
-->

# Introduction

Metabarcoding has made it routine to list the taxa in an environmental sample, but lists alone ignore how those taxa are related. Two samples may share no species and yet consist of close relatives throughout, or share many species and differ deeply in the rest. Phylogenetic diversity measures address this by summarising a community on a tree: Faith's phylogenetic diversity (PD) as the total branch length connecting its members [@faith1992pd], mean pairwise distance (MPD) as the average divergence between them [@webb2002communities], and UniFrac and PhyloSor as the branch length two communities share [@lozupone2005unifrac; @lozupone2007weighted; @bryant2008phylosor]. These measures are only as meaningful as the tree's branch lengths, and they can only be computed for the taxa the tree contains.

For metabarcoding of animals and fungi, the reference units that sequences are identified to are clusters in barcode libraries: Barcode Index Numbers (BINs) in BOLD for the cytochrome *c* oxidase subunit I gene (COI) [@ratnasingham2013bin], and species hypotheses in UNITE for the internal transcribed spacer (ITS) [@koljalg2013sh; @nilsson2019unite]. There are hundreds of thousands of each, most without a place in any published phylogeny.

Building a tree for these units from the barcodes themselves is attractive, since the data are already there, but it fails in predictable ways. COI saturates: beyond a few tens of millions of years, multiple substitutions at the same sites stop pairwise distances from growing, so barcodes cannot resolve or date deep splits. ITS cannot even be aligned across distant fungal groups. Earlier work at our institute ran into exactly these limits. A divide-and-conquer approach that inferred family-level trees and a backbone of exemplars from COI [@bactria], following the logic of earlier large-tree pipelines [@antonelli2017supersmart], did not scale at the backbone, and produced backbone branch lengths that could not be reconciled with those within families. An approach that used embeddings from a DNA language model to resolve a taxonomic tree and assign its branch lengths [@evospaice; @wei2024branchlengths] showed that embedding distances tracked sequence divergence only at short distances; beyond roughly genus level, sequence distance itself had plateaued (Fig. 1).

Elsewhere, the most successful large reference trees for single-gene data have separated the two problems: the deep structure comes from a tree inferred from more informative data, and the single gene is only used to attach sequences to it. Greengenes2 placed more than 21 million 16S sequences onto a genome-based tree of bacteria and archaea [@mcdonald2024greengenes2], using a learned placement model [@jiang2023depp]. For insects, synthesis trees combine phylogenomic backbones with barcodes [@chesters2020insects; @chesters2026lepidoptera]. Placement of COI on fixed backbones is accurate for most species when the backbone is well sampled [@fernando2025placement].

timegraft applies this separation to barcode reference libraries, with two further choices. First, the backbone is not inferred at all but taken from two published syntheses: the branching order of genera from the Open Tree of Life [@hinchliff2015otol; @rees2017ott] and their divergence times from TimeTree [@kumar2022timetree5]. Second, all branch lengths are expressed in time, so that subtrees built independently fit the backbone without rescaling, and diversity is comparable across markers. Barcodes are used where they are informative: to resolve and date relationships within genera, and to place metabarcoding sequences.

Here we describe the method, its implementation as a Snakemake workflow [@molder2021snakemake], and a first application to COI reference data for the genera recorded in western Iberia, and we discuss what the sources of ages contribute and where they fail.

# Methods

## Overview

timegraft builds a dated tree in four stages (Fig. 2): (i) one representative sequence per reference unit; (ii) a dated backbone with one tip per genus; (iii) one dated subtree per genus; and (iv) grafting of the subtrees onto the backbone. Metabarcoding sequences are then matched or placed, and diversity computed. Steps that depend on the marker (filtering, distances, crown ages, dating and placement settings) are configurable plug-ins; COI and ITS settings are provided. The workflow is implemented in Python with DendroPy [@sukumaran2010dendropy] for tree handling and Snakemake for orchestration.

## Reference units

Records are read from a BOLD data package or a table of record identifier, unit, genus and sequence. COI records are kept if they meet a minimum length, contain few ambiguous bases, and have a reading frame without stop codons under the appropriate mitochondrial genetic code. Because a taxonomically broad data set spans several codes, the code is chosen per record: the most specific of the record's genus, family, order, class or phylum that occurs in a mapping derived from the NCBI Taxonomy decides it, with a default for records that match nothing. For each unit, the representative is the medoid record by *k*-mer distance [@ondov2016mash], which requires no alignment. A unit is assigned to the genus that at least 80% of its identified records agree on; units without a genus, or with conflicting identifications, are set aside for placement.

## Dated backbone

Genera are matched by name to the Open Tree Taxonomy within a context taxon (for example Metazoa) to avoid homonyms across kingdoms. The topology is induced from the labelled synthetic tree: a genus that is a node in the synthesis becomes a tip there; otherwise it is attached at the most recent common ancestor (MRCA) of its species in the synthesis, or failing that at its nearest enclosing taxon that is.

Each internal node is dated by asking TimeTree for the age of the MRCA of a few genera from each side of it, keyed by NCBI taxon identifiers. Because TimeTree knows many species but not every genus, each genus is represented together with some of its species; other genera and species are tried when TimeTree knows none of the first picks on one side, and as a last resort the enclosing taxa of a side (from family up to phylum) are used, but only taxa that contain no genus from the other side, so that the MRCA is still the same node. An age is accepted only if TimeTree found taxa on at least two sides of the node.

When a node is not younger than its nearest dated ancestor, equal ages are read as TimeTree not resolving the node, a direct age takes precedence over one obtained through enclosing taxa, and otherwise the ancestor's age is kept. Remaining undated nodes are spaced evenly between their nearest dated ancestor and oldest dated descendant, as in BLADJ [@webb2008phylocom]. Each genus tip then has a stem age *S*.

## Genus subtrees

For each genus with at least two units, the units and up to two units from the nearest genera in the backbone are aligned with MAFFT [@katoh2013mafft], and pairwise distances computed under the Kimura two-parameter model [@kimura1980k2p] (for ITS, *k*-mer distances without alignment). A BioNJ tree [@gascuel1997bionj] is rooted on the outgroup, which is then removed; if the ingroup is not a clade with respect to the outgroup, it is rooted at its midpoint instead.

The crown age *C* of the genus is taken, in order of preference, from a user-supplied table (for example literature values); from TimeTree, as the MRCA of up to twelve of the genus's species, provided this is younger than *S* and at most three times the barcode-based estimate below; or from the barcodes, as the mean root-to-tip distance divided by a substitution rate (by default 0.0115 substitutions per site per million years per lineage, corresponding to the commonly used 2.3% pairwise divergence per million years for arthropod mitochondrial DNA [@brower1994]). *C* is capped at 0.9 *S*. The subtree is made ultrametric with its root at *C* by mean path length [@britton2002mpl] or, optionally, by least-squares dating with LSD2 [@to2016lsd].

## Grafting

Each genus tip of the backbone is replaced by a stem of length *S* − *C* followed by its dated subtree, which keeps the tree ultrametric. Genera with a single unit keep their tip, relabelled with the unit. Genera without units are pruned. There is no global tree search, and genera are processed independently and in parallel.

## Placement and diversity

Units without a genus and metabarcoding amplicon sequence variants (ASVs) are compared with the reference records by global alignment with VSEARCH [@rognes2016vsearch], which accommodates fragments such as the Leray region of COI [@leray2013primers]. An ASV with at least 97% identity to a unit is counted as that unit; equally good matches to several units attach at their MRCA. Below that, down to 80% identity, the sequence is attached on the path from its nearest unit to the root, at the age implied by its divergence and the substitution rate (COI) or halfway along the nearest unit's terminal branch (ITS). Faith's PD, MPD (presence and abundance-weighted), unweighted and weighted UniFrac, and PhyloSor are computed from a single edge-by-sample summary of the tree.

## Validation of the informative window

To measure how deep a marker carries signal, pairs of units are sampled per taxonomic stratum (same genus, same family, same order, deeper), and each barcode distance is compared with the patristic distance in a reference tree inferred from independent data, by rank correlation per stratum and by the reference distance at which the binned barcode distance reaches 90% of its maximum.

## First application

We applied timegraft to COI reference data for the genera recorded in western Iberia. From the BOLD public data package of March 2026 we selected records of Annelida, Mollusca, Arthropoda, Echinodermata, Cnidaria, Bryozoa, Porifera and Chordata from Portugal and Spain, excluding insects, arachnids, springtails, millipedes, centipedes and tetrapods, to obtain a list of genera, and then all COI-5P records of those genera worldwide, at most ten per BIN. The backbone was built from Open Tree of Life synthesis 16.1 and OTT 3.7.3 with Metazoa as context. The run used two processor cores.

# Results

The data comprised 102,794 records of 25,578 BINs in 1,232 genera; the largest genera were the groundwater amphipod *Niphargus* (1,035 BINs) and *Gammarus* (907). The full workflow completed in 20 to 28 minutes. The final tree had 24,894 tips, was ultrametric with its root at 750 million years, and comprised 1,046 genus subtrees and 137 single-unit genera. Of 215 units without a usable genus, 212 were placed next to their nearest unit.

Of the 1,218 genera that kept units after filtering, 996 were nodes in the synthetic tree, 149 were attached at the MRCA of their species and 38 at an enclosing taxon; 35 could not be matched unambiguously in the taxonomy. TimeTree dated 261 of 779 backbone nodes directly and 53 through enclosing taxa (Table 1). A further 57 nodes received the same age as their dated ancestor, meaning TimeTree does not resolve them, and 76 TimeTree ages conflicted with an ancestor's age and were dropped; 332 nodes were interpolated.

Genus crown ages came from TimeTree for 131 genera and from the barcode rate for 915, including 174 for which a TimeTree estimate was rejected. Rejected estimates reached hundreds of millions of years (*Haliclona* 444, *Stenasellus* 289, *Margaritifera* 230 million years), whereas accepted TimeTree crown ages had a median of 12.4 million years against 9.5 million years for rate-based ones. In 194 of the 1,046 genera (19%), the outgroup fell inside the ingroup and midpoint rooting was used.

Conflicts between the two syntheses pointed to errors in either. The largest involved flatfishes, crabs and freshwater snails; one was an outright taxonomic error, in which the gastropod genus *Onchidella* was placed among the sea urchins because two homonymous records had been merged in the taxonomy.

[TODO: Table 1, sources of backbone node ages and genus crown ages.]
[TODO: validation results per phylum, once reference trees are in place.]
[TODO: diversity results for the CBMA samples.]

# Discussion

[TODO: develop. Points to make:]

- **Separating deep and shallow structure works at scale.** Building a tree of 25,000 BINs takes half an hour without any global tree search, and the result is in a unit that is meaningful across markers.
- **The ages are the weak point, not the topology within genera.** TimeTree's coverage of marine invertebrates is patchy, interpolation is crude, and crown ages from TimeTree must be checked because non-monophyletic genera inflate them. Literature crown ages for key groups and a calibrated substitution rate are the obvious next steps.
- **Syntheses check each other.** Disagreements between Open Tree topology and TimeTree ages flagged real errors, which suggests a use for timegraft beyond diversity: auditing both resources.
- **Provenance matters.** Because every age is recorded with its source, users can see which parts of a diversity estimate rest on direct calibrations and which on interpolation or an assumed rate.
- **Limitations.** Genera that are not monophyletic violate the assumption of one subtree per genus; distance-based placement is less accurate than likelihood-based placement [@fernando2025placement; @matsen2010pplacer; @barbera2019epang; @balaban2020apples]; BINs and UNITE species hypotheses are delimited differently, so tips are not equivalent across markers.
- **Extension to ITS** requires extraction of ITS1 and ITS2 [@bengtsson2013itsx], alignment-free distances within genera, and crown ages without a reliable substitution rate; existing placement tools for fungi rely on alignable loci [@carbone2017tbas].

# Data and code availability

timegraft is available at <https://github.com/naturalis/timegraft> [@timegraft] under the Apache 2.0 licence; trees and tables it produces are released under CC0 1.0. [TODO: archive the release and the Iberian tree on Zenodo and cite the DOIs.]

# Acknowledgements

[TODO]

# Figures

**Figure 1.** Saturation of COI. Patristic distance in a tree built from barcode embeddings against Kimura two-parameter distance for 1,000 pairs of BINs; K2P distance plateaus around 0.13 beyond a tree distance of about 0.15 to 0.2. [TODO: regenerate against an independent reference tree with timegraft's validation step.]

**Figure 2.** Overview of timegraft: topology and ages above genus level from published syntheses, a dated subtree per genus from barcodes, grafting, and placement of metabarcoding sequences. [TODO: draw from the flowcharts in docs/science/how-it-works.md.]

**Figure 3.** Sources of ages in the Iberian tree, by depth. [TODO]

# References
