# How timegraft works

timegraft takes topology and ages down to genus level from external sources, builds one dated subtree per genus from barcodes, and grafts the subtrees onto the dated backbone. The result is an ultrametric tree in millions of years with reference units (BINs or UNITE species hypotheses) at the tips. Samples are then matched or placed on it, and diversity is computed.

## Overview

```
records (BOLD or UNITE) ─► filter ─► one representative per unit ─► split by genus ─┐
                                                                                     ├─► per genus: distances, BioNJ, root, date ─► graft ─► place units without a genus ─► timegraft.tre
genera ─► OToL topology ─► TimeTree node ages ─► dated genus backbone (stem ages) ───┘
                                                                                                                                     │
ASVs ─► match to a unit, or place near the nearest unit ─────────────────────────────────────────────────────────────────────────────┴─► Faith's PD, MPD, UniFrac, PhyloSor
```

Steps that differ between markers are plug-ins configured per marker (`markers.coi`, `markers.its` in `config/config.yaml`): the filter, the distance method, the crown age method, the dating method and the placement settings.

## Step by step

### 1. Records and representatives

- **Filter.** COI: minimum length, maximum share of ambiguous bases, and a reading frame free of stop codons under a chosen genetic code. ITS: length and ambiguity (ITSx extraction is expected upstream). No profile HMM is used.
- **Representative per unit.** For every unit, the medoid record by k-mer distance (alignment-free, so it works for any marker), from at most `max_records_per_unit` records.
- **Genus per unit.** The majority genus among a unit's records, if at least `genus_agreement` (default 0.8) of them agree. Units without a genus or with conflicting genera are set aside and placed later (step 5).

### 2. The dated backbone

Either given as a file (a dated, ultrametric Newick tree with genera at the tips) or built:

- **Topology** from OToL release files: genera are matched to OTT within a context taxon (e.g. Metazoa), to resolve homonyms. A genus that is a node in the synthetic tree becomes a tip there; otherwise it is attached at the MRCA of its species in the synthesis, or failing that at its nearest enclosing taxon in the synthesis.
- **Node ages** from TimeTree: for each internal node, a few genera from each side are sent along with some of their species, as TimeTree knows many species but not every genus. If TimeTree knows none of them on one side, other genera and species are tried, and as a last resort enclosing taxa (family up to phylum) that contain no genus from the other side, so that the MRCA is still the same node. An age is only accepted if TimeTree found taxa from at least two sides.
- **Conflicts.** If a node's age is not younger than its nearest dated ancestor's: equal ages mean TimeTree does not resolve the node; otherwise a direct age outranks one obtained via enclosing taxa, and between equals the ancestor is kept.
- **Interpolation.** Undated nodes are spaced evenly between their nearest dated ancestor and their oldest dated descendant (as in BLADJ).
- **Genus crown ages** from TimeTree (the MRCA of up to twelve species of the genus), checked in step 3.

Each genus tip of the backbone has a stem age S: when it split from its sister group.

### 3. One dated subtree per genus

Genera are processed in parallel shards. For a genus with at least two units:

1. **Outgroup**: up to `n_outgroup` units from the nearest genera in the backbone that have units.
2. **Distances** between the units and the outgroup: K2P on a within-genus alignment (MAFFT) for COI, k-mer (Mash) distances without alignment for ITS.
3. **BioNJ** tree, **rooted** on the outgroup, which is then removed. If the ingroup is not a clade with respect to the outgroup, the ingroup is midpoint-rooted instead.
4. **Crown age C**, in order of preference:
   - the user's `crown_ages` table;
   - the TimeTree estimate, if it is younger than S and, for COI, at most `timetree_max_ratio` times the rate-based age;
   - the marker's method: for COI the mean root-to-tip distance divided by a substitution rate; for ITS a stated fraction of S.
   C is always capped at `cap_fraction` (0.9) of S.
5. **Dating**: the subtree is made ultrametric with its root at C, by mean path length (default) or LSD2.

A genus with one unit stays a tip, relabelled with that unit.

### 4. Grafting

Each genus tip is replaced by a stem of length S − C followed by its subtree, which keeps the tree ultrametric. Backbone genera without units are pruned. There is no tree search in this step.

### 5. Placement

Units without a genus, and metabarcoding ASVs, are compared with the reference records by VSEARCH global alignment (which handles fragments such as the Leray window of COI):

- **ASVs** at or above `match_identity` (0.97) count as the unit they match; equal best hits in several units attach at their MRCA.
- **Below that, down to `place_identity` (0.80)**, the sequence is attached on the path from its nearest unit to the root, at an age implied by its distance: (1 − identity) / (2 × rate) for COI, or halfway along the nearest unit's terminal edge for ITS.
- **Anything below** is reported as unplaced, with its reads.

### 6. Diversity

From the tree and a table of units (or placed ASVs) per sample: Faith's PD (including the path to the root), MPD (presence and abundance-weighted), unweighted and normalised weighted UniFrac, and PhyloSor. All metrics come from one edge-by-sample summary, so no tree is pruned per sample.

### 7. Validation

Pairs of units are sampled per taxonomic stratum (same genus, same family, same order, deeper) and each distance method is compared with the patristic distance in a trusted reference tree: rank correlation per stratum, binned medians, and the reference distance at which the barcode distance plateaus. This is how the informative window of a marker is measured. It needs a reference tree from independent data; a tree built from the same barcodes tests nothing.

## Assumptions

- **The backbone is right.** Topology from the OToL synthesis and ages from TimeTree are taken as given above genus level. Errors there propagate to PD, and PD is dominated by deep branches.
- **Genera are clades.** The pipeline builds one subtree per genus and grafts it onto one tip. Non-monophyletic genera break this (see weaknesses).
- **Barcodes are informative within a genus.** The default treats the genus as the informative window; the validation step exists to check this per clade.
- **A clock within each genus.** Mean path length and LSD2 both assume roughly clock-like evolution within a genus subtree.
- **A COI substitution rate** (default 0.0115 substitutions per site per million years per lineage, about 2.3% pairwise) when no crown age is known. It also sets the attachment age of placed ASVs and the threshold for rejecting TimeTree crown ages.
- **One representative per unit** stands in for the unit's variation.
- **One genetic code per run** for the COI reading-frame check.

## Known weaknesses

From the first real run ([results](results/2026-11-iberia.md)) and from the design:

- **TimeTree coverage of marine invertebrates is patchy.** Many polychaete, amphipod and some echinoderm genera are unknown to it. Undated nodes are interpolated, and even spacing tends to make shallow undated nodes too old (e.g. *Hediste*/*Nereis* at 291 Ma in a small test).
- **TimeTree crown ages are often inflated** for non-monophyletic genera: the MRCA of all species under a name can be hundreds of millions of years old (*Haliclona* 444 Ma). The ratio check rejects these, but it depends on the assumed COI rate.
- **OToL contains misplacements.** In the first run the gastropod *Onchidella* sat among the sea urchins, because OTT merged it with an echinoderm homonym. Conflicts between TimeTree and the backbone (`conflict_ratio` in `node_ages.tsv`) help find such cases, but nothing corrects them automatically.
- **Rooting with the sister genus often fails**: in about one genus in five, the outgroup fell inside the ingroup and midpoint rooting was used. Likely causes are non-monophyletic genera and saturation over the distance to the sister genus.
- **Distance-based trees and placement** are less accurate than likelihood-based ones (cf. EPA-ng against APPLES); the gain is scale and marker independence.
- **The single genetic code** is wrong for some groups in a mixed set (vertebrates use table 2, echinoderms and flatworms table 9, cnidarians table 4); with table 5, vertebrate sequences pass because its stop codons are a subset.
- **Placement thresholds** (`match_identity`, `place_identity`) are not calibrated; at 80% identity the rate-based attachment age is already about 9 Ma.
- **Units are not equivalent across markers.** BINs and UNITE species hypotheses are delimited differently, so PD across markers is comparable in units of time but not in the granularity of the tips.
