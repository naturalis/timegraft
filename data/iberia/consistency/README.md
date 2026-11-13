# Consistency layer: evidence on the Iberian run

Material for the consistency layer ([#46](https://github.com/naturalis/timegraft/issues/46)) on the Iberian run. Each issue of the epic adds what it needs here, so that the benchmark of [#56](https://github.com/naturalis/timegraft/issues/56) can compare every step with the same starting point.

## `backbone/` (#47)

The backbone of `iberia-run-tolerant.tar.gz` rebuilt with the evidence columns of [#47](https://github.com/naturalis/timegraft/issues/47) (described in `docs/technical/provenance.md`): `genera.tsv`, `node_ages.tsv`, `crown_ages.tsv`, `metadata.json` and `backbone.tre`.

How it was made, on 12 November 2026: `workflow/scripts/build_backbone.py` was run on `units/representatives.tsv` of the archive, the records of `records.tsv.gz` without those in its `units/rejected.tsv` (the workflow's `units/filtered.tsv`), OTT 3.7.3, synthesis 16.1 with its `annotations.json`, and `timetree_cache.json` as it was, with `upgrade_cache: true`. The run upgraded all 2,507 cache entries to version 2; TimeTree (database version 20260820) gave the same age and found the same taxa for every one of them. The upgraded cache is now `../timetree_cache.json`.

Checks that nothing changed:

- The same script at `main` (commit `f4f3d7f`), on the same inputs, reproduces the backbone files of the archive byte for byte.
- Against that, `backbone.tre` is identical, and the earlier columns of `node_ages.tsv` (779 nodes), `genera.tsv` (1,218 genera) and `crown_ages.tsv` (312 genera) are identical in content and order; only columns were added.

The new columns reproduce the counts of the holistic analysis: of the 314 TimeTree node ages used, 73 rest on one estimate, 89 on two or three and 152 on four or more; 218 of 779 backbone nodes have no source tree in the synthesis; 38 genera carry more than one NCBI id in OTT; 20 genera attached at the MRCA of their species sit on a clade of more than 50 genera. For *Onchidella*, `genera.tsv` shows OTT's lineage (Laganidae, Echinoidea) against that of all 82 of its records (Onchidiidae, Gastropoda, Mollusca), and `node_ages.tsv` shows that the two dropped ages above it were TimeTree's age of NCBI 33213, Bilateria.

## `backbone-col/` (#48)

The same backbone with the Catalogue of Life as identifier hub ([#48](https://github.com/naturalis/timegraft/issues/48)), at the default settings: release 316691 (COL26.10 XR), `col.join: col`, no species. Made on 12 November 2026 like `backbone/`, with `timetree_cache.json`, which now also holds the ChecklistBank answers: rebuilding from it needs no network access and gives the same files.

What COL resolved, for the 1,218 genera:

- 1,123 genera at `high` confidence, 87 `medium`, 8 `low` (taxa that come from XR only, such as woodlice and bathynellaceans); none left unresolved. All matched with the lineage of their records as context.
- ChecklistBank calls 57 of these names ambiguous, 47 of them only because synonyms of the name exist. For the other 10, several accepted genera of the name exist and the record lineage chose between them (*Mitra*: Mitridae against a ciliate; *Lysidice*: Eunicidae against Eulepethidae and a legume).
- The NCBI id derived from COL is the one OTT uses for 1,126 genera; for 21 it is another of the ids OTT merged into the genus; for one (*Cuvierina*) OTT holds a different id; for 70 one side has none.

What that changed:

- **Three genera added to the backbone.** OTT has two genera of these names in Metazoa, so they were left `ambiguous` before. COL's OTT id chose between them: *Cepola* (Cepolidae), *Micropterus* (Centrarchidae) and *Sergia* (Sergestidae). With `col.join: ott`, 36 of the 779 earlier nodes change as a result, in the genera they hold or in age.
- **48 genera sent to TimeTree with COL's NCBI id.** For 21 it was another of OTT's merged ids. These include *Onchidella*: 36932, the gastropod, instead of 262330, the sand dollar *Peronella*. For 27, OTT had no NCBI id. This changes 30 node ages, and no crown ages. The largest changes are among bryozoans, which gained NCBI ids: *Caberea* + *Notoplites* is now dated by TimeTree at 108.8 Ma instead of interpolated at 214.2 Ma, and the nodes above it move by up to 90 Ma. The others are a few crustacean nodes: one moves from 482 to 601 Ma and another from 465 to 535 Ma. The node joining *Onchidella* to the echinoderms is no longer dated through higher taxa at 515 Ma and is left to interpolation (495 Ma): with the right id, TimeTree is now asked about a gastropod, while the genus still sits among the sea urchins in OTT (#49, #50).

Whether these moves are improvements is for the benchmark of #56 to judge; `col.join: ott` gives the identifiers as before.

## `concordance/` (#49)

Placement concordance ([#49](https://github.com/naturalis/timegraft/issues/49)) on the backbone of `backbone-col/`, at the default settings (`flag_level: class`, `max_share: 0`, `min_hits: 5`, `min_flags: 2`): `concordance.tsv`, and `genera.tsv` and `metadata.json` with the concordance columns filled. The backbone itself is unchanged: `backbone.tre`, `node_ages.tsv`, `crown_ages.tsv` and every other column of `genera.tsv` are identical to those in `backbone-col/`. The barcode hits are those of `units/screened.tsv` in `iberia-run-tolerant.tar.gz`; the NCBI lineages of merged ids were read from ChecklistBank and are now in `timetree_cache.json`.

**Result.** One genus flagged, *Onchidella*, which fails all three checks. Two genera for review: *Arbacia* and *Magelona*. 34 genera cannot be checked, mostly because they are not in OTT. The other 1,181 pass.

**Each check on its own**, and why the genera it fails on are not misplacements, apart from *Onchidella*:

| Check | Evaluable | Fails at default | Distribution and false positives |
| --- | --- | --- | --- |
| Reference lineage (COL, else records) | 1,184 | *Onchidella* (Onchidiidae against Laganidae: they meet only at Bilateria) | 1,115 agree, 59 part at family level, 2 at superfamily level and 7 at order level. All 68 are classification differences: COL splits Gobionidae and Tincidae from Cyprinidae, places Andvakiidae differently, and uses other sponge orders (*Achaeta*, *Antigonia*, *Gobio*, *Phakellia*, *Paratimea*, *Telmatactis*, *Tinca*). With the records' own lineage as the reference instead of COL, the picture is the same (63 at family level, 6 at order level), plus *Trachysalambria*, a penaeid shrimp that BOLD files under Pilumnidae. With COL as the reference that BOLD error disappears |
| Merged ids | 35 of the 38 genera with more than one NCBI id | *Onchidella* (36932, a gastropod, and 262330, the sand dollar *Peronella*) | 29 of the merged pairs are the same taxon, 3 part at family level, *Patella* at superfamily level and *Myosotella* at order level. These are synonyms and lumped genera. 3 ids have no lineage in ChecklistBank |
| Barcode neighbours | 676 genera with at least 5 cross-genus hits | *Onchidella* (16 hits, all pulmonates; sisters are echinoderms), *Arbacia* (its sister in the backbone is *Onchidella*), *Magelona* (11 hits in molluscs; its COI looks molluscan, probably a marker artefact) | *Plumularia* has the lowest share above zero, 0.25: 24 of its hits are the brittle star *Ophiothrix*, the epibiont contamination noted in the analysis. Its other 8 hits are hydroids, so it passes |

The false positives of the three checks do not overlap. That is why two failures are required, as the holistic analysis proposed.

**Sensitivity.** *Onchidella* is still the only flagged genus whether `flag_level` is `superfamily`, `order`, `class` or `phylum`, and whether `max_share` is 0 or 0.25. Only the number of genera for review changes: 13 or 14 at superfamily level, 10 or 11 at order level, and 2 or 3 at class or phylum level. The run has a single real misplacement, so it cannot tell how many misplacements the checks miss. The error injection of #56 is meant for that.

## `reanchor/` (#50)

The backbone with flagged genera re-anchored before dating ([#50](https://github.com/naturalis/timegraft/issues/50)), at the defaults (`mode: apply`, `min_severity: class`, no curated table). It contains `backbone.tre`, `node_ages.tsv`, `crown_ages.tsv`, `genera.tsv`, `concordance.tsv` and `metadata.json`. The new TimeTree queries are in `timetree_cache.json`.

**What moved.** Only *Onchidella*, the one flagged genus. Its reference family, Onchidiidae (from COL), is a node of the synthesis (ott474706), and *Onchidella* now hangs there instead of at its own synthesis node among the sea urchins. Its row in `genera.tsv` records the trigger (three of three checks failed, severity 5), the original anchor (ott162893), the target and how it was found (`reference_lineage;taxon_in_synthesis`), and `auto` as the source. No other genus is affected; `mode: report` and a curated `keep` for *Onchidella* both reproduce the backbone of `backbone-col/` exactly.

**What that changed:**

- **Sisters and outgroups.**
  - *Onchidella*'s nearest sister clades are now 325 gastropod genera, the bivalves and the cephalopods. Before, they were *Arbacia*, the other sea urchins and the sea cucumbers.
  - *Arbacia*'s sisters are now *Paracentrotus*, *Psammechinus* and *Sphaerechinus*. Before, its only sister was *Onchidella*.
- **Barcode neighbours.** The check rerun on the corrected backbone (`neighbour_concordance_after`) goes from 0 to 1 for both *Onchidella* and *Arbacia*.
- **Stem ages.**
  - *Arbacia*: 360.0 Ma before, now 207.5 Ma.
  - *Onchidella*: 360.0 Ma before, now 466.4 Ma. The synthesis resolves Onchidiidae only to a 30-way polytomy within the gastropods, which TimeTree dates at 466 Ma (19 estimates). The genus is in the right clade now, but at the depth the synthesis allows, which is a matter for #51 and #52.
- **The two deepest nodes below the root.** Before, they were 750.0 and 749.9 Ma; now they are 682.3 and 682.0 Ma. With *Onchidella* among the sea urchins, TimeTree had dated their node at 721 Ma (Bilateria). That age outranked the higher-taxon age of 682 Ma for the node above, so both were dropped and interpolated near the root. One misplaced genus therefore shifted the deepest backbone ages, which enter the PD of every sample. This is the propagation the holistic analysis describes.
- **Everything else.** 767 of the 780 node ages are unchanged.

## `anchoring/` (#51)

The backbone with majority-clade anchoring ([#51](https://github.com/naturalis/timegraft/issues/51)) at the defaults (`rule: majority`, `share: 0.5`, `min_species: 3`, `prefer_sampled: true`, `min_sampled: 2`), on top of `reanchor/`: `backbone.tre`, `node_ages.tsv`, `crown_ages.tsv`, `genera.tsv`, `concordance.tsv`, `metadata.json`. `rule: mrca` with `prefer_sampled: false` reproduces `reanchor/` exactly.

**Which genera.** The rule applies to the 149 genera that are not a node of the synthesis and are anchored on their species; the 998 that are nodes and the 38 attached at an enclosing taxon are untouched (998: the 996 of the original run, plus *Cepola*, *Micropterus* and *Sergia*, which COL resolved in #48, minus *Onchidella*, which has status `reanchored` since #50). For 145 of the 149 the anchor uses the species sampled in the run; for 4 there were fewer than two of those in the synthesis, so all species of the genus in OTT were used. 62 genera are now anchored on a majority clade. 75 keep the MRCA because no smaller clade holds a majority of their species (`no_smaller_majority`), and 12 keep it because they have fewer than three species in the synthesis (`too_few_species`).

**Outlier sensitivity.** Genera anchored at the base of a clade of more than 50 genera go from 20 to 12, and of more than 10 genera from 61 to 43. Examples, as genera in the anchor clade under the MRCA and then the majority rule:

- *Carassius*: 277 to 2, with *C. gibelio* and *C. cuvieri* outside.
- *Caryophyllia*: 269 to 10.
- *Oculina*: 235 to 3.
- *Squatina*: 117 to 1.
- *Barbus*: 114 to 2.
- *Octopus*: 49 to 10.

`anchor_outside_share` and `anchor_outside_species` give, for each genus, the species left out.

**Effect on ages.**

- 31 stem ages change: 25 decrease and 6 increase.
- Largest decreases: *Etmopterus* 258 to 130 Ma, *Loligo* and *Alloteuthis* 158 to 77 Ma, *Squatina* 258 to 182 Ma, *Pagellus* and *Dentex* 62 to 22 Ma, *Balanus* 94 to 57 Ma, *Raja* and *Dipturus* 98 to 68 Ma.
- The increases are small, at most 13 to 20 Ma (*Sepietta*, *Rondeletiola*).
- 758 of the 780 earlier node ages are unchanged; the backbone now has 784 nodes, because genera that were part of polytomies are separated.
- The total branch length of the backbone falls by 0.25%.
- Placement concordance is unchanged: *Onchidella* is flagged and re-anchored, and *Arbacia* and *Magelona* are for review.

Taking each part on its own:

- Sampled species with the MRCA rule changes one anchor, for *Loligo*.
- The majority rule on all species moves 16 genera and changes 8 stems.

**What it does not fix.** The worst cases in the analysis keep their anchors, because their species are scattered over the synthetic tree rather than drawn out by a few strays: *Mycale* (465 genera, all 8 sampled species spread out), *Munidopsis* (238), *Trachurus* and *Pseudocaranx* (190), *Periclimenes* (170) and *Phakellia* (136). They are marked `no_smaller_majority`. They need either weighting of the conflicting TimeTree ages around them (#52) or a curated decision.

## `agefit/` (#52)

Backbone ages from a weighted fit of all TimeTree calibrations ([#52](https://github.com/naturalis/timegraft/issues/52)), at the defaults, on top of `anchoring/`. The folder holds `backbone.tre`, `node_ages.tsv` and `metadata.json`. `genera.tsv`, `concordance.tsv` and `crown_ages.tsv` are identical to those in `anchoring/`, because the topology does not change. `ages.method: precedence` reproduces `anchoring/` exactly.

**Conflict metrics over the 387 TimeTree calibrations** (from `metadata.json`; a factor is the final age over TimeTree's, or the reverse):

| | Median factor | 90% | Largest | Off by more than 1.5 | RMS of residual / sigma |
| --- | --- | --- | --- | --- | --- |
| Precedence rule | 1.00 | 1.32 | 6.49 | 26 | 1.32 |
| Weighted fit | 1.00 | 1.11 | 3.07 | 5 | 0.38 |

This is what the fit optimises, so it shows that the fit honours the evidence it is given. It does not show that the ages are better.

**Held-out TimeTree ages** (`holdout.py`, 50 random halves of the calibrations hidden; each hidden age is compared with the age the visible half gives after interpolation):

| | Median factor | 90% | Off by more than 1.5 | Off by more than 2 |
| --- | --- | --- | --- | --- |
| Precedence, all hidden | 1.319 | 2.52 | 35.6% | 17.3% |
| Fit, all hidden | 1.318 | 2.49 | 35.0% | 16.9% |
| Precedence, 4+ estimates | 1.210 | 1.88 | 22.5% | 7.2% |
| Fit, 4+ estimates | 1.212 | 1.85 | 21.9% | 6.7% |

Where a node has no calibration of its own, the fit predicts about as well as the precedence rule, and a little better in the tail. Its gain is at the calibrated nodes: it no longer throws away well-supported ages.

**The conflicts listed in the analysis** (TimeTree age and number of estimates; precedence age; fitted age):

| Node | TimeTree | Precedence | Fit |
| --- | --- | --- | --- |
| 27-genus crab clade (ancestor) | 54.9 Ma (4) | 54.9 | 168.6, reported as the one outlier (3.6 sigma) |
| Eubrachyura | 156.8 Ma (6) | 54.8 (dropped) | 166.9 |
| Heterotreme nodes | 195.5 Ma (2) | 30 to 44 (dropped) | 164 to 167 |
| Lymnaeoidea | 136.8 Ma (2) under Hygrophila 50.0 Ma (1) | 39 to 44 (dropped) | 106.8 |
| Soleid node | 63.6 Ma (2) under Soleidae 14.9 Ma (2) | 17.2 (dropped) | 44.3 |
| Sparid node | 62.2 Ma (1) under Sparidae 29.2 Ma (1) | 25.8 (dropped) | 42.2 |
| Hydrobiid node | 33.9 Ma (1) under Hydrobiidae 16.6 Ma (1) | 8.3 (dropped) | 23.6 |

For the crab clade, the fit comes out as the analysis argued: it trusts the descendants. Where both sides rest on one or two estimates, the fit settles between them, which is as much as the evidence allows. Re-querying with disjoint picks (#53) is meant to settle those cases.

**Overall:**

- 118 of the 784 node ages move by more than 10% and 40 by more than a factor 1.5, almost all upwards.
- 136 genus stems move by more than 10%: 130 older and 6 younger.
- The total backbone length grows by 2.8%. The root stays at 750.1 Ma.
- TimeTree crown ages, which are capped at 0.9 of the stem, will therefore be capped less often in the genus trees.
