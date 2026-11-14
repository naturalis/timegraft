# Consistency layer: evidence on the Iberian run

> Each section below was measured on the state of the pipeline at its own issue; later issues changed some of the earlier numbers (for example, the screening fix of #55 changed 175 decisions). The numbers for the pipeline as merged into `main` are in `benchmark/` (#56), mainly `benchmark/results/end_to_end.tsv` and `benchmark/REPORT.md`.

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

## `requery/` (#53)

Conflicting nodes re-queried with disjoint picks and classified ([#53](https://github.com/naturalis/timegraft/issues/53)), at the defaults (`requery.mode: classify`, `z_trigger: 2`, `ratio_trigger: 1.5`, `move: 1.25`, `attempts: 3`), on top of `agefit/`.

- `node_ages.tsv` and `metadata.json` are from that run. The backbone itself is identical to `agefit/`, since classification changes no age.
- `node_ages_apply.tsv` is from the same run with `mode: apply`.
- The re-queries are in `timetree_cache.json`. Rebuilding from it asks TimeTree nothing and gives the same classes.

**22 nodes were re-queried**: 21 as one of the two nodes of a precedence conflict of 1.5 or more (among them the crab clade, which is also the fit's outlier, 3.6 sigma off), and one gadiform node for a fit residual of 2.4 sigma alone (unresolved).

| Class | Nodes | Which |
| --- | --- | --- |
| `under_spanning` | 1 (high) | The 27-genus crab clade: 54.9 Ma at first; with disjoint picks 156.8 Ma (6 estimates), exactly TimeTree's age for Eubrachyura below it. The first picks fell within Portunoidea, as the holistic analysis suspected |
| `topology` | 4 | Same age when asked again, conflict stands. Leuciscine cyprinids, 19.2 Ma against 29.6 from the fit (high). Hydrobiidae, 16.6 against 23.8 (medium). An eel node, 172.4 against 136.9 (high). Asellid isopods, 244.2 against 157.2 (medium) |
| `consistent` | 4 (high) | Partners of a conflict whose other side was at fault, among them Eubrachyura itself |
| `unresolved` | 13 | 9 where TimeTree found the disjoint picks on only one side (among them Lymnaeoidea and Hygrophila, and the heterotreme crab nodes), and 4 where every taxon of a side had been used (among them the sparid and the hydrobiid node) |

So the analysis's two readings are both seen. The crab clade is confirmed as under-spanning. Where TimeTree repeats its age with other picks, the conflict lies between OToL's and TimeTree's trees; these are the nodes to look at for topology. The 13 unresolved cases are what the method cannot reach with the taxa in the run.

**With `mode: apply`**, the crab clade's calibration becomes 156.8 Ma with its own uncertainty. The node is no longer an outlier, and its fitted age moves from 168.6 to 167.7 Ma: the fit had already set the first age aside. No other node moves by more than 1%. Against the calibrations as replaced, calibrations off by more than a factor 1.5 go from 5 to 4.

## `crowns/` (#54)

TimeTree crown ages of the species sampled, and the plausibility test against the barcodes ([#54](https://github.com/naturalis/timegraft/issues/54)), at the defaults (`crown_query.basis: sampled`, `min_sampled: 2`, `single_species: none`; `crown_age.plausibility` with `rate_factor: 10`, `saturation: 0.15`, `young_factor: 3`). The backbone is identical to `agefit/`.

- `crown_ages.tsv` and `metadata.json` are the backbone outputs, with both the sampled and the global crown per genus.
- `crowns_compare.tsv` holds the crown age of each of the 1,039 genera with a subtree, under the four combinations of basis and test.
- `crown_compare.py` makes that table from a genus-stage run. The crown decision depends only on the root-to-tip distance, rate and stem age the run records, together with the crown age table, so the trees need not be rebuilt per combination.
- The genus stage was run as the Iberian configuration does it, with the literature crown ages of `data/crown_ages/literature.tsv` (*Gammarus*, *Niphargus*, *Daphnia*), on shards made from this backbone and `units/screened.tsv` of the archive. All four combinations apply that table first, as the workflow does. Recomputed with the run's own settings, the decisions match the run for all 1,039 genera.
- The new TimeTree queries are in `timetree_cache.json`.

**Basis.** 312 genera get the crown of their sampled species and 21 fall back to the global crown. 262 are sampled with a single species and get no TimeTree crown; 10 of those had a global one before. Of the genera with both answers, the sampled crown differs from the global one for 43, younger for 26 (*Trachurus* 24 to 6 Ma, *Mugil* 57 to 25, *Margaritifera* 230 to 109). Note that "sampled" means the species of the reference records, which for these BOLD data is mostly the genus worldwide (*Gammarus*: 184 species), so the basis changes less than it would for a truly regional reference.

**Extreme crowns**, among the 1,039 genera:

| | TimeTree crown used | Crown > 10 × barcode age | > 20 × | Crown < barcode age / 3 |
| --- | --- | --- | --- | --- |
| Before #54 (global, no test) | 267 | 42 | 7 | 5 |
| Sampled only | 281 | 41 | 10 | 6 |
| Global with the test | 251 | 31 | 0 | 0 |
| Sampled with the test (new default) | 260 | 26 | 0 | 0 |

The sampled crowns on their own let a few more extreme values in: more genera get a TimeTree crown, and TimeTree sometimes knows only a few of the sampled species. *Gammarus* is the clearest case: TimeTree knows 4 of its 184 sampled species, all close, and gives 1.1 Ma against 13.8 Ma from the barcodes. Here the literature table sets its crown instead. With the test, the extremes are gone in both directions.

**Demoted, too old** (15; TimeTree crown of the sampled species, then the new crown age): *Ruditapes* 69 to 9.0 Ma (1.3% divergence), *Haliclona* 444 to 73, *Mycale* 162 to 29, *Hydra* 193 to 40, *Phorbas* 127 to 25, *Tethya* 131 to 27, *Aulopus* 53 to 11, *Cliona*, *Ircinia*, *Spongia*, *Holothuria*, *Margaritifera*, *Mimachlamys*, *Chimaera* and *Illex*. These include all seven genera the holistic analysis found above 20 times the barcode minimum.

**Not used, too young** (6): *Stephanocyathus*, *Conchoderma* (1.0 Ma against 9.3), *Gadomus* (0.1 against 3.6), *Leptaxis*, *Euchaeta* and *Orculella*. Their crowns become the midpoint between the barcode age and the stem cap.

**Overall.** 76 crown ages change, 26 of them older, by a median factor of 1.85. The sum of crown ages falls by 4.8%. The test allows a genus to evolve ten times slower than the rate, so 26 genera remain above 10 times the barcode age. Most are sponges, cnidarians and elasmobranchs, where that is plausible.

## `screening/` (#55)

Screening of the reference units with the extensions of [#55](https://github.com/naturalis/timegraft/issues/55): a wider reference, thresholds per class, a codon comparison with the nearest relative, a choice for mutual pairs, a relative-rate check in the genus subtrees, and a reason for every decision. `benchmark.py` makes every table here except the two `rate_outliers` tables, from `units/representatives.tsv` of `iberia-run-tolerant.tar.gz` and the records of `records.tsv.gz` without those in its `units/rejected.tsv`, with VSEARCH 2.29.1 on two cores. Run on 13 November 2026.

- `screened.tsv.gz` and `calibration.tsv`: `units/screened.tsv` and `units/screening_calibration.tsv` as `screen_units` writes them at the defaults.
- `injection.tsv`, `injection_classes.tsv`: mislabelled units, over a grid of thresholds and per class.
- `numt_simulation.tsv`: simulated nuclear copies against real pairs of nearest relatives.
- `external.tsv`: the wider reference.
- `runtime.tsv`: time per step.
- `rate_outliers.tsv`: the new provenance columns from a genus-stage run, made as for `crowns/` (same backbone, crown table and literature ages, `units/screened.tsv` of the archive). `rate_outliers_set_aside.tsv`: the same with `rate_outlier_action: set_aside`.

**What changes at the defaults.** The rules are those of #38, but the search now fills in what it used to miss. VSEARCH keeps the first ten hits of each unit, and for 12,481 units all ten were of one side: 11,734 only congeners, 747 only units of other genera. Those units are searched again with up to 100 hits (`deep_maxaccepts`). That gives 4,949 units a nearest unit of another genus and 276 a nearest congener that they lacked, and changes 175 decisions: 163 units that were `other_genus` only because their congeners had been crowded out are now kept, and 12 others change in smaller ways (6 kept units become `other_genus`). The defaults now set aside 1,139 units (695 `other_genus`, 444 `isolated`; 4.6%), against 1,294 in the archive; 23,687 are kept and 56 kept with their whole genus. Every row of `screened.tsv.gz` has the thresholds that applied, where its nearest hits came from and its reason in words, so no unit leaves its genus without a stated cause.

**Mislabelled units.** In each of five replicates, 2% of the kept units of genera with at least three units were given the name of another such genus, either of the same family (333 units) or of another order of the same phylum (467). The search is not repeated: only the names change, so each setting is cheap to test.

| `margin` | `mutual_pairs` | Caught, same family | Caught, other order | Correct units set aside with them, per 333 (family) | Set aside without injected errors |
| --- | --- | --- | --- | --- | --- |
| 0 | both | 0.93 | 0.98 | 349 | 7.9% |
| 0.01 | both | 0.91 | 0.98 | 163 | 5.8% |
| **0.02** | **both** | **0.89** | **0.99** | **112** | **4.6%** |
| 0.03 | both | 0.86 | 0.99 | 90 | 4.1% |
| 0.05 | both | 0.79 | 0.99 | 59 | 3.7% |
| 0.02 | further | 0.88 | 0.99 | 25 | 4.2% |

All at `min_identity` 0.80. Errors across orders are nearly always caught; errors within a family are the hard case, and one in nine is missed at the defaults (36 of 333 in one replicate): 29 units whose congener and other-genus identities lie within the margin of each other, and 7 in genera that would otherwise lose all their units. Raising `min_identity` catches a few more within families (0.93 at 0.85) but sets aside 15% of units without injected errors, against 4.6%.

What the error does to correct units was the surprise. 76 of 86 correct units set aside in one replicate had an injected unit as their nearest unit of "another genus", closer than their own congeners: in most cases its nearest true congener. A pair cannot tell which of the two carries the wrong name. `mutual_pairs: further` keeps the unit of the pair that is closer to its own genus. It cuts the collateral by more than three quarters at the cost of one error in a hundred, and sets aside 91 fewer units in the run itself (604 `other_genus` instead of 695). It is not the default yet: what it does to rooting and crown ages is for #56 to judge.

**Thresholds per class.** `calibration.tsv` suggests `min_identity` from 0.80 (most classes) to 0.89 (Demospongiae), from the 1% quantile of congener identity. With the suggestions, recall per class is the same (Teleostei 0.88 to 0.89; the other classes unchanged), and slightly more units are set aside (Gastropoda 4.3% to 5.1%, Branchiopoda 1.7% to 2.5%), so the defaults stay. The classes differ more in how many units are set aside without injected errors than in recall: from 1.1% (Elasmobranchii) and 1.9% (Teleostei) to 9.4% (Polychaeta), 10% (Demospongiae) and 15% (Gymnolaemata). Those high shares are where genera are least natural or least well identified, not where the thresholds are wrong. Recall per class rests on few injected units in some classes (4 sponges, 8 bryozoans per replicate set).

**Nuclear copies.** 2,000 unit sequences were copied with random substitutions at all codon positions, as a nuclear copy accumulates them, and compared with their source. At the defaults (`min_changed_codons: 8`, `max_nonsynonymous_share: 0.5`), copies 2% from their source are flagged `numt_like` in 97% of cases and copies 4% or more from it in nearly all (99.8% and more); copies 1% from it are not, because fewer than eight codons differ. Such copies fall within the BIN of their source and do not form a unit of their own. Among real pairs of nearest relatives, 1.1% are flagged: 265 units, at a median identity of 0.97 to their relative with 17 changed codons. Nuclear copies and sequences with many errors look alike in this test. Of the copies with substitutions, between 86% (1%) and 34% (8%) still pass the reading-frame filter, so the filter alone misses most recent copies.

A first version took the reading frame of the frame with fewest stops and, among equals, the first. Under the echinoderm code, with only two stops, several frames are often open, and eight sea cucumber units (seven *Holothuria*, one *Psolus*) were flagged for changes that are silent in their true frame. Among equally open frames, the one with the fewest amino-acid changes is now taken, and all eight have nearly all their changes silent. Two other *Holothuria* units remain flagged: a pair at 97% identity of which 15 of 21 changed codons alter the amino acid.

245 pairs are out of frame over a stretch (`frameshift`), mostly by a single base in either sequence; the pair cannot tell which. They are reported, not acted on.

**Relative rates.** On the genus stage of `crowns/`, which used the screening of the archive: in 820 genera with at least four units, 412 units in 235 genera lie more than 3.5 robust standard deviations above the median root-to-tip distance of their genus, 116 genera above 5 and 33 above 10. The largest scores come from small genera whose units are otherwise very close (*Malacolimax*, *Baeolidia*, *Hypleurochilus*). 46 of the 412 units (11%) are `numt_like`, against 1.1% of all units: the two diagnostics point at the same sequences far more often than by chance. By default the list is a diagnostic only. With `rate_outlier_action: set_aside` (genus stage rerun as above), the 412 units leave their subtrees and are placed like other set-aside units: mean root-to-tip distance in the 235 genera falls by a median of 4%, and 173 crown ages change, by a median of 2% (8 older, 165 younger) and by more than a quarter in 10 (*Pelagia* 47 to 18 Ma, *Bathytoshia* 16 to 8 Ma, *Myzobdella* 32 to 19 Ma). The other 62 genera take their crown from TimeTree or the literature. The sum of all crown ages falls by 1.1%. Whether that is better is for #56 to judge; `rate_outliers_set_aside.tsv` has the provenance columns of that run.

**A wider reference.** The genera of each family were split at random between the run (16,011 units) and an external reference (8,870 sequences, all within families of the run). Screening the run half alone sets aside 683 units, of the 783 that the screening with all genera sets aside among them: 108 are missed, because the genus they are closest to is not in the run, and 8 are extra. With the reference, 780 of the 783 are found and none is extra; 173 decisions rest on an external hit, as `best_other_scope` records.

**Runtime**, on two cores: the all-against-all search of 24,882 units takes 95 s with ten hits per unit. Searching the 12,481 one-sided units again with up to 100 hits, and keeping the alignments for the codon comparison, brings that to 269 s, nearly three times as long; `deep_maxaccepts` below 100 saves time at the risk of leaving a side empty. Decisions take 0.9 s, the codon comparison 9 s and the calibration 0.1 s. The external reference cost nothing measurable: 185 s for the 16,011 units of the run half, 182 s with 8,870 more sequences in the reference, because the search time is set by the queries more than by the reference. The relative-rate check only summarises root-to-tip distances that the genus stage already computes; that stage took 26 minutes here and 24 for `crowns/`, a difference that time spent on these sums cannot account for. Memory is not an issue: the kept alignments are two per unit.

## `benchmark/` (#56)

The release gate of the whole layer ([#56](https://github.com/naturalis/timegraft/issues/56)): `consistency-layer` against `main`, both run end to end, plus backbone builds with injected misplacements, perturbed TimeTree ages and half the genera dropped. `benchmark/GATE.md` has the metrics and thresholds as declared before the run, `benchmark/REPORT.md` the outcome and the recommendation, `benchmark/README.md` how to rerun it, and `benchmark/results/` the tables.

In short: 25 of 28 must metrics pass. The real misplacement is fixed end to end; 118 of 120 injected merges are found and 117 placed back in their class, with no false flag; conflicts and extreme crown ages fall; the backbone moves less under noise in TimeTree (node ages that move by more than a factor 1.5: 21% to 9%) and when half the genera are dropped; sample rankings by PD are stable (0.97 to 0.99); the run takes 3% to 9% longer and is reproducible offline. The three failures: after an injected merge is corrected, the species of the wrong taxon still go to TimeTree and the genus returns to its family, not to its place in it (two metrics), and peak memory is 1.61 times that of `main` at 1.3 GB (one metric, with a threshold set too tight). The report recommends merging with the defaults on and lists the follow-ups.

