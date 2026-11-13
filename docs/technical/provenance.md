# Evidence and provenance of the backbone

The backbone outputs record not only the ages and placements timegraft uses, but also the evidence behind them: what TimeTree said, how each genus was matched in OTT and anchored in the synthetic tree, and what lineage its reference records give it. This is the first step of the consistency layer ([#46](https://github.com/naturalis/timegraft/issues/46), [#47](https://github.com/naturalis/timegraft/issues/47)): later steps score placements and ages against each other, and need these columns to do so.

None of these columns changes a result. With the same inputs and the same cache, the backbone tree, every node age and its source, the genus matches and the TimeTree crown ages are identical to those of earlier versions; only columns have been added.

## TimeTree cache, version 2

TimeTree answers are cached in `backbone_build.cache` under keys `mrca:<id>+<id>+...` (NCBI ids, sorted).

| Version | Entry | Written by |
| --- | --- | --- |
| 1 | `{"found": [...], "age": ...}` | timegraft before #47 |
| 2 | `{"v": 2, "found": [...], "age": ..., ...}` with the fields below | timegraft from #47 |

`found` and `age` mean the same in both versions, and the dating rules read only those two. Version 2 adds:

| Field | From TimeTree | Meaning |
| --- | --- | --- |
| `missing` | `missing_ids` | Queried ids TimeTree does not know |
| `mrca_ncbi` | `taxon_id` | NCBI id of the node TimeTree dated; compare with the OTT node to check that TimeTree dated the intended split |
| `mrca_ttid` | `mrca_ttid` | TimeTree's own id of that node |
| `estimates` | `time_estimates` | The individual time estimates, in Ma |
| `ci_low`, `ci_high` | `precomputed_ci_low`, `precomputed_ci_high` | Confidence interval of `age`; empty where TimeTree gives none (a single estimate) |
| `adjusted_age` | `adjusted_age` | TimeTree's adjusted age, where it gives one |
| `studies` | `study_data` | Per study: first author, year, number of trees and its age estimate (titles are left out) |
| `error` | `error` | TimeTree's message when it has no age, or the HTTP status of a failed query |
| `api_version` | `version` | Version of the TimeTree database that answered |
| `upgraded_from`, `fresh_age`, `fresh_found` | | See below |

**Upgrading a version 1 cache.** With `backbone_build.upgrade_cache: true` (the default), every version 1 entry that a run uses is asked again, once, for its evidence. The cached `found` and `age` are kept, so the backbone does not change with the cache format; `upgraded_from` is set to 1, and where TimeTree now answers differently the new answer is kept beside the old one as `fresh_age` or `fresh_found`. If TimeTree cannot be reached, the entry stays at version 1 and is tried again next run. With `upgrade_cache: false` nothing is asked again and the evidence columns of those queries say `legacy_cache`.

To build a backbone with TimeTree's current answers instead, start from an empty cache.

Version 1 code reads a version 2 cache without trouble, since it only reads `found` and `age`; a cache can therefore be shared between versions.

## `backbone/node_ages.tsv`

The first ten columns are as before. Added:

| Column | Meaning |
| --- | --- |
| `timetree_attempt` | Which set of picks gave TimeTree's age: `0` for the first, `1` to `3` for later attempts, `higher_taxa` for enclosing taxa; empty if TimeTree gave none |
| `conflict_with` | For a node in conflict (`conflict_ratio`), the node it conflicted with: the dated ancestor, or for a dropped `timetree_higher_taxa` age, the descendant that outranked it. Both describe the conflict as the precedence rule sees it, whichever `ages.method` is used |
| `timetree_evidence` | `full` (version 2 cache entry), `legacy_cache` (version 1, only `timetree_found_ids` known), `none`; empty if TimeTree gave no age |
| `timetree_found_ids` | The queried NCBI ids TimeTree found |
| `timetree_mrca_ncbi` | NCBI id of the node TimeTree dated |
| `timetree_n_estimates`, `timetree_n_studies` | Number of time estimates and of studies behind the age |
| `timetree_estimate_min`, `timetree_estimate_max` | Range of the estimates; their ratio is the spread used in the [holistic analysis](https://github.com/naturalis/timegraft/issues/46) |
| `timetree_ci_low`, `timetree_ci_high` | TimeTree's confidence interval |
| `timetree_api_version` | Version of the TimeTree database that answered |
| `synth_node` | Label of the node of the synthetic tree behind the backbone node (`ott…` or `mrcaott…ott…`); static source only |
| `synth_supported_by`, `synth_conflicts_with`, `synth_resolves`, `synth_partial_path_of`, `synth_terminal` | Number of source trees of the synthesis per kind of annotation, from its `annotations.json`. `synth_supported_by` of 0 means the node rests on the taxonomy alone. Empty without `backbone_build.synthesis_support` |
| `calibration_sigma`, `calibration_weight`, `fitted_age`, `calibration_residual`, `calibration_z`, `calibration_outlier`, `precedence_age`, `precedence_source` | The weighted fit of backbone ages; see [below](#backbone-ages-from-a-weighted-fit) |
| `conflict_class`, `requery_age` | Reserved for the re-query of conflicting nodes ([#53](https://github.com/naturalis/timegraft/issues/53)); empty |

The TimeTree columns describe the query in `query_ids`, the one whose age was used (or dropped, for `dropped_inconsistent` nodes).

## `backbone/genera.tsv`

The first four columns (`genus`, `status`, `ott_id`, `ncbi_id`) are as before. Added, for the static source (with `source: api` the OTT and anchor columns stay empty):

| Column | Meaning |
| --- | --- |
| `ott_resolution` | `col` where OTT had several genera of this name in context and COL chose between them ([#48](https://github.com/naturalis/timegraft/issues/48)); otherwise empty |
| `ott_candidates` | OTT genera of this name, anywhere |
| `ott_candidates_in_context` | Of these, those within `backbone_build.context` |
| `ott_candidate_ids` | For an `ambiguous` genus, the OTT ids it could not choose between |
| `ott_flags`, `ott_sourceinfo` | Flags and source ids of the matched OTT taxon, as in `taxonomy.tsv` |
| `ott_ncbi_ids`, `ott_worms_ids`, `ott_gbif_ids` | All ids of that source in `ott_sourceinfo`. More than one NCBI id means OTT merged records: mostly synonyms, but sometimes a wrong merge, as for *Onchidella* ([#20](https://github.com/naturalis/timegraft/issues/20)). `ncbi_id`, used for TimeTree, is the first |
| `ott_family`, `ott_order`, `ott_class`, `ott_phylum` | Nearest enclosing OTT taxon of each rank |
| `anchor_node` | Node of the synthetic tree the genus is a tip at (see `status`) |
| `anchor_clade_genera` | Genera in the clade of that node; large for a genus attached at the MRCA of its species when one species lies far away ([#51](https://github.com/naturalis/timegraft/issues/51)) |
| `n_species_ott`, `n_species_in_synthesis` | Species of the genus in OTT, and those of them in the synthetic tree |
| `anchor_rule`, `anchor_fallback`, `anchor_basis`, `anchor_species`, `anchor_species_inside`, `anchor_outside_share`, `anchor_outside_species`, `anchor_mrca_node`, `anchor_mrca_clade_genera` | How a genus that is not a node of the synthesis was anchored on its species; see [below](#anchoring-genera-that-are-not-a-node-of-the-synthesis) |
| `records_n` | Reference records of the genus |
| `records_family`, `records_order`, `records_class`, `records_phylum` and their `_share` | The most common name per rank among those records, and the share of records naming that rank that agree. This is the reference lineage that placement concordance compares with OTT ([#49](https://github.com/naturalis/timegraft/issues/49)) |
| `col_*`, `join_ncbi_id`, `join_ncbi_source` | Catalogue of Life as identifier hub; see [below](#catalogue-of-life-as-identifier-hub) |
| `lineage_distance`, `merged_id_distance`, `neighbour_concordance`, `concordance_flag` | Placement concordance ([#49](https://github.com/naturalis/timegraft/issues/49)); see [below](#placement-concordance), with the details in `backbone/concordance.tsv` |
| `placement_action` | Reserved for re-anchoring flagged genera ([#50](https://github.com/naturalis/timegraft/issues/50)); empty |

## `backbone/crown_ages.tsv`

`genus`, `crown_age` and `n_species` are as before. Added: `query_ids`, the species queried, and the same `timetree_*` columns as in `node_ages.tsv`. With `source: api` the table has no rows, as before.

## Backbone ages from a weighted fit

Before [#52](https://github.com/naturalis/timegraft/issues/52), a TimeTree age that was not younger than that of its nearest dated ancestor was dropped: the ancestor won, by its position in the tree, however weak its evidence. TimeTree ages that are too young usually come from queries that did not span the node, and those can then override better supported, older ages below them; the 27-genus crab clade at 54.9 Ma above Eubrachyura at 156.8 Ma is the example in the holistic analysis.

With `ages.method: fit` (the default) every TimeTree age is a soft constraint instead, and all backbone ages are fitted together (`timegraft.agefit`):

- **On a log scale,** so that a factor two weighs the same at 5 and at 500 Ma.
- **Weighted by the evidence.** The uncertainty of each age (`calibration_sigma`, the standard deviation of its log) comes from TimeTree's confidence interval where it gives one (taken as 95%), else from the spread of the study estimates (the standard error of their mean log), else `sigma_single` (0.5, about a factor 1.65) for a single estimate. It is never below `sigma_floor` (0.1), and is multiplied by `higher_taxa_factor` (2) for ages through enclosing taxa and by `taxonomy_only_factor` (1.5) for nodes no source tree of the synthesis supports. The weight is 1/sigma².
- **Under a hard constraint:** every age is younger than that of its nearest calibrated ancestor by at least `min_gap` (1%) per edge between them, so no branch has length zero. The fit is an exact weighted isotonic regression on the tree of calibrated nodes.
- **Robust.** A calibration more than `huber_k` (1.5) sigma from the fit is down-weighted (Huber weights, refitted until stable), so that one wild age does not drag its neighbours; beyond `outlier_z` (3) sigma it is reported as an outlier.

A configured root age (`root_age`, where TimeTree has none) is fixed. A TimeTree age equal to that of its nearest calibrated ancestor still means TimeTree does not resolve the node, and is not used. Nodes without a calibration are interpolated between the fitted ages, as before.

`ages.method: precedence` gives the rule before #52. Whichever method is used, `node_ages.tsv` also gives the precedence ages, so the two can be compared during the validation period.

| Column | Meaning |
| --- | --- |
| `source` | With the fit: `timetree_fit` or `timetree_higher_taxa_fit` for a calibrated node (with `_outlier` appended for outliers), `configured`, `interpolated` or `interpolated_unresolved_in_timetree`. With precedence, as before |
| `calibration_sigma`, `calibration_weight` | Uncertainty of the calibration as above, and its final weight after down-weighting |
| `fitted_age` | The fitted age (equal to `age`) |
| `calibration_residual`, `calibration_z` | log(fitted / TimeTree), and that over sigma |
| `calibration_outlier` | `yes` beyond `outlier_z` |
| `precedence_age`, `precedence_source` | Age and source of the node under the precedence rule, for comparison |

`metadata.json` gives conflict metrics over all calibrations, for the final ages and for the precedence ages: the median, 90th percentile and largest factor between age and TimeTree's, the number off by more than a factor 1.5 (`severe`), the root mean square of residuals over sigma, and the number of outliers.

## Catalogue of Life as identifier hub

Each genus is matched once against a fixed release of the Catalogue of Life on ChecklistBank ([#48](https://github.com/naturalis/timegraft/issues/48)), with the lineage of its reference records (`records_*`) and the kingdom of `backbone_build.context` as classification context, so that homonyms resolve to the right taxon. From the COL usage, the NCBI and OTT ids of the same name are found through ChecklistBank's names index, and every candidate's classification is compared with COL's: a candidate is only taken if it shares a taxon between phylum and tribe with COL, and the one that shares most wins.

COL is used for identifiers only. The topology still comes from the OTT and synthesis release files, and COL changes a join in two places only:

- **Which NCBI id goes to TimeTree** (`backbone_build.col.join: col`, the default). Where COL's NCBI id is one of the ids OTT merged into the genus but not the one OTT lists first, or OTT has none, the COL id is used, if COL's match is of `high` or `medium` confidence. A COL id that OTT does not hold is never used, only reported (`col_ncbi_vs_ott` = `different`): that is a disagreement for placement concordance ([#49](https://github.com/naturalis/timegraft/issues/49)) to weigh. Species ids are not changed. `join: ott` restores the earlier behaviour.
- **OTT homonyms.** Where OTT has more than one genus of the name within the context and the genus would be left `ambiguous`, the OTT id COL leads to is taken, if it is one of them (`ott_resolution` = `col`).

| Column | Meaning |
| --- | --- |
| `col_id` | COL id of the accepted taxon (for a synonym, of its accepted name) |
| `col_release` | ChecklistBank dataset key of the release, e.g. 316691 for COL26.10 XR |
| `col_match_type` | ChecklistBank's match type: `exact`, `variant` (same name, other or no authorship), `canonical`, `ambiguous`, `higherrank` (only a higher taxon fits the context), `none` |
| `col_name`, `col_rank`, `col_status` | The matched usage; `col_status` is `accepted` or a kind of synonym |
| `col_accepted_id` | For a synonym, the accepted taxon |
| `col_context` | Context the match was made with: `records` (record lineage and kingdom), `kingdom` (the record lineage did not fit COL, so only the kingdom was used), `none` |
| `col_source`, `col_sector` | `base` for a taxon from the expert-curated base release, `xr` for one added in the extended release; the COL sector it came from |
| `col_family`, `col_order`, `col_class`, `col_phylum` | COL's classification |
| `col_alternatives` | For an `ambiguous` match, the COL ids it could not choose between. The genus is left unresolved |
| `col_confidence` | `high`: matched with the record lineage, accepted, from the base release, one NCBI candidate. `medium`: matched only with the kingdom, a synonym, or the NCBI id needed the classification to choose between homonyms. `low`: added in XR only. Empty: unresolved. Only `high` and `medium` change a join; later scores can weigh the rest |
| `col_ncbi_id`, `col_ott_id` | NCBI and OTT ids derived from the COL taxon. The OTT ids come from OTT 3.7.2 in ChecklistBank and are compared with, not used instead of, the release in use |
| `col_ncbi_candidates`, `col_ott_candidates` | All ids of the name in the names index |
| `col_ncbi_path`, `col_ott_path` | `names_index` (one candidate, classification agrees), `names_index+classification` (several, one agrees best), `disagrees` (none shares a taxon with COL), `ambiguous` (several agree equally; left unresolved), `none` |
| `col_ncbi_vs_ott` | How `col_ncbi_id` relates to OTT's ids: `same` as the one OTT uses, `merged` (another id OTT holds), `different`, or empty |
| `col_ott_vs_ott` | Whether `col_ott_id` is the OTT id of the match: `same`, `different`, or empty |
| `join_ncbi_id`, `join_ncbi_source` | The NCBI id the genus was sent to TimeTree with, and whether it is OTT's (`ott`) or COL's (`col`) |

With `backbone_build.col.species: true`, every species name of the records is resolved too, into `backbone/col_species.tsv` (`genus`, `species`, the COL match and confidence, the NCBI id derived from it, and the NCBI id OTT has for the name, `ott_ncbi_id`). This is one query per name and off by default; it is meant for crown ages from the sampled species ([#54](https://github.com/naturalis/timegraft/issues/54)). Otherwise the file has only its header.

**Reproducibility.** `col.release` is a ChecklistBank dataset key, by default 316691 (COL26.10 XR, 16 October 2026). An alias such as `3LXR` also works but follows the latest release; the key it resolved to is in `metadata.json` and the log warns. Every ChecklistBank answer is cached in `backbone_build.cache` under keys starting with `col:` and the release key, so a run with the same release and cache gives the same result without network access.

## Anchoring genera that are not a node of the synthesis

A genus that is not a node of the synthetic tree (in Iberia, 149 of 1,218) is anchored on its species there. The MRCA of all of them, the rule before [#51](https://github.com/naturalis/timegraft/issues/51), lets a single stray species (a misidentification, a species that belongs elsewhere, a regional subset of a global genus) pull the anchor, and the stem age with it, deep into the tree: *Carassius* hung at the base of a clade of 277 genera because two of its five species sit elsewhere. With `anchoring.rule: majority` (the default) the anchor is instead the smallest clade holding more than `anchoring.share` (default half) of the species; clades holding more than half of the same species are nested, so that clade is unique. `share` must therefore be at least 0.5, and below 1, since no clade holds more than all of them; other values are refused.

**Which species.** Those sampled in the run (the species names of the genus's records found in OTT), where at least `min_sampled` (default 2) of them are in the synthesis: the genus subtree is built from those units, so that is the clade it should hang from. Otherwise all species of the genus in OTT.

**Fallbacks.** With fewer than `min_species` (default 3) species in the synthesis there is no majority to speak of, and the MRCA is used (`too_few_species`). Where no clade smaller than the MRCA holds a majority, because the species are scattered, the MRCA is used too (`no_smaller_majority`). With a single species, the genus hangs from its parent, as before. Genera that are a node of the synthesis, or that are only attached at an enclosing taxon, are not affected. `anchoring.rule: mrca` with `prefer_sampled: false` gives the rule before #51.

| Column | Meaning |
| --- | --- |
| `anchor_rule` | `majority` or `mrca`: which rule set the anchor. Empty for genera that are a node of the synthesis or attached at an enclosing taxon |
| `anchor_fallback` | Why the MRCA was used under the majority rule: `too_few_species`, `no_smaller_majority`; empty otherwise |
| `anchor_basis` | `sampled` (the run's species) or `all` (all species in OTT) |
| `anchor_species`, `anchor_species_inside` | Species used, and those in the clade anchored at |
| `anchor_outside_share`, `anchor_outside_species` | The share of them left outside that clade, and their names (at most ten) |
| `anchor_mrca_node`, `anchor_mrca_clade_genera` | The MRCA of all of them and the genera in its clade: where the earlier rule would have put the genus, against `anchor_node` and `anchor_clade_genera` for where it is now |

The status of a genus anchored on a majority clade is `attached_at_majority`. For a re-anchored genus (#50) these columns describe its original anchoring.

## Placement concordance

Is each genus where OToL puts it? Three checks, each independent of OToL's placement and of the others, test that before any age is derived from the placement ([#49](https://github.com/naturalis/timegraft/issues/49)). The checks themselves change nothing; what is done with a flagged genus is described [below](#re-anchoring-flagged-genera) ([#50](https://github.com/naturalis/timegraft/issues/50)).

| Check | Question | Distance |
| --- | --- | --- |
| Reference lineage | Where in OTT is the family that a reference classification gives the genus? The reference is COL's classification where COL resolved the genus, otherwise the most common family among its records (`reference_source`); without a family in OTT, the order, then the class | Rank level at which the genus's OTT lineage and that taxon's part |
| Merged ids | Where OTT merged more than one NCBI id into the genus (`ott_ncbi_ids`), how far apart are they in NCBI Taxonomy (read from ChecklistBank)? | Rank level at which the widest pair parts |
| Barcode neighbours | Do the genera the genus's units come closest to outside the genus (`best_other_genus` in `units/screened.tsv`) share a phylum with its sisters in the backbone? Phyla from COL, else from the records | Share of those hits in a phylum of the sisters, with at least `min_hits` hits |

A rank level is 0 where the lineages agree (the genus lies within the reference family, or the ids are of one taxon), then 1 family, 2 superfamily, 3 order, 4 class, 5 phylum and 6 kingdom: the highest rank among the taxa below the last shared one, so two phyla that meet at an unranked node such as Bilateria still part at phylum level. Where a reference name is held by more than one OTT taxon, the closest is used, so a homonym never makes a genus look misplaced.

A lineage or merged-id check fails when the level is at least `concordance.flag_level` (default `class`); the neighbour check fails when the share is at most `max_share` (default 0, no hit at all in a sister phylum). Each check alone fails for reasons other than misplacement (classification differences, slow mitochondrial evolution in sponges, contaminated units), and those reasons do not coincide; a genus is therefore flagged when `min_flags` (default 2) checks fail.

`backbone/concordance.tsv` has one row per genus, sorted for review (most failed checks first, then the widest split):

| Column | Meaning |
| --- | --- |
| `rank` | Position in the review order |
| `genus`, `status`, `ott_id`, `ott_family`, `ott_class`, `ott_phylum` | The genus and its OTT placement |
| `reference_source`, `reference_taxon` | `col` or `records`, and the reference taxon used, as `rank:name` |
| `lineage_mrca`, `lineage_split`, `lineage_distance`, `lineage_flag` | Last OTT taxon the placement and the reference taxon share, the rank they part at, its level, and whether the check failed (`yes`) |
| `merged_ncbi_ids`, `merged_mrca`, `merged_split`, `merged_id_distance`, `merged_flag` | The same for the NCBI ids OTT merged into the genus; empty where it has one |
| `neighbour_hits`, `neighbour_phyla`, `sister_phyla`, `neighbour_concordance`, `neighbour_flag` | Units with a nearest unit in another genus, the phyla of those genera with counts, the phyla of the backbone sisters, the share in a sister phylum, and whether the check failed |
| `n_evaluable`, `n_flags`, `severity` | Checks that could be made and that failed; the widest split level found |
| `concordance_flag` | `flag` (at least `min_flags` checks failed), `review` (one failed), empty (none failed) or `not_evaluable` (no check could be made, e.g. the genus is not in OTT) |

The four columns in `genera.tsv` are `lineage_distance`, `merged_id_distance`, `neighbour_concordance` and `concordance_flag`. How these settings behave on the Iberian run, with every genus each check fails on, is in `data/iberia/consistency/README.md`.

## Re-anchoring flagged genera

A genus that placement concordance flags sits in the wrong clade of the synthetic tree, and its stem age, the ages above it, its own outgroup and those of its neighbours all follow from that placement. It is therefore moved before any of these are derived ([#50](https://github.com/naturalis/timegraft/issues/50)); the TimeTree queries, the interpolation and the outgroups are then built from the corrected backbone as for any other genus, and the enclosing taxa used for TimeTree queries become those of the new place.

**When.** With `reanchor.mode: apply` (the default), a genus is moved when its `concordance_flag` is `flag` and its `severity` is at least `reanchor.min_severity` (default `class`). `mode: report` makes the same decisions and records them (`would_reanchor`) without changing anything; `mode: off` skips them (written unquoted, YAML reads `off` as false, which is taken to mean the same). Re-anchoring needs `concordance.enabled`.

**Where to.** To the most specific taxon of the genus's reference lineage that OTT holds (the reference family, else order, else class), taking only an OTT taxon of that name whose own lineage contains a higher name of the reference lineage. A name that OTT holds only elsewhere, or more than once without exactly one copy agreeing, is passed over for the next rank, so a homonym is never chosen; if no rank qualifies, the genus stays where it is and is recorded as `flagged_no_target`. A target given as `ott:<id>` in the curated table must exist in the loaded taxonomy. The genus is attached in the synthetic tree at that taxon's node; where the taxon is not a node of the synthesis, at the MRCA of its members there; failing that, at its nearest ancestor that is. The genus's own synthesis node and species are left out of that search, since they are what is misplaced.

**Curated table.** `reanchor.overrides` names a TSV with columns `genus`, `action`, `target`, `source` and `note`. `action` is `keep` (never move this genus, whatever the checks say) or `reanchor` (move it to `target`, given as `rank:name` such as `family:Onchidiidae` or as `ott:<id>`, flagged or not). `source` and `note` record who decided and why. `config/placement_overrides.tsv` is an empty table to start from.

Columns in `genera.tsv`:

| Column | Meaning |
| --- | --- |
| `placement_action` | `reanchored`; `would_reanchor` (mode `report`); `kept_override` (flagged, but `keep` in the table); `flagged_no_target` (flagged, but no taxon of the reference lineage could be found or placed); empty for every other genus |
| `original_status`, `original_anchor` | The `status` and synthesis node of the genus before the move. After a move, `status` is `reanchored:<node>` and `anchor_node` is the new node |
| `reanchor_target` | The taxon moved to, as `rank:name (ott<id>)` |
| `reanchor_how` | How it was chosen (`reference_lineage` or `override`) and how it was placed (`taxon_in_synthesis`, `mrca_of_members`, `ancestor_in_synthesis`), separated by `;`. For `flagged_no_target`: `no_agreeing_target` (OTT holds the reference names, but none uniquely with a lineage that agrees; `reanchor_target` is the last one tried), `no_target` (OTT holds none of them) or `override_not_found` |
| `reanchor_trigger` | The concordance result that triggered it: checks failed and made, severity and the three metrics |
| `reanchor_source`, `reanchor_note` | `auto`, or `override:<source>` with the note from the table |

`ott_family` to `ott_phylum` keep OTT's own lineage, which is what was wrong. In `concordance.tsv`, `placement_action` repeats the action, and `neighbour_concordance_after` gives the barcode-neighbour check again on the corrected topology, for every genus, so that a correction can be seen to make the neighbours agree (for the genus and for its former sisters).

## `backbone/metadata.json`

New. What the backbone was built from: the format version of each output and of the cache, the timegraft version, the source, OTT and synthesis releases, the matching context, whether synthesis support was read, the age method with the conflict metrics of the final and the precedence ages, the anchoring, concordance and re-anchoring settings with the number of genera per status and the number of genera per `concordance_flag` and per `placement_action`, the COL release (key, alias, version) with counts of genera by confidence and crosswalk path and of genera joined with COL's NCBI id, and a summary of the TimeTree evidence (queries per evidence level, TimeTree database versions, and how many cache entries are at version 1, were upgraded, or now get a different answer from TimeTree).

## Migration note for downstream consumers

- **Read columns by name.** Every table keeps its earlier columns in their earlier order, and new columns follow them. Code that reads by name, as everything in timegraft does, needs no change. Code that reads by position keeps working for the old columns; code that checks the number of columns does not.
- **Reserved columns are empty** until the issue that fills them lands on the `consistency-layer` branch. Do not treat an empty value as zero.
- **The cache can be shared** between old and new versions. The first run of a new version on an old cache asks TimeTree again for every query it uses (for the Iberian run, a few minutes at four requests at a time) and writes version 2 entries; set `upgrade_cache: false` to avoid this, for example offline.
- **Synthesis support needs one more download**, `annotations.json` of the synthesis release (61 MB for opentree16.1), placed next to `labelled_supertree.tre`. Set `synthesis_support: false` to do without; the `synth_*` support columns are then empty.
- **COL needs ChecklistBank** on the first run with a release (for 1,200 genera, about six requests each), and changes which NCBI id some genera are sent to TimeTree with, so some backbone ages change ([#48](https://github.com/naturalis/timegraft/issues/48)). `col.join: ott` keeps the identifiers as before while still recording COL's; `col.enabled: false` skips COL altogether.
- **`genera.tsv` format 3** (#48) replaces the three reserved `col_*` columns of format 2 with the filled ones above, adds `ott_resolution` and the `join_*` columns, and keeps the earlier columns and their order.
- **`genera.tsv` format 4** (#49) fills `lineage_distance`, `merged_id_distance`, `neighbour_concordance` and `concordance_flag`; `backbone/concordance.tsv` is new. The build now also reads `units/screened.tsv`, so it runs after screening.
- **`genera.tsv` format 5 and `concordance.tsv` format 2** (#50) add the re-anchoring columns above. With the default `reanchor.mode: apply` the backbone changes wherever a genus is flagged: the genus has a new place and `status` (`reanchored:<node>`), and the ages around its old and new place change. `mode: report` keeps the earlier backbone.
- **`genera.tsv` format 6** (#51) adds the `anchor_*` columns above, and a new status, `attached_at_majority`. With the default `anchoring.rule: majority` and `prefer_sampled: true`, genera that are not a node of the synthesis can have a different anchor, and the backbone and its ages change around them.
- **`node_ages.tsv` format 3** (#52) fills the calibration columns and adds `calibration_sigma`, `calibration_z`, `calibration_outlier`, `precedence_age` and `precedence_source`. With the default `ages.method: fit`, `source` takes new values (`timetree_fit`, `timetree_higher_taxa_fit`, with `_outlier`), `dropped_inconsistent` and `interpolated_after_drop` no longer occur, and backbone ages change wherever TimeTree ages conflicted. Code that selected TimeTree-dated nodes by `source in ("timetree", "timetree_higher_taxa")` should accept the `_fit` forms too.
- **The build reads the filtered records** (`units/filtered.tsv`) for the reference lineage, besides the representatives, so it runs again when the records change.
