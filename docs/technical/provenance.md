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
| `conflict_with` | For a node in conflict (`conflict_ratio`), the node it conflicted with: the dated ancestor, or for a dropped `timetree_higher_taxa` age, the descendant that outranked it |
| `timetree_evidence` | `full` (version 2 cache entry), `legacy_cache` (version 1, only `timetree_found_ids` known), `none`; empty if TimeTree gave no age |
| `timetree_found_ids` | The queried NCBI ids TimeTree found |
| `timetree_mrca_ncbi` | NCBI id of the node TimeTree dated |
| `timetree_n_estimates`, `timetree_n_studies` | Number of time estimates and of studies behind the age |
| `timetree_estimate_min`, `timetree_estimate_max` | Range of the estimates; their ratio is the spread used in the [holistic analysis](https://github.com/naturalis/timegraft/issues/46) |
| `timetree_ci_low`, `timetree_ci_high` | TimeTree's confidence interval |
| `timetree_api_version` | Version of the TimeTree database that answered |
| `synth_node` | Label of the node of the synthetic tree behind the backbone node (`ott…` or `mrcaott…ott…`); static source only |
| `synth_supported_by`, `synth_conflicts_with`, `synth_resolves`, `synth_partial_path_of`, `synth_terminal` | Number of source trees of the synthesis per kind of annotation, from its `annotations.json`. `synth_supported_by` of 0 means the node rests on the taxonomy alone. Empty without `backbone_build.synthesis_support` |
| `calibration_weight`, `fitted_age`, `calibration_residual` | Reserved for the weighted fit of backbone ages ([#52](https://github.com/naturalis/timegraft/issues/52)); empty |
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
| `records_n` | Reference records of the genus |
| `records_family`, `records_order`, `records_class`, `records_phylum` and their `_share` | The most common name per rank among those records, and the share of records naming that rank that agree. This is the reference lineage that placement concordance compares with OTT ([#49](https://github.com/naturalis/timegraft/issues/49)) |
| `col_*`, `join_ncbi_id`, `join_ncbi_source` | Catalogue of Life as identifier hub; see [below](#catalogue-of-life-as-identifier-hub) |
| `lineage_distance`, `merged_id_distance`, `neighbour_concordance`, `concordance_flag` | Reserved for placement concordance ([#49](https://github.com/naturalis/timegraft/issues/49)); empty |
| `placement_action` | Reserved for re-anchoring flagged genera ([#50](https://github.com/naturalis/timegraft/issues/50)); empty |

## `backbone/crown_ages.tsv`

`genus`, `crown_age` and `n_species` are as before. Added: `query_ids`, the species queried, and the same `timetree_*` columns as in `node_ages.tsv`. With `source: api` the table has no rows, as before.

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

## `backbone/metadata.json`

New. What the backbone was built from: the format version of each output and of the cache, the timegraft version, the source, OTT and synthesis releases, the matching context, whether synthesis support was read, the COL release (key, alias, version) with counts of genera by confidence and crosswalk path and of genera joined with COL's NCBI id, and a summary of the TimeTree evidence (queries per evidence level, TimeTree database versions, and how many cache entries are at version 1, were upgraded, or now get a different answer from TimeTree).

## Migration note for downstream consumers

- **Read columns by name.** Every table keeps its earlier columns in their earlier order, and new columns follow them. Code that reads by name, as everything in timegraft does, needs no change. Code that reads by position keeps working for the old columns; code that checks the number of columns does not.
- **Reserved columns are empty** until the issue that fills them lands on the `consistency-layer` branch. Do not treat an empty value as zero.
- **The cache can be shared** between old and new versions. The first run of a new version on an old cache asks TimeTree again for every query it uses (for the Iberian run, a few minutes at four requests at a time) and writes version 2 entries; set `upgrade_cache: false` to avoid this, for example offline.
- **Synthesis support needs one more download**, `annotations.json` of the synthesis release (61 MB for opentree16.1), placed next to `labelled_supertree.tre`. Set `synthesis_support: false` to do without; the `synth_*` support columns are then empty.
- **COL needs ChecklistBank** on the first run with a release (for 1,200 genera, about six requests each), and changes which NCBI id some genera are sent to TimeTree with, so some backbone ages change ([#48](https://github.com/naturalis/timegraft/issues/48)). `col.join: ott` keeps the identifiers as before while still recording COL's; `col.enabled: false` skips COL altogether.
- **`genera.tsv` format 3** (#48) replaces the three reserved `col_*` columns of format 2 with the filled ones above, adds `ott_resolution` and the `join_*` columns, and keeps the earlier columns and their order.
- **The build reads the filtered records** (`units/filtered.tsv`) for the reference lineage, besides the representatives, so it runs again when the records change.
