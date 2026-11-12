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
| `col_id`, `col_release`, `col_match_type` | Reserved for Catalogue of Life XR as identifier hub ([#48](https://github.com/naturalis/timegraft/issues/48)); empty |
| `lineage_distance`, `merged_id_distance`, `neighbour_concordance`, `concordance_flag` | Reserved for placement concordance ([#49](https://github.com/naturalis/timegraft/issues/49)); empty |
| `placement_action` | Reserved for re-anchoring flagged genera ([#50](https://github.com/naturalis/timegraft/issues/50)); empty |

## `backbone/crown_ages.tsv`

`genus`, `crown_age` and `n_species` are as before. Added: `query_ids`, the species queried, and the same `timetree_*` columns as in `node_ages.tsv`. With `source: api` the table has no rows, as before.

## `backbone/metadata.json`

New. What the backbone was built from: the format version of each output and of the cache, the timegraft version, the source, OTT and synthesis releases, the matching context, whether synthesis support was read, and a summary of the TimeTree evidence (queries per evidence level, TimeTree database versions, and how many cache entries are at version 1, were upgraded, or now get a different answer from TimeTree).

## Migration note for downstream consumers

- **Read columns by name.** Every table keeps its earlier columns in their earlier order, and new columns follow them. Code that reads by name, as everything in timegraft does, needs no change. Code that reads by position keeps working for the old columns; code that checks the number of columns does not.
- **Reserved columns are empty** until the issue that fills them lands on the `consistency-layer` branch. Do not treat an empty value as zero.
- **The cache can be shared** between old and new versions. The first run of a new version on an old cache asks TimeTree again for every query it uses (for the Iberian run, a few minutes at four requests at a time) and writes version 2 entries; set `upgrade_cache: false` to avoid this, for example offline.
- **Synthesis support needs one more download**, `annotations.json` of the synthesis release (61 MB for opentree16.1), placed next to `labelled_supertree.tre`. Set `synthesis_support: false` to do without; the `synth_*` support columns are then empty.
- **The build reads the filtered records** (`units/filtered.tsv`) for the reference lineage, besides the representatives, so it runs again when the records change.
