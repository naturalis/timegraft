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
