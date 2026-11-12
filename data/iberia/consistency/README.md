# Consistency layer: evidence on the Iberian run

Material for the consistency layer ([#46](https://github.com/naturalis/timegraft/issues/46)) on the Iberian run. Each issue of the epic adds what it needs here, so that the benchmark of [#56](https://github.com/naturalis/timegraft/issues/56) can compare every step with the same starting point.

## `backbone/` (#47)

The backbone of `iberia-run-tolerant.tar.gz` rebuilt with the evidence columns of [#47](https://github.com/naturalis/timegraft/issues/47) (described in `docs/technical/provenance.md`): `genera.tsv`, `node_ages.tsv`, `crown_ages.tsv`, `metadata.json` and `backbone.tre`.

How it was made, on 12 November 2026: `workflow/scripts/build_backbone.py` was run on `units/representatives.tsv` of the archive, the records of `records.tsv.gz` without those in its `units/rejected.tsv` (the workflow's `units/filtered.tsv`), OTT 3.7.3, synthesis 16.1 with its `annotations.json`, and `timetree_cache.json` as it was, with `upgrade_cache: true`. The run upgraded all 2,507 cache entries to version 2; TimeTree (database version 20260820) gave the same age and found the same taxa for every one of them. The upgraded cache is now `../timetree_cache.json`.

Checks that nothing changed:

- The same script at `main` (commit `f4f3d7f`), on the same inputs, reproduces the backbone files of the archive byte for byte.
- Against that, `backbone.tre` is identical, and the earlier columns of `node_ages.tsv` (779 nodes), `genera.tsv` (1,218 genera) and `crown_ages.tsv` (312 genera) are identical in content and order; only columns were added.

The new columns reproduce the counts of the holistic analysis: of the 314 TimeTree node ages used, 73 rest on one estimate, 89 on two or three and 152 on four or more; 218 of 779 backbone nodes have no source tree in the synthesis; 38 genera carry more than one NCBI id in OTT; 20 genera attached at the MRCA of their species sit on a clade of more than 50 genera. For *Onchidella*, `genera.tsv` shows OTT's lineage (Laganidae, Echinoidea) against that of all 82 of its records (Onchidiidae, Gastropoda, Mollusca), and `node_ages.tsv` shows that the two dropped ages above it were TimeTree's age of NCBI 33213, Bilateria.
