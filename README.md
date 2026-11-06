# timegraft

Dated barcode trees for phylogenetic diversity: per-genus subtrees from BOLD or UNITE, grafted onto a dated OToL/TimeTree backbone.

timegraft builds a time-calibrated tree with barcode reference units at the tips: BINs for COI, UNITE species hypotheses for ITS. Branch lengths are in millions of years, so phylogenetic diversity (Faith's PD, MPD, UniFrac, PhyloSor) is comparable across sites, studies and markers.

## The idea

Short barcodes cannot resolve deep relationships: COI saturates and ITS cannot be aligned beyond genus level. timegraft therefore never estimates deep structure from barcodes.

- **The backbone is given.** Topology and ages down to genus level come from external sources (Open Tree of Life, TimeTree, curated synthesis trees). Each genus tip carries its stem age S.
- **Barcodes only fill in each genus.** Within a genus, reference units are joined with BioNJ, rooted on units from the sister genus, and dated so that the genus crown age C lies below S.
- **Grafting is a fixed attachment.** Each genus tip is replaced by a stem of length S minus C followed by its dated subtree. There is no global tree search, and every genus is processed independently.

Everything that differs between markers sits behind four plug-ins, configured per marker in `config/config.yaml`:

| Plug-in | COI | ITS |
| --- | --- | --- |
| Reference data and tip unit | BOLD records; BINs | UNITE records; species hypotheses |
| Filter | Length, ambiguity, reading frame under the group's genetic code | Length, ambiguity (ITSx extraction upstream) |
| Distances within a genus | K2P on aligned sequences | k-mer (Mash) distance, no alignment |
| Crown age when no known age is given | Substitution-rate prior, capped below S | Stated fraction of S |

Within each genus, the rooted subtree is made ultrametric with its root at the crown age, by mean path length (built in, the default) or by LSD2 (`dating.method: lsd2`; the `build_lsd2` rule compiles a pinned version of [LSD2](https://github.com/tothuhien/lsd2) once, in a conda environment with a C++ compiler).

Known crown ages (for instance from TimeTree) always take precedence. The provenance table records where every crown age came from.

## Status

Milestone 1, working on synthetic data:

- [x] Record filtering, one medoid representative per unit, genus assignment with conflict detection
- [x] Sister-genus outgroups from the backbone, BioNJ, rooting, mean-path-length dating, grafting
- [x] Sharded execution so genera are processed in parallel
- [x] Toy data set, unit tests and an end-to-end workflow test
- [x] Backbone assembly from static OToL release files (or the OToL API), dated with TimeTree node ages, the rest interpolated
- [x] Genus crown ages from TimeTree, via species within each genus
- [x] LSD2 (least-squares dating) as an alternative to mean path length; treePL is not packaged for conda, LSD2 is built from source by the workflow
- [ ] Embedding distances as a third distance method
- [x] Matching ASVs to units (VSEARCH) and placing unmatched ASVs on the tree near their nearest unit
- [x] Placement of reference units without a genus, next to their nearest unit
- [x] Per-sample PD metrics: Faith's PD, MPD (presence and abundance-weighted), unweighted and weighted UniFrac, PhyloSor
- [x] Validation of the informative window: barcode distances against a reference tree, per taxonomic stratum, with a saturation plot
- [ ] Using the validation to set the informative window per clade automatically

## Usage

```bash
pip install -e ".[dev]"
snakemake --cores 4                 # runs the toy example in config/config.yaml
snakemake --cores 4 --use-conda     # same, with the conda environment per rule
pytest
```

### Inputs

- `records`: tab-separated, with columns `record_id`, `unit_id`, `genus`, `sequence` (optionally `family`, `order` and other ranks, used by the validation).
- Or `bcdm.path`: a BOLD data package in BCDM format (TSV or Parquet; Parquet is read with DuckDB and is much faster on a full dump), filtered by marker, taxa and country into a record table. The same filter runs on its own, which is handy for cutting a regional subset out of a full dump: `python -m timegraft.bcdm BOLD_Public.parquet records.tsv --taxon phylum=Annelida,Mollusca --country Portugal,Spain`. Use `--exclude class=Insecta` to drop taxa, `rank=@file` to read names from a file (one per line), and `--max-per-unit 10` to keep at most ten records per BIN, which is plenty for choosing a representative and keeps the table small. For `align: none` the sequences must already be aligned; set `align: mafft` otherwise.
- `backbone`: Newick, ultrametric, ages in Ma, one tip per genus. Alternatively set `backbone_build.enabled: true` to build it for all genera in the records (see below).

### Building the backbone

With `backbone_build.enabled: true`:

1. **Topology.** By default (`source: static`) from Open Tree of Life release files: the OTT taxonomy and the labelled synthetic tree, downloaded once into `resources/otol` by the `download_otol` rule (about 140 MB). Genera are matched by name within `context` (e.g. `Metazoa`, `Fungi`) to avoid homonyms. A genus that is not a node in the synthesis is attached at the MRCA of its species, or else at its nearest enclosing taxon that is. `source: api` uses the OToL web API instead.
2. **Node ages.** Each internal node is dated with a TimeTree MRCA age, from a few genera on each side of it, sent with some of their species, since TimeTree knows many species but not every genus. Other genera and species are tried when TimeTree knows none of the first picks. As a last resort, enclosing taxa (family up to phylum) stand in for a side, but only taxa that contain no genus from the other side, so the MRCA is still the same node. When a node's TimeTree age is not younger than that of its nearest dated ancestor, equal ages are read as TimeTree not resolving the node; otherwise a direct TimeTree age outranks one obtained via enclosing taxa, and between equals the ancestor is kept. `node_ages.tsv` gives the ratio for every conflict: large ratios often point to a genus misplaced in the synthetic tree (in the first real run, a gastropod genus that OTT had merged with an echinoderm homonym). Remaining nodes are interpolated.
3. **Genus crown ages.** The TimeTree MRCA of up to twelve species per genus. Such an estimate is far too old when the genus is not monophyletic or holds misplaced species, so it is rejected when it is not younger than the stem age or, for COI, more than `timetree_max_ratio` (default 3) times the age implied by the barcodes and the substitution rate. `provenance.tsv` shows the TimeTree value and whether it was used. A `crown_ages` table of your own always takes precedence.

Needs network access to `files.opentreeoflife.org` (or `api.opentreeoflife.org`) and `api.timetree.org`. All responses are cached in `backbone_build.cache`. The build writes `backbone/genera.tsv` (how each genus was matched and placed), `backbone/node_ages.tsv` (each node's age and its source) and `backbone/crown_ages.tsv`.

Known limitation: TimeTree coverage of marine invertebrates is patchy (e.g. many polychaetes and amphipods), so some nodes are interpolated, and interpolation spaces nodes evenly, which tends to make shallow undated nodes too old.
- `crown_ages` (optional): tab-separated, `genus`, `crown_age`.
- `samples` (optional): tab-separated, `sample`, `unit_id`, `abundance`. When given, the workflow also computes diversity metrics.
- `asvs` and `asv_table` (optional, instead of `samples`): ASV sequences (FASTA) and a table `sample`, `asv_id`, `reads`. Each ASV is matched to a unit at `match_identity` or above; below that, down to `place_identity`, it is attached on the path from its nearest unit to the root, at the age implied by its distance (COI: identity and the substitution rate; ITS: halfway along the nearest unit's terminal edge). Needs VSEARCH.

### Outputs (under `outdir`)

- `timegraft.tre`: the dated tree, tips are unit identifiers. Units without a usable genus are placed next to their nearest unit (`units_placed.tsv`); `grafted.tre` is the tree before that step.
- `provenance.tsv`: per backbone genus, whether it received a subtree, a single unit or was pruned, with stem age, crown age and its source, the rooting method and the outgroup used.
- `placement/placements.tsv`: per ASV, matched, placed or unplaced, with the unit, identity and attachment age; `placement/timegraft_placed.tre`: the tree with placed ASVs added.
- `diversity/alpha.tsv`: per sample, richness, Faith's PD (including the path to the root), MPD and abundance-weighted MPD, in Ma.
- `diversity/beta.tsv`: per pair of samples, unweighted and weighted (normalised) UniFrac and PhyloSor.
- `diversity/units_not_in_tree.tsv`: sample units that are not in the tree, so they are visible rather than silently dropped.
- `validation/` (when `validation.reference_tree` is set): sampled pairs, rank correlations per stratum, binned medians, the reference distance at which each barcode distance plateaus, and `saturation.png`.
- `units/rejected.tsv`, `units/unassigned.tsv`, `genera/missing_from_backbone.tsv`: records and units that did not make it into the tree, and why.

## Layout

```
config/config.yaml        settings and marker plug-ins
workflow/Snakefile        entry point
workflow/rules/           units, per-genus subtrees, grafting
workflow/scripts/         thin wrappers around the package
src/timegraft/            the package: filters, distances, BioNJ, rooting, dating, grafting
tests/                    unit tests, workflow test, toy data
```

## Licence

Code: Apache-2.0. The trees and tables timegraft produces are released into the public domain (CC0 1.0). The design builds on lessons from [Bactria](https://github.com/naturalis/barcode-constrained-phylogeny) and [evospaice](https://github.com/naturalis/evospaice).
