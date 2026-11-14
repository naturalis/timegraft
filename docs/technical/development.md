# Developing timegraft

## Layout

```
config/config.yaml            settings, including one block per marker
workflow/Snakefile            entry point; decides which optional rules run
workflow/rules/*.smk          tools, units, backbone, genus, graft, placement, diversity, validation
workflow/scripts/*.py         thin wrappers: read inputs, call the package, write outputs
workflow/envs/*.yaml          conda environments
src/timegraft/                the package (all logic lives here, so it can be tested without Snakemake)
tests/                        unit tests, a toy workflow test, small fixtures
docs/technical/                using and developing timegraft
docs/science/                  background, method, limitations, results
docs/manuscript/               draft preprint (Markdown + BibTeX)
```

Rule scripts import the package from `src/` through the `src` parameter that every rule passes (`params.src`), so the workflow runs without installing the package into each conda environment.

## Modules

| Module | Responsibility |
| --- | --- |
| `seqio` | TSV and FASTA reading and writing |
| `bcdm` | BOLD data packages (TSV, Parquet via DuckDB) to record tables |
| `filters` | marker filters (COI reading frame; ITS length and ambiguity) |
| `geneticcode` | NCBI translation tables, translation, the genetic code per record, and building the taxon-to-code mapping from the NCBI taxdump |
| `units` | medoid representative and genus per unit |
| `distances` | distance methods (`k2p`, vectorised; `kmer`); the registry `METHODS` |
| `bionj` | BioNJ on a distance matrix |
| `rooting` | rooting genus trees: tolerant (minimum violation) or strict outgroup rule, MAD and midpoint; other plausible roots ([#34](https://github.com/naturalis/timegraft/issues/34)) |
| `rooting_diagnostics` | after a run: why outgroup rooting failed, midpoint and MAD against the outgroup root, outgroup and distance variants ([#19](https://github.com/naturalis/timegraft/issues/19)) |
| `dating` | crown age within bounds, mean path length, LSD2, chronos (optional) |
| `rates` | substitution rate per genus from its taxonomy |
| `screening` | units that do not fit their genus (`screen`), from an all-against-all search |
| `calibrate` | rates per class from genera of known crown age (`python -m timegraft.calibrate`) |
| `build` | one dated subtree per genus |
| `backbone` | reading a dated backbone, stem ages, sister genera and outgroup units, sharding |
| `otol_static` | OTT taxonomy and synthetic tree from release files; genus topology; species and lineages; match, anchor and synthesis-support evidence |
| `opentree` | the same via the OToL web API |
| `col` | Catalogue of Life on ChecklistBank: matching with classification context, crosswalk to NCBI and OTT ids through the names index |
| `concordance` | placement concordance: rank-level distances between lineages, barcode neighbours against backbone sisters, flags |
| `reanchor` | moving flagged genera before dating: targets from the reference lineage, synthesis nodes, curated overrides |
| `agefit` | weighted, robust fit of backbone ages under monotonicity: calibration uncertainty, tree isotonic regression, conflict metrics |
| `requery` | re-querying conflicting nodes with disjoint picks: triggers, picks, classification |
| `timetree` | TimeTree MRCA queries with caching and concurrency; cache format version 2 with the evidence behind each age ([provenance](provenance.md)) |
| `backbone_build` | dating a genus topology with TimeTree; conflicts; interpolation; the columns of `node_ages.tsv` and `genera.tsv` |
| `graft` | grafting subtrees onto the backbone |
| `placement` | VSEARCH search, classification, attaching sequences on the tree |
| `diversity` | Faith's PD, MPD, UniFrac, PhyloSor |
| `validation` | pairs per stratum, correlations, binned medians, plateau onset, plot |
| `webapi` | HTTP with retries and a JSON cache |

## Plug-in points

Marker-specific behaviour is selected by name from the marker block in the configuration:

- **Filter**: `filters.CHECKS` maps a marker name to a check function returning `None` or a reason.
- **Distance**: `distances.METHODS` maps a name to a pairwise function; `NEEDS_ALIGNMENT` says whether it needs aligned input; `distances.matrix` can take a vectorised fast path (as for `k2p`).
- **Crown age**: `dating.crown_age` with `method` `rate` or `fraction`; bounds from the stem age, TimeTree and the crown age table; the plausibility test of TimeTree crowns in `dating.timetree_plausibility`.
- **Dating**: `mpl` or `lsd2`, chosen in `build.build_genus`.
- **Placement height**: `rate` or `terminal`, in `placement.place_all`.

## Adding a distance method

1. Write a function `name(a: str, b: str, **options) -> float`, or better a matrix function for speed, and register it in `distances.METHODS` and `NEEDS_ALIGNMENT`.
2. A method that needs precomputed data (for example vectors per unit) should read them from a file passed as an option, made in a separate rule, rather than compute them on the fly.
3. Add it to `validation.methods` first and compare it with `k2p` and `kmer` before using it for trees.

Embeddings from DNA language models were considered and decided against ([#25](https://github.com/naturalis/timegraft/issues/25)): they separate species, but beyond short distances they are no better than alignment-based or k-mer distances, and they need a GPU.

## Adding a marker (e.g. ITS)

A marker is a block under `markers` plus, where needed, new plug-in functions. For ITS, see [issue #15](https://github.com/naturalis/timegraft/issues/15) and the [roadmap](roadmap.md); the main pieces are a UNITE reader (species hypotheses as units), ITSx extraction, a check of k-mer distances within genera, crown ages without a rate, `context: Fungi` for the backbone, and a fungal test set.

## Testing

```bash
pytest                      # all tests, including the toy workflow (skipped if snakemake or vsearch is missing)
snakemake --lint            # workflow style checks
```

- `tests/data/toy/`: a synthetic data set (made by `tests/data/make_toy.py`): five genera, one without units, a unit without a genus, records that must be rejected, ASVs covering every placement outcome, and three samples.
- `tests/data/otol/`: a tiny OTT taxonomy and synthetic tree, with a homonym, a genus only present through its species, and one missing altogether.
- Network calls are tested with fake fetch functions; nothing in the tests needs network access.

CI (GitHub Actions) installs VSEARCH, builds LSD2, runs the lint, a dry run and the tests.

### Benchmarks on real data

Unit tests do not show what a change does to a real tree. Changes that can move ages or placements are measured on the Iberian run, and the evidence is committed next to the data: `data/iberia/consistency/README.md` has a section per issue of the consistency layer, and `data/iberia/consistency/benchmark/` holds the release gate of [#56](https://github.com/naturalis/timegraft/issues/56): metrics declared in advance (`GATE.md`), a harness that runs the workflow scripts unchanged on perturbed copies of the inputs (injected OTT errors, noisy TimeTree ages, fewer genera), and the report. Rerun it, or its next version ([#74](https://github.com/naturalis/timegraft/issues/74)), before changing a default; see its README.

## Documentation

- Technical documentation in `docs/technical`, scientific in `docs/science`, the manuscript in `docs/manuscript`.
- Diagrams are Mermaid code blocks, which GitHub renders. Regenerate the workflow graph with `snakemake --rulegraph mermaid-js` after changing rules, and update `running.md`.

## Conventions

- All logic in `src/timegraft`, tested there; rule scripts stay thin.
- Every decision that changes an age or a placement is recorded in an output table (provenance, node ages, placements), so results can be traced. A change to an output format gets a migration note in `docs/technical/provenance.md`.
- An experiment never patches the pipeline: errors and noise go into copies of its inputs, as in the benchmark.
- Unit identifiers such as `BOLD:AAA0001` contain colons: never hand them to external tools as names; use neutral names and map back (as `build`, `placement` and `dating` do).
- UK English in documentation.
