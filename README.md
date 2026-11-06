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

Known crown ages (for instance from TimeTree) always take precedence. The provenance table records where every crown age came from.

## Status

Milestone 1, working on synthetic data:

- [x] Record filtering, one medoid representative per unit, genus assignment with conflict detection
- [x] Sister-genus outgroups from the backbone, BioNJ, rooting, mean-path-length dating, grafting
- [x] Sharded execution so genera are processed in parallel
- [x] Toy data set, unit tests and an end-to-end workflow test
- [ ] Backbone assembly: OToL induced subtree dated with TimeTree ages
- [ ] treePL as an alternative to mean-path-length dating
- [ ] Embedding distances as a third distance method
- [ ] Placement of units without a genus, and of unmatched metabarcoding sequences
- [ ] Per-sample pruning and PD metrics
- [ ] Validation of the informative window per marker and clade

## Usage

```bash
pip install -e ".[dev]"
snakemake --cores 4                 # runs the toy example in config/config.yaml
snakemake --cores 4 --use-conda     # same, with the conda environment per rule
pytest
```

### Inputs

- `records`: tab-separated, with columns `record_id`, `unit_id`, `genus`, `sequence`. For `align: none` the sequences must already be aligned; set `align: mafft` otherwise.
- `backbone`: Newick, ultrametric, ages in Ma, one tip per genus.
- `crown_ages` (optional): tab-separated, `genus`, `crown_age`.

### Outputs (under `outdir`)

- `timegraft.tre`: the dated tree, tips are unit identifiers.
- `provenance.tsv`: per backbone genus, whether it received a subtree, a single unit or was pruned, with stem age, crown age and its source, the rooting method and the outgroup used.
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
