# Using timegraft

## Installation

```bash
git clone https://github.com/naturalis/timegraft.git
cd timegraft
pip install -e ".[dev]"          # the package, Snakemake, pytest, matplotlib
```

External tools, needed depending on the settings:

| Tool | Needed for | Get it |
| --- | --- | --- |
| MAFFT | `align: mafft` (real COI data) | conda (`bioconda::mafft`), or `--use-conda` |
| VSEARCH | matching and placing ASVs, placing units without a genus | conda (`bioconda::vsearch`), or `--use-conda` |
| LSD2 | `dating.method: lsd2` | built by the `build_lsd2` rule (needs `--use-conda` or a C++ compiler) |

The simplest route is to let Snakemake create the environments: `snakemake --use-conda`. On macOS, check that the tools are available for your architecture.

## Running

```bash
snakemake --cores 4                                  # the toy example in config/config.yaml
snakemake --cores 4 --use-conda                      # with a conda environment per rule
snakemake --cores 8 --configfile config/config.yaml my_run.yaml   # your own settings on top
pytest                                               # unit tests and the toy workflow
```

A second config file only needs the keys that differ. For the first real run this was:

```yaml
records: /path/to/records.tsv
outdir: results/iberia
crown_ages:
samples:
asvs:
asv_table:
n_shards: 8
backbone_build:
  enabled: true
  context: Metazoa
  cache: resources/api_cache.json
validation:
  reference_tree:
markers:
  coi:
    align: mafft
```

## Inputs

### Reference records

Either a record table (`records`), tab-separated with at least `record_id`, `unit_id`, `genus`, `sequence` (optionally `family`, `order` and other ranks, used by the validation), or a BOLD data package (`bcdm.path`, TSV or Parquet), filtered by marker, taxa and country.

The same filter runs on its own, which is the practical way to cut a subset out of a full BOLD dump:

```bash
# 1. which genera occur in the region
python -m timegraft.bcdm BOLD.parquet regional.tsv \
    --taxon phylum=Annelida,Mollusca,Arthropoda,Echinodermata,Cnidaria,Bryozoa,Porifera \
    --exclude class=Insecta,Arachnida,Collembola,Diplopoda,Chilopoda \
    --country Portugal,Spain
cut -f3 regional.tsv | tail -n +2 | grep -v '^$' | sort -u > genera.txt

# 2. all BINs of those genera, worldwide, at most ten records per BIN
python -m timegraft.bcdm BOLD.parquet records.tsv --taxon genus=@genera.txt --max-per-unit 10
```

A list of genera from the samples themselves is a better filter than a country: a regional BOLD subset also contains terrestrial and freshwater taxa.

### Backbone

Either `backbone`: a dated, ultrametric Newick tree in Ma with one tip per genus; or `backbone_build.enabled: true`, which builds it from OToL release files (downloaded once by `download_otol`, about 140 MB) and TimeTree. Network access is needed to `files.opentreeoflife.org` and `api.timetree.org`. All TimeTree answers are cached in `backbone_build.cache`; keep that file between runs.

### Optional

- `crown_ages`: `genus`, `crown_age` (Ma). Always used as given. The place to put literature values for groups TimeTree covers poorly.
- `samples`: `sample`, `unit_id`, `abundance`.
- `asvs` (FASTA) and `asv_table` (`sample`, `asv_id`, `reads`), instead of `samples`.
- `validation.reference_tree`: a trusted tree whose tips are unit ids.

## Key settings

| Setting | Meaning |
| --- | --- |
| `marker` | which block under `markers` to use (`coi`, `its`) |
| `n_shards` | parallel jobs for the genus subtrees |
| `n_outgroup` | units from the sister genus used for rooting |
| `genus_agreement` | share of a unit's records that must agree on its genus |
| `prune_empty` | drop backbone genera without units |
| `place_unassigned_units` | place units without a usable genus |
| `markers.<m>.filter` | length, ambiguity, genetic code (COI) |
| `markers.<m>.align` | `none` or `mafft` |
| `markers.<m>.distance` | `k2p` (aligned) or `kmer` (with `k`) |
| `markers.<m>.crown_age` | `rate` (with `rate`) or `fraction`; `cap_fraction`; `timetree_max_ratio` |
| `markers.<m>.dating` | `mpl` or `lsd2` (with `seq_length`) |
| `markers.<m>.placement` | `match_identity`, `place_identity`, `query_cov`, `height` (`rate` or `terminal`) |

## Outputs (under `outdir`)

| File | Content |
| --- | --- |
| `timegraft.tre` | the dated tree; tips are unit ids (`grafted.tre` is the same before units without a genus were placed) |
| `provenance.tsv` | per backbone genus: subtree, single unit or pruned; stem and crown age; where the crown age came from (and the TimeTree value considered); rooting and dating method; outgroup |
| `backbone/genera.tsv` | how each genus was matched in OTT and placed in the synthesis |
| `backbone/node_ages.tsv` | per backbone node: age, its source, TimeTree's raw answer, the ids queried, and a conflict ratio |
| `backbone/crown_ages.tsv` | TimeTree crown age estimates, before checking |
| `units/rejected.tsv`, `units/unassigned.tsv`, `genera/missing_from_backbone.tsv`, `units_placed.tsv` | what did not make it into the tree directly, and why |
| `placement/placements.tsv`, `placement/timegraft_placed.tre` | per ASV: matched, placed or unplaced; the tree with placed ASVs |
| `diversity/alpha.tsv`, `diversity/beta.tsv`, `diversity/units_not_in_tree.tsv` | diversity per sample and per pair of samples |
| `validation/` | pairs, correlations per stratum, binned medians, plateau onset, `saturation.png` |

## Run time

The first real run (102,794 records, 1,232 genera, 25,578 BINs; the largest genus has 1,035 BINs) took 20 to 28 minutes on 2 cores: about 5 minutes for representatives, 10 for the backbone (mostly TimeTree, fast on later runs thanks to the cache) and 13 for the genus subtrees with MAFFT.
