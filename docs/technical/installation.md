# Installation

## The package and Snakemake

```bash
git clone https://github.com/naturalis/timegraft.git
cd timegraft
pip install -e ".[dev]"     # timegraft, Snakemake, pytest, matplotlib
```

Python 3.10 or later. The package depends on NumPy, DendroPy and DuckDB.

## External tools

| Tool | Needed for | How to get it |
| --- | --- | --- |
| MAFFT | `align: mafft` (real COI data, which is not pre-aligned) | `conda install -c bioconda mafft`, or `--use-conda` |
| VSEARCH | matching and placing ASVs; placing units without a genus | `conda install -c bioconda vsearch`, or `--use-conda` |
| LSD2 | `dating.method: lsd2` | built from source by the `build_lsd2` rule, in a conda environment with a C++ compiler |
| curl, tar | downloading the Open Tree of Life release files | usually present; `--use-conda` provides them |

The simplest route is to let Snakemake create one conda environment per rule:

```bash
snakemake --cores 4 --use-conda
```

The environments are defined in `workflow/envs/`. LSD2 is not packaged for conda (the bioconda package called `lsd` is unrelated), so the workflow compiles a pinned commit once into `resources/bin/`. treePL is not packaged either and is not supported.

## macOS

Everything runs on macOS. Two things to watch:

- Shell examples in this documentation use only commands that behave the same in BSD and GNU versions; if you adapt them, note that `paste -sd, -` needs the trailing `-` on macOS.
- On Apple silicon, check that MAFFT and VSEARCH from conda are available for your architecture (`conda install` will say if not).

## Network access

Needed only when building the backbone (`backbone_build.enabled: true`), and not at all when the cache already holds every answer (`data/iberia/timetree_cache.json` does for the Iberian run):

| Host | For |
| --- | --- |
| `files.opentreeoflife.org` | Open Tree of Life release files, once: the synthesis tree (31 MB), its `annotations.json` (61 MB; only with `backbone_build.synthesis_support`) and the taxonomy (about 420 MB unpacked) |
| `api.timetree.org` | TimeTree node ages and the evidence behind them (answers are cached on disk) |
| `api.checklistbank.org` | Catalogue of Life names, classifications and identifiers, with `backbone_build.col.enabled` (the default; answers are cached in the same file) |
| `api.opentreeoflife.org` | only with `backbone_build.source: api` |
| `github.com` | building LSD2 |

## Memory and time

The Iberian run (about 100,000 records, 25,000 units, 1,200 genera) takes about 50 minutes on two cores with a filled cache, and at most 1.3 GB of memory, in the backbone step ([benchmark](../../data/iberia/consistency/benchmark/REPORT.md)).

## Checking the installation

```bash
pytest                    # unit tests and an end-to-end run on the toy data
snakemake --lint          # workflow checks
```

The end-to-end test is skipped if Snakemake or VSEARCH is missing.
