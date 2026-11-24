# Installation

## With conda (recommended)

`environment.yml` in the root of the repository defines one conda environment, named `timegraft`, with everything needed to run the workflow, the tests and the benchmark: Python, the package itself (installed in editable mode from your checkout, so changes after a `git pull` take effect without reinstalling), Snakemake, MAFFT, VSEARCH, curl and tar (for downloading the Open Tree of Life release files), pytest, matplotlib and Biopython.

Create it once, from the root of your checkout:

```bash
git clone https://github.com/naturalis/timegraft.git
cd timegraft
conda env create -f environment.yml
```

Activate it in every new shell before running anything:

```bash
conda activate timegraft
```

When a `git pull` changes `environment.yml`, bring the environment up to date (`--prune` removes packages the file no longer lists):

```bash
conda env update -n timegraft -f environment.yml --prune
```

To start again from scratch: `conda env remove -n timegraft`, then `conda env create -f environment.yml`. `mamba` accepts the same commands, and is faster at solving.

With the environment active, the workflow runs as it is (`snakemake --cores 4`); `--use-conda` is not needed, but still works (see [below](#external-tools)). Two optional dating methods are the exception and need `--use-conda`: `dating.method: lsd2`, whose `build_lsd2` rule compiles LSD2 with git, make and a C++ compiler, and `dating.method: chronos`, which needs R and ape. The default, `mpl`, needs neither.

## With pip

```bash
git clone https://github.com/naturalis/timegraft.git
cd timegraft
pip install -e ".[dev]"     # timegraft, Snakemake, pytest, matplotlib
```

Python 3.10 or later. The package depends on NumPy, DendroPy and DuckDB. This installs neither MAFFT nor VSEARCH: install them yourself, or run the workflow with `--use-conda`.

## External tools

| Tool | Needed for | How to get it |
| --- | --- | --- |
| MAFFT | `align: mafft` (real COI data, which is not pre-aligned) | `conda install -c bioconda mafft`, or `--use-conda` |
| VSEARCH | matching and placing ASVs; placing units without a genus | `conda install -c bioconda vsearch`, or `--use-conda` |
| LSD2 | `dating.method: lsd2` | built from source by the `build_lsd2` rule, in a conda environment with a C++ compiler |
| curl, tar | downloading the Open Tree of Life release files | usually present; `--use-conda` provides them |

The `timegraft` conda environment above has MAFFT and VSEARCH. Without it, let Snakemake create one conda environment per rule:

```bash
snakemake --cores 4 --use-conda
```

The environments are defined in `workflow/envs/`. LSD2 is not packaged for conda (the bioconda package called `lsd` is unrelated), so the workflow compiles a pinned commit once into `resources/bin/`. treePL is not packaged either and is not supported.

## macOS

Everything runs on macOS. Two things to watch:

- Shell examples in this documentation use only commands that behave the same in BSD and GNU versions; if you adapt them, note that `paste -sd, -` needs the trailing `-` on macOS.
- On Apple silicon, MAFFT comes from conda-forge (bioconda has no build for `osx-arm64`) and VSEARCH from bioconda; `environment.yml` lists conda-forge first, so `conda env create` finds both. Keep that channel order if you add packages.

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

The Iberian run (about 100,000 records, 25,000 units, 1,200 genera) takes about 50 minutes on two cores with a filled cache, and at most 1.3 GB of memory, in the backbone step ([benchmark](../../data/iberia/consistency/benchmark/REPORT.md)). Since [#72](https://github.com/naturalis/timegraft/issues/72) the backbone step reads `annotations.json` incrementally and peaks at 0.91 GB instead of 1.19 GB.

## Checking the installation

```bash
pytest                    # unit tests and an end-to-end run on the toy data
snakemake --lint          # workflow checks
```

The end-to-end test is skipped if Snakemake or VSEARCH is missing.
