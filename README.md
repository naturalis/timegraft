# timegraft

[![CI](https://github.com/naturalis/timegraft/actions/workflows/ci.yml/badge.svg)](https://github.com/naturalis/timegraft/actions/workflows/ci.yml)
[![Licence: Apache-2.0](https://img.shields.io/badge/licence-Apache--2.0-blue.svg)](LICENSE)
[![Output data: CC0](https://img.shields.io/badge/output%20data-CC0%201.0-lightgrey.svg)](https://creativecommons.org/publicdomain/zero/1.0/)
[![Python ≥ 3.10](https://img.shields.io/badge/python-%E2%89%A5%203.10-3776ab.svg?logo=python&logoColor=white)](pyproject.toml)
[![Snakemake ≥ 8](https://img.shields.io/badge/snakemake-%E2%89%A5%208-039475.svg)](https://snakemake.readthedocs.io)
[![Open issues](https://img.shields.io/github/issues/naturalis/timegraft.svg)](https://github.com/naturalis/timegraft/issues)
[![Last commit](https://img.shields.io/github/last-commit/naturalis/timegraft.svg)](https://github.com/naturalis/timegraft/commits/main)

**Dated phylogenies of DNA barcode reference libraries, for phylogenetic diversity of metabarcoding samples.**

timegraft builds a time-calibrated tree with barcode reference units at the tips (BINs for COI, UNITE species hypotheses for ITS), places metabarcoding sequences on it, and computes phylogenetic diversity per sample. Deep relationships and ages come from Open Tree of Life and TimeTree; barcodes are used only within genera, where they carry reliable signal. Branch lengths are in millions of years, so diversity is comparable across sites, studies and markers.

[**Documentation**](docs/README.md) · [**Quick start**](#quick-start) · [**How it works**](docs/science/how-it-works.md) · [**Roadmap**](docs/technical/roadmap.md) · [**Issues**](https://github.com/naturalis/timegraft/issues)

## Features

- Reads BOLD data packages (TSV or Parquet) or plain record tables
- Builds a dated genus-level backbone from Open Tree of Life release files and TimeTree
- One dated subtree per genus, in parallel, then grafted onto the backbone without any global tree search
- Marker plug-ins for filtering, distances, crown ages and dating; COI and ITS settings included
- Matches metabarcoding ASVs to reference units and places unmatched ones on the tree
- Faith's PD, MPD, UniFrac and PhyloSor per sample and per pair of samples
- Validation of barcode distances against a reference tree, per taxonomic stratum
- Every age and placement traceable to its source in output tables

## Quick start

```bash
git clone https://github.com/naturalis/timegraft.git
cd timegraft
pip install -e ".[dev]"
snakemake --cores 4 --use-conda     # runs the toy example in config/config.yaml
```

The dated tree is written to `results/toy/timegraft.tre`. To run on your own data, see [running timegraft](docs/technical/running.md).

## Documentation

| | |
| --- | --- |
| **Technical** | [Installation](docs/technical/installation.md) · [Configuration](docs/technical/configuration.md) · [Running and outputs](docs/technical/running.md) · [Development](docs/technical/development.md) · [Roadmap](docs/technical/roadmap.md) |
| **Scientific** | [Background](docs/science/background.md) · [Prior art](docs/science/prior-art.md) · [How it works](docs/science/how-it-works.md) · [Assumptions and limitations](docs/science/limitations.md) · [First results](docs/science/results/2026-11-iberia.md) |
| **Manuscript** | [Draft preprint](docs/manuscript/manuscript.md) |

## Citing

A preprint is in preparation ([draft](docs/manuscript/manuscript.md)). Until then, please cite this repository.

## Licence

Code: [Apache-2.0](LICENSE). Trees and tables produced by timegraft: [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).

Developed at the [Naturalis Biodiversity Center](https://www.naturalis.nl). The design builds on [Bactria](https://github.com/naturalis/barcode-constrained-phylogeny) and [evospaice](https://github.com/naturalis/evospaice).
