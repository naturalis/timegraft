# timegraft documentation

The documentation has three parts.

## Technical: using and developing timegraft

| Document | What it covers |
| --- | --- |
| [Installation](technical/installation.md) | Python package, Snakemake, external tools, conda |
| [Configuration](technical/configuration.md) | The configuration file: inputs, marker plug-ins, every setting |
| [Running and outputs](technical/running.md) | Running the workflow, the Snakemake DAG, output files, run times |
| [Evidence and provenance](technical/provenance.md) | The evidence kept with the backbone, the TimeTree cache format, and a migration note for code that reads the outputs |
| [Development](technical/development.md) | Code layout, modules, plug-in points, extending, testing, conventions |
| [Roadmap](technical/roadmap.md) | Where the project is going, linked to the open issues |

## Scientific: what timegraft does and why

| Document | What it covers |
| --- | --- |
| [Background](science/background.md) | The problem, the users, and the earlier attempts at Naturalis |
| [Prior art](science/prior-art.md) | Related work elsewhere and how timegraft relates to it |
| [How it works](science/how-it-works.md) | The method explained without code, with flowcharts |
| [Assumptions and limitations](science/limitations.md) | What the method takes for granted, and where it is weak |
| [First results](science/results/2026-11-iberia.md) | The first run on real data: genera of western Iberia (numbers since superseded) |
| [Screening reference units](science/results/2026-11-screening.md) | How many units do not fit their genus, what screening does about it, and what it changes |
| [Crown ages: rates, bounds and penalised likelihood](science/results/2026-11-dating.md) | Why crown ages are handled as a range between bounds, and why penalised likelihood and fitted rates are not used |
| [Ages for undated backbone nodes](science/results/2026-11-interpolation.md) | Pure birth on clade sizes against even spacing, tested on held-out TimeTree ages and on subsampled genera |
| [Benchmark of the consistency layer](../data/iberia/consistency/benchmark/REPORT.md) | The release gate of the checks and corrections of epic #46: two full runs compared, injected errors and noise, and what each part contributes |

The results pages record what was found at the time; each says where later work changed its numbers. The current numbers are those of the benchmark. `data/iberia/consistency/README.md` records, per issue of the consistency layer, what it did on the Iberian data.

## Manuscript

[docs/manuscript](manuscript/README.md) holds the draft of a preprint describing timegraft, as Markdown with a BibTeX file. It is a separate activity from the documentation and may lag behind or run ahead of the code.
