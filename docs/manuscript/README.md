# Manuscript

A draft preprint describing timegraft. This is a separate activity from the documentation: it describes the method as it stands when the manuscript is submitted, which may differ from the code at any given time.

| File | Content |
| --- | --- |
| `manuscript.md` | The text, in Pandoc Markdown with citations as `[@key]` |
| `references.bib` | BibTeX references |

## Building

With [Pandoc](https://pandoc.org) 2.11 or later:

```bash
cd docs/manuscript
pandoc manuscript.md --citeproc -o manuscript.pdf     # needs a LaTeX engine
pandoc manuscript.md --citeproc -o manuscript.docx    # for co-authors who prefer Word
```

A journal style can be added with `--csl` and a CSL file from the [Zotero style repository](https://www.zotero.org/styles).

## Status

- **Draft**: introduction, methods and the first results are written; discussion is in outline; figures and tables are placeholders.
- **Authors**: to be agreed.
- **References**: entries with a `note` starting "Check:" have details that still need verifying against the source.
- **Numbers** in the results come from the first run on real data ([docs/science/results/2026-11-iberia.md](../science/results/2026-11-iberia.md)) and will change as the method is calibrated (see the [roadmap](../technical/roadmap.md)).

Open points are marked `[TODO]` in the text.
