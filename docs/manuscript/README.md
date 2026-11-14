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
- **Numbers** in the results come from the archived Iberian run of 12 November 2026 ([data/iberia](../../data/iberia/README.md)), which predates the consistency layer.
- **Out of date since the consistency layer** ([#46](https://github.com/naturalis/timegraft/issues/46), merged November 2026):
  - *Methods, Reference units*: screening now searches one-sided hit lists again and records codon diagnostics.
  - *Methods, Dated backbone*: COL as identifier hub; majority-clade anchoring instead of the MRCA; placement checks and re-anchoring; a weighted fit of ages instead of the precedence rule; re-queries.
  - *Methods, Genus subtrees*: TimeTree crowns of the sampled species, the plausibility test; long branches listed.
  - *Results*: the counts of set-aside units, backbone node sources, crown sources and rooting, and the *Onchidella* paragraph (now corrected by the pipeline).
  - A new reference run (#77) should be the source of the revised numbers; the benchmark ([report](../../data/iberia/consistency/benchmark/REPORT.md)) could become an evaluation section.

Open points are marked `[TODO]` in the text.

## Figures

Each figure in `figures/` is drawn by a script next to it. For Figure 2:

```
cd figures
python fig2.py                                    # writes fig2.svg
rsvg-convert -f pdf -o fig2.pdf fig2.svg
rsvg-convert -d 300 -p 300 -o fig2.png fig2.svg
```
