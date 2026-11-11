# Crown ages from the literature

`literature.tsv` holds genus crown ages from published studies, for groups that TimeTree covers poorly ([#17](https://github.com/naturalis/timegraft/issues/17)). Give it to the workflow as `crown_ages` (as `config/iberia.yaml` does).

| Column | Meaning |
| --- | --- |
| `genus` | Genus name as in the records |
| `crown_age` | Point age in Ma. Used as given, in place of the age from the barcodes |
| `crown_age_min`, `crown_age_max` | Bounds in Ma on the crown age, for sources that give a limit and no estimate (a fossil gives a minimum). Only used when `crown_age` is empty |
| `reference`, `doi` | Where the age comes from |
| `note` | What the age refers to in the source, and anything a user should know |

Rules for adding a genus:

- Take the age of the node that joins the living species of the genus (the crown), not the split from its sister group (the stem).
- Record what the number is in the `note`: which node, which analysis, and whether an interval was available.
- One row per genus. Where studies disagree, choose one and say why in the `note`.
- An age that is not below the stem age of the genus in the backbone is cut to `cap_fraction` of the stem age (`crown_source` `table_capped` in `provenance.tsv`). That points to a backbone node that is too young, which is worth reporting ([#16](https://github.com/naturalis/timegraft/issues/16)).

## Status

Three genera so far, all checked against the cited source in November 2026. This is a start: which studies to trust is a scientific choice, and the table should be filled for the genera that matter most in the samples at hand.

| Genus | Age | In the Iberian run |
| --- | --- | --- |
| *Daphnia* | at least 145 Ma | Used as a minimum; the stem age is 179 Ma |
| *Gammarus* | 61 Ma | Above the stem age in the backbone (53 Ma), so cut to 0.9 of it |
| *Niphargus* | 87 Ma | Above the stem age in the backbone (23 Ma), so cut to 0.9 of it |

Searched for and not added, because no crown age could be taken from the sources at hand: *Alpheus*, *Proasellus*, *Helix*, *Littorina*, *Patella*, *Carcinus*, *Mytilus* and the lumbricid genera.
