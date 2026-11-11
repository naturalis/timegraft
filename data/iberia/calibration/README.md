# Substitution rates calibrated on the Iberian genera (not used)

`rates.tsv` holds COI rates per class derived from the Iberian genera with `python -m timegraft.calibrate` ([#18](https://github.com/naturalis/timegraft/issues/18)): substitutions per site per million years, per lineage, under the Kimura two-parameter model. **`config/iberia.yaml` does not use them.** They are kept as a record of an approach that did not work, and why.

| Column | Meaning |
| --- | --- |
| `rank`, `name` | The taxon; the row `all` is the pooled value |
| `rate` | Trimmed median of the implied rates of the taxon's genera |
| `n`, `n_trimmed` | Genera used, and genera left out for being more than three times above or below the median |
| `q1`, `q3` | Quartiles of the implied rates used |
| `untrimmed_median` | Median before trimming |

For a genus with a known crown age, the implied rate is the mean root-to-tip distance of its rooted subtree over that age. Known ages here are the 262 TimeTree crown ages below the stem age.

## Why they are not used

- **They are not rates.** The mean root-to-tip distance hardly grows with crown age: its median is 0.04 for crown ages below 5 Ma and 0.13 above 80 Ma, so the implied rate falls from 0.016 to 0.001 (`../dating-benchmark/`). A median over genera then gives the distance of a typical genus over the crown age of a typical genus, about 0.004, which is two to seven times slower than COI rates calibrated on splits of a few million years (0.007 to 0.026 per lineage).
- **They push every genus to the same age.** With these rates the median crown age of the Iberian genera went from 9 to 25 Ma, and 129 genera reached 0.9 of their stem age (5 before), as often where the stem age comes from TimeTree as where it is interpolated.
- **The calibration set is narrow.** 128 of the 262 genera are teleosts. Marine invertebrates have few fossils and will keep having few calibration points.

One result is worth keeping: elasmobranchs come out about three times slower than teleosts, in line with what is known of shark mitochondrial DNA. A slower rate for a class can be set in `crown_age.rates`.

To rebuild the table from the archived run (`../iberia-run-tolerant.tar.gz`, unpacked into `results/iberia`):

```bash
python -m timegraft.calibrate --provenance results/iberia/provenance.tsv \
    --taxonomy results/iberia/genera/taxonomy.tsv \
    --crown-ages data/crown_ages/literature.tsv --out data/iberia/calibration/rates.tsv
```
