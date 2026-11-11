# Substitution rates calibrated on the Iberian genera

`rates.tsv` holds the COI rates used by `config/iberia.yaml` ([#18](https://github.com/naturalis/timegraft/issues/18)): substitutions per site per million years, per lineage, under the Kimura two-parameter model.

| Column | Meaning |
| --- | --- |
| `rank`, `name` | The taxon; the row `all` is the pooled value, used as `crown_age.rate` |
| `rate` | Trimmed median of the implied rates of the taxon's genera |
| `n`, `n_trimmed` | Genera used, and genera left out for being more than three times above or below the median |
| `q1`, `q3` | Quartiles of the implied rates used |
| `untrimmed_median` | Median before trimming |

For a genus with a known crown age, the implied rate is the mean root-to-tip distance of its rooted subtree over that age. Known ages are point ages in `data/crown_ages/literature.tsv` (none qualified here, as both are above the stem age) and TimeTree crown ages below the stem age: 262 genera. Rates are given per class with at least ten genera after trimming.

To rebuild from the archived run (`../iberia-run-tolerant.tar.gz`, unpacked into `results/iberia`):

```bash
python -m timegraft.calibrate --provenance results/iberia/provenance.tsv \
    --taxonomy results/iberia/genera/taxonomy.tsv \
    --crown-ages data/crown_ages/literature.tsv --out data/iberia/calibration/rates.tsv
```

## What to keep in mind

- **The calibration set is narrow.** 128 of the 262 genera are teleosts. Classes without a rate of their own (all annelids, bivalves, echinoderms, cnidarians and sponges) get the pooled value, which fishes dominate.
- **The classes differ.** Elasmobranchs evolve about three times more slowly than teleosts, in line with what is known of shark mitochondrial DNA.
- **A constant rate is a simplification.** In these genera the mean root-to-tip distance grows far less than in proportion to crown age, so the implied rate falls from 0.016 for crown ages below 5 Ma to 0.001 above 80 Ma (`../dating-benchmark/`). The rates here are those of a typical genus, with a crown age of 10 to 40 Ma. Young genera will come out too old and old genera too young.
- **TimeTree crown ages are themselves uncertain**, and inflated for genera that are not natural groups. Trimming removes the worst of these; the untrimmed medians are given for comparison.
