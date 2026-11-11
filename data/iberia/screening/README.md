# Screening units: analysis behind #38

Scripts and tables behind [docs/science/results/2026-11-screening.md](../../../docs/science/results/2026-11-screening.md) ([#38](https://github.com/naturalis/timegraft/issues/38)).

| File | Contents |
| --- | --- |
| `unit_matches.py`, `unit_matches.tsv.gz` | Exploratory, before the screening step existed: the nearest unit of the same genus and of another genus for every representative, and the rank the two genera share |
| `prune_impact.py` | Exploratory: within-genus distances after pruning suspect units from the rooted genus trees of `../dating-benchmark/trees.tsv.gz` |
| `compare_runs.py` | Two runs that differ in screening: decisions, placements, distances, crown ages, agreement with TimeTree |
| `sample_pd.py` | Faith's PD of the same random samples on two trees of the same units |

The run before screening is the archive at the commit before this change; the run after is `../iberia-run-tolerant.tar.gz`. With both unpacked:

```bash
python data/iberia/screening/compare_runs.py BEFORE/provenance.tsv AFTER/provenance.tsv \
    AFTER/units_placed.tsv AFTER/units/screened.tsv
python data/iberia/screening/sample_pd.py BEFORE/timegraft.tre AFTER/timegraft.tre AFTER/units/screened.tsv
```

The all-against-all search for `unit_matches.py` (two minutes on two cores):

```bash
vsearch --usearch_global reps.fa --db reps.fa --self --id 0.80 --maxaccepts 10 --maxrejects 64 \
    --query_cov 0.7 --blast6out hits.tsv --threads 2
```
