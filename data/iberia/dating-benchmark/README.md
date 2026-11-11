# Dating benchmark: penalised likelihood, rates and bounds

Scripts and tables behind [docs/science/results/2026-11-dating.md](../../../docs/science/results/2026-11-dating.md) ([#36](https://github.com/naturalis/timegraft/issues/36)). They use the tolerant Iberian run as it was before crown ages became bounded (`provenance_before.tsv`; the run's other outputs are in the repository history, `data/iberia/iberia-run-tolerant.tar.gz` at commit `ea8edb3`).

| File | Contents |
| --- | --- |
| `export_trees.py` | Rebuilds the rooted, undated tree of every genus with the code path of `build_genus`, and writes `trees.tsv`, `taxonomy.tsv` and `outgroup_ages.tsv` |
| `trees.tsv.gz` | Per genus: units, stem age, rooting, mean root-to-tip distance (`mrt`), zero and negative branches, and the rooted tree in substitutions per site |
| `outgroup_ages.tsv.gz` | Per genus and outgroup unit: the backbone divergence age and the mean K2P distance to the genus |
| `taxonomy.tsv` | Family, order, class and phylum per genus |
| `chronos_bounds.R`, `A_bounds.tsv` | `chronos` with the crown bounded between half and twice the TimeTree crown age, three random starts, three models (105 genera) |
| `chronos_fixed.R`, `C_fixed.tsv`, `C_fixed_trees.tsv.gz`, `C_big.tsv` | `chronos` with the crown fixed, five settings, on 200 random genera of 4 to 60 units, and run times for three larger genera |
| `summarise.py` | The numbers quoted for the two `chronos` tests, including the comparison with mean path length |
| `holdout.py` | Crown age estimators against TimeTree crown ages, saturation of the distance to the outgroup, and the implied rate by crown age |
| `provenance_before.tsv` | `provenance.tsv` of the run before the change; its `crown_source` says which TimeTree crown ages the earlier rule accepted |

## Rerunning

With R (ape), MAFFT and timegraft installed, from the repository root:

```bash
B=data/iberia/dating-benchmark; W=$(mktemp -d)
mkdir -p $W/run && tar xzf <the archive at commit ea8edb3> -C $W/run
python $B/export_trees.py $W/run data/iberia/records.tsv.gz $W/out        # about 25 minutes on 2 cores
Rscript $B/chronos_bounds.R $W/out $B/provenance_before.tsv 300           # about 20 minutes
Rscript $B/chronos_fixed.R $W/out                                         # about 10 minutes
python $B/summarise.py $W/out
python $B/holdout.py $W/out $B/provenance_before.tsv
```

To check the numbers from the stored tables without rerunning, unpack `trees.tsv.gz` and `outgroup_ages.tsv.gz` into a folder with the other tables; `summarise.py` also needs the `newick` column of `C_fixed_trees.tsv.gz` joined to `C_fixed.tsv`.

`chronos` draws its starting ages at random, so the bounded test gives slightly different crown ages on each run; that is the point of the test.
