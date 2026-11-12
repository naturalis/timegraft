# Backbone interpolation: pure birth against even spacing

Scripts and tables behind [docs/science/results/2026-11-interpolation.md](../../../docs/science/results/2026-11-interpolation.md) ([#16](https://github.com/naturalis/timegraft/issues/16)). They use the backbone of the Iberian run in `data/iberia/iberia-run.tar.gz` (`backbone/backbone.tre`, `node_ages.tsv`, `genera.tsv` and `crown_ages.tsv`), whose undated nodes were spaced evenly.

| File | Contents |
| --- | --- |
| `richness.py` | Induces the backbone topology again from the Open Tree of Life files, checks that it is the run's, and writes `richness.tsv` |
| `richness.tsv` | Per backbone node (named as in `node_ages.tsv`): genera below it in the backbone, genera and tips below the corresponding node of the synthetic tree (opentree16.1, OTT 3.7.3, context Metazoa) |
| `holdout.py` | Hides TimeTree ages (each alone, and in 50 random halves) and compares the interpolated ages with them, per method, clade size and age |
| `stability.py` | Keeps 75%, 50% and 25% of the genera at random (50 times each) and measures how far the ages of undated nodes move |
| `effect.py` | Ages of the undated nodes, genus stems and total branch length with each method, and stems younger than the TimeTree crown age of their genus |

## Rerunning

With timegraft installed, from the repository root:

```bash
B=data/iberia/backbone-interpolation; W=$(mktemp -d)
mkdir -p $W/run && tar xzf data/iberia/iberia-run.tar.gz -C $W/run backbone
# only to recompute richness.tsv: the release files as the workflow downloads them (about 450 MB)
mkdir -p $W/otol/opentree16.1 $W/otol/ott3.7.3
curl -sSfL -o $W/otol/opentree16.1/labelled_supertree.tre \
  https://files.opentreeoflife.org/synthesis/opentree16.1/opentree16.1_tree/labelled_supertree/labelled_supertree.tre
curl -sSfL https://files.opentreeoflife.org/ott/ott3.7.3/ott3.7.3.tgz | tar xzf - -C $W/otol ott3.7.3/taxonomy.tsv
python $B/richness.py $W/run $W/otol $B/richness.tsv     # about 20 seconds
python $B/holdout.py $W/run $B/richness.tsv              # about 10 seconds
python $B/stability.py $W/run $B/richness.tsv            # about 30 seconds
python $B/effect.py $W/run $B/richness.tsv
```

The random draws are seeded, so the numbers are reproduced exactly.
