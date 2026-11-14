# Ages for backbone nodes that TimeTree does not date (November 2026)

TimeTree dates only part of the backbone: in the Iberian run 261 of 779 nodes directly and 53 through enclosing taxa, leaving 465 nodes without an age. These were spaced evenly between their nearest dated ancestor and their oldest dated descendant, as in BLADJ. [#16](https://github.com/naturalis/timegraft/issues/16) asked for something better. Literature ages will never cover all of these nodes (collecting them is what TimeTree already does), and barcodes cannot help above genus level, least of all k-mer distances for ITS, so some model is needed. This page reports the tests of a pure-birth model against even spacing and what was changed. Scripts and tables are in [`data/iberia/backbone-interpolation`](../../../data/iberia/backbone-interpolation/README.md).

In short: a pure-birth model is as accurate as even spacing only if it is given the size of the whole clade, not the number of genera that happen to be in the data. Given the number of genera in the clade in the Open Tree of Life synthesis, it is about as accurate as even spacing, and much less sensitive to which genera are sampled. It is now the default.


> **Since then** ([#52](https://github.com/naturalis/timegraft/issues/52)): conflicting TimeTree ages are no longer resolved by precedence (dropping an age not younger than its dated ancestor) but fitted together, weighted by their evidence; interpolation of undated nodes is unchanged and comes after the fit. In the current Iberian run 387 of 784 nodes have a TimeTree age and 397 are interpolated ([benchmark](../../../data/iberia/consistency/benchmark/REPORT.md)).
## The model

Under a pure-birth (Yule) process with speciation rate λ, a clade of *n* lineages is expected to have its crown at (*H<sub>n</sub>* − 1)/λ before the present, where *H<sub>n</sub>* is the *n*-th harmonic number, so roughly at ln *n*. An undated node with *n* lineages, below a dated node with *n<sub>a</sub>* and above a dated node (or a tip) with *n<sub>f</sub>*, is placed at the fraction (*h*(*n*) − *h*(*n<sub>f</sub>*)) / (*h*(*n<sub>a</sub>*) − *h*(*n<sub>f</sub>*)) of the way up from the lower age, with *h*(*n*) = *H<sub>n</sub>* − 1 and *h* = 0 for a tip. This amounts to fitting λ to each pair of dated nodes, so every TimeTree age is kept and ages still decrease towards the tips. Nothing in it depends on the barcodes, so it is the same for COI and ITS.

What *n* should count is the question. Three options were tested:

- **Genera in the backbone** (`pb_sampled`): the tips below the node in our tree.
- **Genera in the clade** (`pb_genera`): the genera of all tips below the corresponding node of the synthetic tree (opentree16.1 with OTT 3.7.3), whether or not they are in the data.
- **Tips of the clade** (`pb_species`): all tips below that node of the synthetic tree, mostly species.

## Tests

### 1. Held-out TimeTree ages

Each node dated directly by TimeTree was hidden, alone (leave one out, 261 nodes) and in random halves (50 replicates, 6,500 predictions; this leaves longer undated chains, closer to the poorly covered parts of the tree), and interpolated again from the remaining ages. Bias is the median ratio of interpolated to TimeTree age; error is the median of the larger over the smaller of the two.

| Method | Bias, leave one out | Bias, halves | Error, leave one out | Error, halves | Within a factor 2, halves | Bias on two-genus nodes, halves |
| --- | --- | --- | --- | --- | --- | --- |
| Even spacing | 1.00 | 1.00 | 1.21 | 1.26 | 88% | 0.88 |
| Pure birth, genera in the backbone | 0.97 | 0.99 | 1.26 | 1.29 | 84% | 0.72 |
| Pure birth, genera in the clade | 1.01 | 1.03 | 1.27 | 1.30 | 83% | 0.92 |
| Pure birth, tips of the clade | 1.06 | 1.15 | 1.27 | 1.33 | 84% | 1.39 |

- **Counting only the genera in the data makes shallow nodes far too young.** A pair of genera in the data usually stands for a much larger clade, and a two-genus node came out at 0.72 of its TimeTree age.
- **Counting species makes them too old** (1.39), the opposite error: ln *n* grows too slowly with clade age, in line with the well-known weak relation between clade age and species richness, and species counts in the taxonomy reflect description effort as much as diversification.
- **Counting genera in the clade is close to unbiased** at every depth (0.92 on two-genus nodes, 1.04 to 1.10 at other sizes), but scatters a little more than even spacing.
- **Even spacing is not worse than expected.** It made held-out shallow nodes slightly too young (0.88), not too old as the issue supposed; the 291 Ma for *Hediste* and *Nereis* came from a small test with few dated nodes. Where TimeTree dates most of the neighbours, a node is usually halfway between two dated nodes, and the midpoint is a good guess.

These nodes are those TimeTree can date, which leans towards vertebrates and other well-studied groups, while the nodes that need interpolation lie mostly among invertebrates. The test measures the methods where the truth is known, not where they are used.

### 2. Dependence on which genera are sampled

Even spacing places a node by how many undated nodes lie between the dated ones, and that number depends on which genera are in the data. The same clade can then get a different age in a run on another region or another marker. To measure this, a random share of the Iberian genera was kept (50 replicates per share), the backbone pruned to them with the TimeTree ages of the surviving nodes, and the undated nodes interpolated again. The table gives the 90th percentile of the factor by which a node's age moved against the full backbone.

| Method | 75% of genera kept | 50% kept | 25% kept |
| --- | --- | --- | --- |
| Even spacing | 1.34 | 1.59 | 1.78 |
| Pure birth, genera in the backbone | 1.45 | 1.77 | 2.11 |
| Pure birth, genera in the clade | 1.01 | 1.20 | 1.34 |
| Pure birth, tips of the clade | 1.04 | 1.28 | 1.52 |

With pure birth on the genera in the clade, a node's age is a property of the clade: it moves only when a dated node nearby drops out. The median node does not move at all (factor 1.00, 1.00 and 1.01), against 1.00, 1.05 and 1.14 with even spacing.

## What was changed

- Undated backbone nodes are placed by the pure-birth model with the number of genera in the clade in the synthetic tree (`backbone_build.interpolation: pure_birth`, the default). The count is recorded per node in `node_ages.tsv` (`clade_genera`).
- Even spacing remains available (`interpolation: even`), and is used when no clade sizes are available, as with `backbone_build.source: api`, since pure birth on the genera in the data alone is the worst option tested.
- Where the clade sizes do not increase strictly towards the root (which the synthetic tree should not allow), the node is spaced evenly.

The choice rests on the second test, not the first: the two methods are about equally accurate, but only pure birth on the clade gives the same age to the same clade whatever else is in the data, which matters when results from different regions or markers are compared.

## Effect on the Iberian backbone

With the run's TimeTree ages kept (`effect.py`):

| | Even spacing | Pure birth, genera in the clade |
| --- | --- | --- |
| Undated nodes, age relative to even spacing, median (quartiles) | 1 | 1.06 (0.93 to 1.37) |
| Median stem age of the 788 genera with an interpolated stem | 91.1 Ma | 99.5 Ma |
| Total branch length of the backbone | 181,573 Ma | 192,540 Ma (+6%) |
| Stems younger than the genus crown age from TimeTree, of 126 | 23 | 19 |

Ages change in both directions, depending on the size of the clade rather than the length of the chain: *Nereis* and *Platynereis* go from 95 to 39 Ma, *Hediste* from 190 to 267 Ma and *Marphysa* from 82 to 205 Ma. No independent ages for these nodes have been checked yet.

## What this does not solve

- **The model is crude.** Pure birth ignores extinction and differences in rate between clades, and counts of genera depend on how finely each group has been split into genera. It is used because something has to fill the gaps, not because it is right.
- **Literature ages are still better where they exist.** A curated table of node ages for groups TimeTree covers poorly ([#17](https://github.com/naturalis/timegraft/issues/17)) would take precedence, and would also give independent ages to test this model on invertebrates.
- **Some TimeTree ages conflict with the literature**, which is not a matter of interpolation: the stems of *Gammarus* (53 Ma) and *Niphargus* (23 Ma) are dated by TimeTree, below their published crown ages of 61 and 87 Ma.
