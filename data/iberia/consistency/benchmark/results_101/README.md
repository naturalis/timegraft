# Species found by name inside the target (#101): the injection part of the gate

[#70](https://github.com/naturalis/timegraft/issues/70) attaches a re-anchored genus within its target taxon, on the species of its matched OTT taxon that the synthesis places inside the target. Where OTT gave the genus's name to another taxon without its species, the matched taxon holds none of them, and the genus went to the target's node: in the `inject_name` variant of the injections, none of the 120 genera was placed within its target (`../results_70/README.md`, "What this does not show"). [#101](https://github.com/naturalis/timegraft/issues/101) looks for the genus's species by name inside the target in that case, and dates the genus with their NCBI ids.

## Method

As for #70 (`../results_70/README.md`): only the default arm of the three injections (`inject`, `inject_own`, `inject_name`), five replicates each, was rebuilt, from the cache of its build in #73 and with the same injected taxonomy; everything else comes from the benchmark of #70 and #73. The rebuild was made twice:

1. with the placement by name alone (the genus still dated with the ids of its matched taxon), column `placement_by_name`;
2. with the placement by name and the NCBI ids of the species found by name, as merged, column `with_101`.

```
PYTHONPATH=$B:src python $B/backbone.py --repo . --units $BENCH73/branch/results/iberia/units \
    --otol $BENCH73/branch/resources/otol --cache <copy of exp/K/sN/default/cache.json of #73> \
    --out exp/K/sN/default --taxonomy $BENCH73/exp/K/sN/taxonomy.tsv
PYTHONPATH=$B:src python $B/evaluate.py --version 2 $BENCH101 out
```

with `B=data/iberia/consistency/benchmark`. The second rebuild asked TimeTree for the nodes whose query ids changed; those 401 answers were added to `../timetree_cache.json.gz` with `merge_caches.py`, and a rebuild of `inject_name` seed 3 from it with `--offline` gives the same backbone, node ages, crown ages and `genera.tsv`. Only metrics 1a to 1m are reported, for the reasons given in `../results_70/README.md`.

The Iberian backbone was built three times each with `main`, with this branch, and with this branch and `--set backbone_build.reanchor.species_by_name=false`, offline from the cache of #73's clean build, on the full OTT 3.7.3 taxonomy, measuring the build time and the peak memory of the process (`iberia_runtime.tsv`).

## Result

`injection_metrics.tsv`:

| Metric | With #70 | Placement by name | With #101 |
| --- | --- | --- | --- |
| 1f corrected genera with a stem within 1.5 | 119 of 120 | 120 of 120 | 120 of 120 |
| 1g other stems moved by more than 10% (legacy 248) | 4 | 3 | 3 |
| 1j `inject_own`: stem within 1.5 (other stems moved) | 119 of 120 (4) | 120 of 120 (3) | 120 of 120 (3) |
| **1m `inject_name`: stem within 1.5 (report)** | **86 of 120** | **106 of 120** | **116 of 120** |
| `inject_name`: other stems moved (legacy 184) | 98 | 86 | **16** |

Every other injection metric is unchanged and passes.

`within_target.tsv` says where each moved genus went. In `inject_name`, 114 of the 120 injected genera are now placed within their target, all on species found by name (93 at the majority clade, 16 at the MRCA, 5 at the parent of one); 3 stay at the target's node, because the synthesis has none of their species inside it, and 3 at the enclosing node, because the target is not in the synthesis. *Onchidella* (the five "not injected" rows) stays at the family node.

**Placement alone is not enough.** After the first rebuild, 14 `inject_name` genera were placed within their target and still missed 1.5. In all 14, every NCBI id was dropped by the identifier check of #69: in the name-only case the matched taxon holds the other animal's ids, so TimeTree had no id for the genus, and its stem was interpolated. The species found by name inside the target carry the genus's own NCBI ids, and dating with them (still through the check of #69) brings 1m to 116 and the other stems moved from 86 to 16. Of the four that still miss, *Spicara* has none of its species inside Sparidae in the synthesis and stays at the family node; *Myliobatis*, *Hygophum* and *Hemicycla* are placed within their target and keep 11 to 22 ids, with stems 1.56 to 1.64 times those without the error: they are dated with species ids only, as the genus id of the matched taxon belongs to the other animal.

**`inject` and `inject_own` change in one genus.** *Camachoaglaja* (seed 4, the one genus that missed 1f and 1j with #70, by a factor 2.3) went to the target's node with #70 and now goes to the MRCA of two species found by name inside Aglajidae. In OTT 3.7.3 itself its genus taxon is barren: its two species, *C. mariagordae* and *C. binter*, hang under *Chelidonura*, next to it (flag `sibling_higher`). So the merge moved no species into the host, and #70 found none of them inside the target. This is the name search doing what it is for, on a real gap in OTT, and it accounts for the whole change in 1f, 1g and 1j. The issue expected these two variants to be unchanged.

**The Iberian run is unchanged.** `backbone.tre`, `node_ages.tsv`, `crown_ages.tsv`, `genera.tsv` and `concordance.tsv` are byte-identical between `main`, this branch, and this branch without the name search. *Onchidella* has a single OTT taxon, the merged one, so no species named *Onchidella* lies outside it, and it stays at the family node (`reanchor_how` ends in `target_node`).

**Time and memory.** `iberia_runtime.tsv`: the build takes 25 to 28 s with the name search, 24 to 28 s without and 25 to 26 s on `main`, with a peak of 908 to 909 MB in every run. The pass over the taxonomy alone takes 11.4 and 12.3 s with the index and 12.4 and 11.3 s without: the difference is within the variation between runs. The index holds 2,321 species, for 139 of the 1,218 genera of the run (species named after a genus of the run but outside its matched taxon, mostly homonyms elsewhere).

## What this does not show

- **Name-only cases in real data.** Like the merge injections, `inject_name` changes only the taxonomy, so the synthesis still places the genus's species in their own family, where the name search finds them. A real name-only error in OTT might also be followed by the synthesis.
- **Synonyms.** Species under the genus's name only are found. Where OTT holds the genus's species under a synonym of the genus (a different genus name), they are not, and the genus goes to the target's node as before.
