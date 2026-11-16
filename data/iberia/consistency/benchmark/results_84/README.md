# Why the stem ages of #56 changed (#84)

Metric 6b of the gate of [#56](https://github.com/naturalis/timegraft/issues/56) asked what share of the stem ages that differ by more than 10% between the legacy and the default run has a recorded cause. It passed with 175 of 175, but only by counting any fitted age on the path from the stem to the root as a cause; on a stricter reading 159 of 175 were explained. [#84](https://github.com/naturalis/timegraft/issues/84) records per backbone node where its final age comes from (`age_source`, `interpolated_from` and `query_genera` in `node_ages.tsv`, format 5), and `compare_backbones.py` uses these to follow each changed stem down to the nodes whose ages actually moved it.

## Method

Both backbones were rebuilt on `issue-84-node-causes` with the harness, offline, from the cache of the #73 gate rerun:

```
PYTHONPATH=$B python3 $B/backbone.py --repo . --units data/iberia/units --otol $OTOL --cache api_cache.json \
    --out legacy --offline --overlay $B/legacy.yaml
PYTHONPATH=$B python3 $B/backbone.py --repo . --units data/iberia/units --otol $OTOL --cache api_cache.json \
    --out default --offline
python3 $B/compare_backbones.py legacy default results_84 --factor 1.10
```

with `B=data/iberia/consistency/benchmark`. Both backbones are identical to the gate's `legacy_clean` and `default_clean`: the new columns change no age.

For every genus whose stem age differs by more than the factor, the script starts at the stem node of the second run and compares it with the node with the same genera in the first. A calibrated node changes because its query changed (an `identifier` of one of its `query_genera`, other `picks`, or a query on enclosing taxa, `higher_taxa`), because TimeTree answered differently (`timetree_answer`), or because the fit moved it (`fit`). A node calibrated in one run only adds `calibration_gained` or `calibration_lost`, which says what happened but not why, so the script also names the reason the other run has no calibration: TimeTree had no age for any of its picks (`no_answer`, or `identifier` if a genus of the node was sent with another id), the runs sent other queries (`identifier`, `picks`, `higher_taxa`), TimeTree answered the same query differently (`timetree_answer`), TimeTree's age was that of an ancestor because its MRCA lies above the node (`unresolved`), or the precedence rule dropped the age that the fit keeps (`fit`). A stem whose only cause is a gained or lost calibration counts as not traced. An interpolated node moves with the anchors it was placed between (`interpolated_from`), which the script follows, or else with clade richness (`richness`) or the number of undated nodes between its anchors (`topology`). A node with no counterpart in the first run adds `topology`, names the genera that moved, and is compared with the node of the first run that holds the same genera.

## Result

`stem_causes.tsv`: one row per changed stem, with both ages, the factor, the causes, whether a cause was found (`traced`) and the chain of nodes that led there.

- **175 stems** differ by more than 10%, the same genera as `results/stem_changes.tsv`; **all 175 are traced**.
- Causes, counted per stem (a stem can have several): fit 138, calibration gained 124, topology 31, identifier 20, picks 20, richness 7, calibration lost 7, no answer 4, unresolved 4, higher taxa 1. Every gained or lost calibration comes with the reason the other run lacks it.
- The most common chains: a calibration that the precedence rule dropped and the fit keeps (96), the fit alone (24), a changed identifier that made a node datable (10), and new topology with other picks (8). Where a genus moved, a smaller clade can lose its calibration because TimeTree has no age for the new picks (4, such as *Bathysolea* and *Hypleurochilus*) or resolves them only to an ancestor (4, such as *Caryophyllia*).
- Examples: *Deosergestes* (180 to 8 Ma) through topology, identifier and calibration gained, as described in #56; *Alloteuthis* (158 to 76 Ma) through richness alone, its clade of 8 genera becoming 3; *Balanus* (94 to 57 Ma) through topology and the picks of its query.

Metric 6b of the gate (version 2, [#74](https://github.com/naturalis/timegraft/issues/74)) can count a stem as explained where `traced` is `yes`, instead of where a fitted age lies anywhere on its path.
