# Re-anchoring within the target family (#70): the injection part of the gate

[#70](https://github.com/naturalis/timegraft/issues/70) attaches a re-anchored genus within its target taxon, at the majority clade of its own species that the synthesis places inside it, instead of at the target's node. On the injected errors of the gate, that is where the remaining stem misses (1f) and most of the other stems moved (1g) came from (`../REPORT.md`, "Two things the correction does not yet do", point 2). This folder has the injection metrics of gate version 2 (1a to 1m) after #70, next to those with #69 alone (the rerun of #73, `../results_73/`).

## Method

Only the builds that #70 can change were redone: the default arm of the three injections (`inject`, `inject_own`, `inject_name`), five replicates each, with `issue-70-within-target`. Everything else is the benchmark of #73:

- the legacy arms of the injections: re-anchoring is off in `legacy.yaml`, so #70 does not touch them;
- `default_clean`: rebuilt with #70 and byte-identical to `main` in `backbone.tre`, `node_ages.tsv` and `crown_ages.tsv` (only the new columns of `genera.tsv` and `concordance.tsv` differ), because the one genus moved in the Iberian run, *Onchidella*, has none of its species inside Onchidiidae in the synthesis and stays at the family node;
- `legacy_clean`: the rebuild of #84, identical to that of #73 apart from the columns of #84.

Each rebuilt arm started from the cache of its build in #73, with the same injected taxonomy:

```
PYTHONPATH=$B:src python $B/backbone.py --repo . --units $BENCH73/branch/results/iberia/units \
    --otol $BENCH73/branch/resources/otol --cache <copy of exp/K/sN/default/cache.json of #73> \
    --out exp/K/sN/default --taxonomy $BENCH73/exp/K/sN/taxonomy.tsv
PYTHONPATH=src python $B/evaluate.py $BENCH70 out           # gate version 2
```

with `B=data/iberia/consistency/benchmark`. Only metrics 1a to 1m are reported here: the other metrics of the evaluation were not rerun (they come from #73), and 7c and 7d compare clean backbones of #84 and #70 with runs of #73 that predate their columns. The full gate, with the downstream metrics of #81, follows once #81 is done.

## Result

`injection_metrics.tsv`:

| Metric | With #69 alone | With #70 |
| --- | --- | --- |
| 1d flagged genera back in their own class | 119 of 120 | 120 of 120 |
| **1f corrected genera with a stem within 1.5** | 97 of 119 (0.815) | **119 of 120 (0.992)** |
| **1g other stems moved by more than 10%** (legacy 248; at most 124) | 142, fail | **4, pass** |
| 1j `inject_own`: stem within 1.5 | 97 of 119 | 119 of 120 |
| 1m `inject_name`: stem within 1.5 (report) | 86 of 120 | 86 of 120 |

Every other injection metric is unchanged and passes. The other stems moved are 98 in `inject_name` against 184 for legacy, as before.

`within_target.tsv` says where each moved genus went. Of the 120 injected genera in `inject` and `inject_own`, 95 were attached at the majority clade of their species inside the target, 15 at the MRCA of two such species, 4 at the parent of one, 3 at the target's node because the synthesis has none of their species inside it, and 3 at the nearest enclosing node because the target itself is not in the synthesis (neither as a node nor through members), where the genus is not placed within that larger clade. *Onchidella* (the five "not injected" rows) stays at the family node in every replicate.

## What this does not show

- **`inject_name` is unchanged**: none of its 120 genera is placed within its target (117 stay at the target's node, 3 at the enclosing node). In that variant only the name moves: the genus's species stay under its own taxon (renamed), so the taxon the pipeline matched holds only the host's species, and none of them is inside the target. Finding the genus's species there would need a search by species name, which the issue mentions as an option and this change does not make.
- **The merge injections are kinder than real merges.** `inject.py` changes the taxonomy only, so the synthesis still places the genus's species in their own family, where #70 finds them. Where OTT itself merged two genera, the synthesis follows OTT: *Onchidella*'s species are not inside Onchidiidae, and #70 does nothing for it. The gain above is the gain where the synthesis knows better than the taxonomy.
- **Fewer than three species inside the target.** The issue proposed keeping the target's node when fewer than `anchoring.min_species` species are inside it. This change uses their MRCA (or the parent of a single one), as ordinary anchoring does since #51; that concerns 19 of the 120 injected genera.
