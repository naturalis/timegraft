#!/usr/bin/env bash
# The detection metrics of the gate (#56) again, for #71: the default with concordance.lineage_alone_level
# at class, and a control of the same code with it off (the rule of the gate), on the clean run and on both
# injection variants: inject (merged NCBI ids, metrics 1b-1e) and inject_name (name only). Evaluate with
# evaluate_lineage_alone.py. Exits non-zero if any build failed (listed in exp/failed.txt); rerunning it
# retries only those.
#
#   run_lineage_alone.sh BENCH REPO OTOL [JOBS]
#
# BENCH: working folder, laid out as for run_experiments.sh (branch/results/iberia/units from a default run,
# branch/resources/api_cache.json); gets exp/. REPO: the checkout whose pipeline is run.
set -euo pipefail
BENCH=$1; REPO=$2; OTOL=$3; JOBS=${4:-2}
H=$(cd "$(dirname "$0")" && pwd)
U=$BENCH/branch/results/iberia/units; E=$BENCH/exp; CACHE=$BENCH/branch/resources/api_cache.json
mkdir -p "$E"; rm -f "$E/failed.txt"
export PYTHONPATH=$H:$REPO/src
OFF="--set backbone_build.concordance.lineage_alone_level=off"
bb() {  # bb NAME UNITS CACHE_SOURCE [backbone.py options]: one backbone, with its own copy of the cache
  local name=$1 units=$2 cache=$3; shift 3
  if [ -s "$E/$name/seconds.txt" ]; then return 0; fi
  rm -rf "$E/$name"; mkdir -p "$E/$name"; cp "$cache" "$E/$name/cache.json"
  if ! python "$H/backbone.py" --repo "$REPO" --units "$units" --otol "$OTOL" --cache "$E/$name/cache.json" \
      --out "$E/$name" "$@" > "$E/$name.log" 2>&1; then
    echo "FAILED $name (log: $E/$name.log)" >> "$E/failed.txt"; rm -rf "$E/$name"; return 1
  fi
  mv "$E/$name.log" "$E/$name/stdout.log"
}
export -f bb; export E H REPO OTOL U CACHE OFF
bb default_clean "$U" "$CACHE" & p1=$!
bb control_clean "$U" "$CACHE" $OFF & p2=$!
wait $p1 && wait $p2
for s in 1 2 3 4 5; do
  for mode in merge name; do
    mkdir -p "$E/$mode/s$s"
    [ -s "$E/$mode/s$s/manifest.tsv" ] || python "$H/inject.py" "$OTOL/ott3.7.3/taxonomy.tsv" \
        "$OTOL/opentree16.1/labelled_supertree.tre" "$E/default_clean/genera.tsv" "$s" \
        "$E/$mode/s$s/taxonomy.tsv" "$E/$mode/s$s/manifest.tsv" 12 "$mode"
  done
done
{
  for s in 1 2 3 4 5; do
    for mode in merge name; do
      echo "bb $mode/s$s/default $U $CACHE --taxonomy $E/$mode/s$s/taxonomy.tsv"
      echo "bb $mode/s$s/control $U $CACHE --taxonomy $E/$mode/s$s/taxonomy.tsv $OFF"
    done
  done
} | xargs -P "$JOBS" -I{} bash -c "{}" || true
if [ -e "$E/failed.txt" ]; then cat "$E/failed.txt" >&2; exit 1; fi
