#!/usr/bin/env bash
# Experiments added after the first results of the gate, to explain them (#56). Not part of the gate.
#   run_posthoc.sh BENCH REPO OTOL [JOBS]
set -euo pipefail
BENCH=$1; REPO=$2; OTOL=$3; JOBS=${4:-2}
H=$(cd "$(dirname "$0")" && pwd)
U=$BENCH/branch/results/iberia/units; E=$BENCH/exp; CACHE=$BENCH/branch/resources/api_cache.json
export PYTHONPATH=$H:$REPO/src
bb() {
  local name=$1 units=$2 cache=$3; shift 3
  if [ -s "$E/$name/metadata.json" ]; then return 0; fi
  mkdir -p "$E/$name"; cp "$cache" "$E/$name/cache.json"
  python "$H/backbone.py" --repo "$REPO" --units "$units" --otol "$OTOL" --cache "$E/$name/cache.json" \
      --out "$E/$name" "$@" > "$E/$name/stdout.log" 2>&1 || echo "FAILED $name" >> "$E/failed.txt"
}
export -f bb; export E H REPO OTOL U CACHE
for s in 1 2 3 4 5; do
  for mode in name own; do
    mkdir -p "$E/inject_$mode/s$s"
    [ -s "$E/inject_$mode/s$s/manifest.tsv" ] || python "$H/inject.py" "$OTOL/ott3.7.3/taxonomy.tsv" \
        "$OTOL/opentree16.1/labelled_supertree.tre" "$E/default_clean/genera.tsv" "$s" \
        "$E/inject_$mode/s$s/taxonomy.tsv" "$E/inject_$mode/s$s/manifest.tsv" 12 "$mode"
  done
done
{
  for s in 1 2 3 4 5; do
    for mode in name own; do
      echo "bb inject_$mode/s$s/default $U $CACHE --taxonomy $E/inject_$mode/s$s/taxonomy.tsv"
      echo "bb inject_$mode/s$s/legacy $E/legacy_units $CACHE --overlay $H/legacy.yaml --taxonomy $E/inject_$mode/s$s/taxonomy.tsv"
    done
  done
} | xargs -P "$JOBS" -I{} bash -c "{}"
