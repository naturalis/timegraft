#!/usr/bin/env bash
# The injection part of the gate (#56) again, for #69: the builds of run_experiments.sh that metrics 1b-1g
# need, and a control arm of the same code with reanchor.check_ids off (with its own clean build), so that the
# effect of the identifier check is measured on one machine with one cache. Evaluate with evaluate_injection.py.
# Exits non-zero if any build failed (listed in exp/failed.txt); rerunning it retries only those.
#
#   run_injection.sh BENCH REPO OTOL [JOBS]
#
# BENCH: working folder, laid out as for run_experiments.sh (branch/results/iberia/units from a default run,
# branch/resources/api_cache.json); gets exp/. REPO: the checkout whose pipeline is run.
set -euo pipefail
BENCH=$1; REPO=$2; OTOL=$3; JOBS=${4:-2}
H=$(cd "$(dirname "$0")" && pwd)
U=$BENCH/branch/results/iberia/units; E=$BENCH/exp; CACHE=$BENCH/branch/resources/api_cache.json
mkdir -p "$E"; rm -f "$E/failed.txt"
export PYTHONPATH=$H:$REPO/src
bb() {
  local name=$1 units=$2 cache=$3; shift 3
  if [ -s "$E/$name/seconds.txt" ]; then return 0; fi
  rm -rf "$E/$name"; mkdir -p "$E/$name"; cp "$cache" "$E/$name/cache.json"
  if ! python "$H/backbone.py" --repo "$REPO" --units "$units" --otol "$OTOL" --cache "$E/$name/cache.json" \
      --out "$E/$name" "$@" > "$E/$name.log" 2>&1; then
    echo "FAILED $name (log: $E/$name.log)" >> "$E/failed.txt"; rm -rf "$E/$name"; return 1
  fi
  mv "$E/$name.log" "$E/$name/stdout.log"
}
export -f bb; export E H REPO OTOL U CACHE
[ -s "$E/legacy_units/screened.tsv" ] || python "$H/screen.py" --repo "$REPO" --units "$U" --out "$E/legacy_units" --overlay "$H/legacy.yaml"
for f in representatives.tsv filtered.tsv; do ln -sf "$U/$f" "$E/legacy_units/$f"; done
bb legacy_clean "$E/legacy_units" "$CACHE" --overlay "$H/legacy.yaml" & p1=$!
bb default_clean "$U" "$CACHE" & p2=$!
wait $p1 && wait $p2
bb default_clean_control "$U" "$CACHE" --set backbone_build.reanchor.check_ids=false
for s in 1 2 3 4 5; do
  mkdir -p "$E/inject/s$s"
  [ -s "$E/inject/s$s/manifest.tsv" ] || python "$H/inject.py" "$OTOL/ott3.7.3/taxonomy.tsv" \
      "$OTOL/opentree16.1/labelled_supertree.tre" "$E/default_clean/genera.tsv" "$s" \
      "$E/inject/s$s/taxonomy.tsv" "$E/inject/s$s/manifest.tsv"
done
{
  for s in 1 2 3 4 5; do
    echo "bb inject/s$s/legacy $E/legacy_units $CACHE --overlay $H/legacy.yaml --taxonomy $E/inject/s$s/taxonomy.tsv"
    echo "bb inject/s$s/default $U $CACHE --taxonomy $E/inject/s$s/taxonomy.tsv"
    echo "bb inject_control/s$s/default $U $CACHE --taxonomy $E/inject/s$s/taxonomy.tsv --set backbone_build.reanchor.check_ids=false"
  done
} | xargs -P "$JOBS" -I{} bash -c "{}" || true
[ -e "$E/failed.txt" ] && { cat "$E/failed.txt" >&2; exit 1; } || true
