#!/usr/bin/env bash
# The backbone experiments of the gate (#56): clean builds, injected misplacements, perturbed ages,
# dropped genera, and the reproducibility builds. Run after the two full runs (see README.md).
# Since version 2 of the gate (#74) this includes the two injection variants, which version 1 ran
# separately (run_posthoc.sh); builds that already exist are not repeated.
# Exits non-zero if any build failed (listed in exp/failed.txt); rerunning it retries only those.
#
#   run_experiments.sh BENCH REPO OTOL [JOBS] [offline]
#
# BENCH: working folder; holds branch/results/iberia (the full default run) and gets exp/.
# REPO:  checkout of consistency-layer whose pipeline is run. OTOL: folder with the OTT and synthesis files.
# Every build starts from the benchmark cache (timetree_cache.json.gz: the TimeTree and ChecklistBank answers
# of all earlier builds, #87) together with the cache of the full run. With "offline", no build may ask
# TimeTree or ChecklistBank anything: a question the caches cannot answer makes the build fail.
set -euo pipefail
BENCH=$1; REPO=$2; OTOL=$3; JOBS=${4:-2}; OFFLINE=${5:-}
H=$(cd "$(dirname "$0")" && pwd)
U=$BENCH/branch/results/iberia/units
E=$BENCH/exp
SEEDS="1 2 3 4 5"
mkdir -p "$E"
rm -f "$E/failed.txt"  # failures of an earlier attempt are retried below
export PYTHONPATH=$H:$REPO/src
CACHE=$E/start_cache.json
# made again on every invocation, so that answers added to the benchmark cache since the last one are used
python "$H/merge_caches.py" "$CACHE" "$H/timetree_cache.json.gz" "$BENCH/branch/resources/api_cache.json"
case "$OFFLINE" in ""|offline) ;; *) echo "the fifth argument is 'offline' or nothing, not '$OFFLINE'" >&2; exit 2;; esac

bb() {  # bb NAME UNITS CACHE_SOURCE [backbone.py options]: build one backbone with its own copy of the cache
  # A build counts as done only when backbone.py has written seconds.txt, its last step. A failed build
  # is recorded in failed.txt, its partial outputs are removed, and the function returns non-zero.
  local name=$1 units=$2 cache=$3; shift 3
  if [ -s "$E/$name/seconds.txt" ]; then return 0; fi
  rm -rf "$E/$name"; mkdir -p "$E/$name"
  cp "$cache" "$E/$name/cache.json"
  if ! python "$H/backbone.py" --repo "$REPO" --units "$units" --otol "$OTOL" --cache "$E/$name/cache.json" \
      --out "$E/$name" ${OFFLINE:+--offline} "$@" > "$E/$name.log" 2>&1; then
    echo "FAILED $name (log: $E/$name.log)" >> "$E/failed.txt"
    rm -rf "$E/$name"
    return 1
  fi
  mv "$E/$name.log" "$E/$name/stdout.log"
}
export -f bb; export E H REPO OTOL U CACHE OFFLINE

# 0. the screening with the layer off, and the two clean backbones
[ -s "$E/legacy_units/screened.tsv" ] || python "$H/screen.py" --repo "$REPO" --units "$U" --out "$E/legacy_units" --overlay "$H/legacy.yaml"
for f in representatives.tsv filtered.tsv; do ln -sf "$U/$f" "$E/legacy_units/$f"; done
# the default screening once more, for the reproducibility check
[ -s "$E/default_units_again/screened.tsv" ] || python "$H/screen.py" --repo "$REPO" --units "$U" --out "$E/default_units_again"
bb legacy_clean "$E/legacy_units" "$CACHE" --overlay "$H/legacy.yaml" & p1=$!
bb default_clean "$U" "$CACHE" & p2=$!
wait $p1 && wait $p2  # every later build depends on these two: stop here if either failed
# ages of the default without re-queries, which change no age: the reference of the offline noise runs
bb default_norequery "$U" "$E/default_clean/cache.json" --set backbone_build.ages.requery.mode=off
# every query either arm asks, for the offline runs
python - "$E/legacy_clean/cache.json" "$E/default_clean/cache.json" "$E/union_cache.json" <<'PY'
import json, sys
a, b = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
json.dump({**a, **b}, open(sys.argv[3], "w"))
PY

# 1. inputs of the experiments
for s in $SEEDS; do
  mkdir -p "$E/inject/s$s" "$E/noise/s$s" "$E/sub/s$s"
  [ -s "$E/inject/s$s/manifest.tsv" ] || python "$H/inject.py" "$OTOL/ott3.7.3/taxonomy.tsv" \
      "$OTOL/opentree16.1/labelled_supertree.tre" "$E/default_clean/genera.tsv" "$s" \
      "$E/inject/s$s/taxonomy.tsv" "$E/inject/s$s/manifest.tsv"
  [ -s "$E/noise/s$s/cache.json" ] || python "$H/perturb_cache.py" "$E/union_cache.json" "$s" "$E/noise/s$s/cache.json"
  [ -s "$E/sub/s$s/units/genera.txt" ] || python "$H/subsample.py" "$U" "$s" "$E/sub/s$s/units"
  for mode in name own; do  # the two injection variants (metrics 1h to 1m of version 2)
    mkdir -p "$E/inject_$mode/s$s"
    [ -s "$E/inject_$mode/s$s/manifest.tsv" ] || python "$H/inject.py" "$OTOL/ott3.7.3/taxonomy.tsv" \
        "$OTOL/opentree16.1/labelled_supertree.tre" "$E/default_clean/genera.tsv" "$s" \
        "$E/inject_$mode/s$s/taxonomy.tsv" "$E/inject_$mode/s$s/manifest.tsv" 12 "$mode"
  done
done

# 2. the builds, JOBS at a time
{
  for s in $SEEDS; do
    echo "bb inject/s$s/legacy $E/legacy_units $CACHE --overlay $H/legacy.yaml --taxonomy $E/inject/s$s/taxonomy.tsv"
    echo "bb inject/s$s/default $U $CACHE --taxonomy $E/inject/s$s/taxonomy.tsv"
    echo "bb noise/s$s/legacy $E/legacy_units $E/noise/s$s/cache.json --overlay $H/legacy.yaml --offline"
    echo "bb noise/s$s/default $U $E/noise/s$s/cache.json --set backbone_build.ages.requery.mode=off --offline"
    echo "bb sub/s$s/legacy $E/sub/s$s/units $CACHE --overlay $H/legacy.yaml"
    echo "bb sub/s$s/default $E/sub/s$s/units $CACHE"
    for mode in name own; do
      echo "bb inject_$mode/s$s/default $U $CACHE --taxonomy $E/inject_$mode/s$s/taxonomy.tsv"
      echo "bb inject_$mode/s$s/legacy $E/legacy_units $CACHE --overlay $H/legacy.yaml --taxonomy $E/inject_$mode/s$s/taxonomy.tsv"
    done
  done
  # one feature of the layer switched back at a time (report only)
  echo "bb ablate/identifiers_from_ott $U $CACHE --set backbone_build.col.join=ott"
  echo "bb ablate/no_reanchoring $U $CACHE --set backbone_build.reanchor.mode=report"
  echo "bb ablate/mrca_anchoring $U $CACHE --set backbone_build.anchoring.rule=mrca --set backbone_build.anchoring.prefer_sampled=false"
  echo "bb ablate/precedence_ages $U $CACHE --set backbone_build.ages.method=precedence"
  echo "bb ablate/global_crowns $U $CACHE --set backbone_build.crown_query.basis=global"
  # reproducibility: the default again from the cache it left, with and without the network
  echo "bb repro/again $U $E/default_clean/cache.json"
  echo "bb repro/offline $U $E/default_clean/cache.json --offline"
} | tr '\n' '\0' | xargs -0 -n 1 -P "$JOBS" bash -c 'eval "$1"' _ || true  # xargs runs every job; failures are collected below
if [ -e "$E/failed.txt" ]; then
  echo "some builds failed; rerun this script to retry them:" >&2
  cat "$E/failed.txt" >&2
  exit 1
fi
