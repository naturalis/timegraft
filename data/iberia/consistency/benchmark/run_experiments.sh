#!/usr/bin/env bash
# The backbone experiments of the gate (#56): clean builds, injected misplacements, perturbed ages,
# dropped genera, and the reproducibility builds. Run after the two full runs (see README.md).
#
#   run_experiments.sh BENCH REPO OTOL [JOBS]
#
# BENCH: working folder; holds branch/results/iberia (the full default run) and gets exp/.
# REPO:  checkout of consistency-layer whose pipeline is run. OTOL: folder with the OTT and synthesis files.
set -euo pipefail
BENCH=$1; REPO=$2; OTOL=$3; JOBS=${4:-2}
H=$(cd "$(dirname "$0")" && pwd)
U=$BENCH/branch/results/iberia/units
E=$BENCH/exp
CACHE=$BENCH/branch/resources/api_cache.json
SEEDS="1 2 3 4 5"
mkdir -p "$E"
export PYTHONPATH=$H:$REPO/src

bb() {  # bb NAME UNITS CACHE_SOURCE [backbone.py options]: build one backbone with its own copy of the cache
  local name=$1 units=$2 cache=$3; shift 3
  if [ -s "$E/$name/metadata.json" ]; then return 0; fi
  mkdir -p "$E/$name"
  cp "$cache" "$E/$name/cache.json"
  python "$H/backbone.py" --repo "$REPO" --units "$units" --otol "$OTOL" --cache "$E/$name/cache.json" \
      --out "$E/$name" "$@" > "$E/$name/stdout.log" 2>&1 || echo "FAILED $name" >> "$E/failed.txt"
}
export -f bb; export E H REPO OTOL U CACHE

# 0. the screening with the layer off, and the two clean backbones
[ -s "$E/legacy_units/screened.tsv" ] || python "$H/screen.py" --repo "$REPO" --units "$U" --out "$E/legacy_units" --overlay "$H/legacy.yaml"
for f in representatives.tsv filtered.tsv; do ln -sf "$U/$f" "$E/legacy_units/$f"; done
# the default screening once more, for the reproducibility check
[ -s "$E/default_units_again/screened.tsv" ] || python "$H/screen.py" --repo "$REPO" --units "$U" --out "$E/default_units_again"
bb legacy_clean "$E/legacy_units" "$CACHE" --overlay "$H/legacy.yaml" &
bb default_clean "$U" "$CACHE" &
wait
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
} | xargs -P "$JOBS" -I{} bash -c "{}"
[ ! -e "$E/failed.txt" ] || { echo "some builds failed:"; cat "$E/failed.txt"; }
