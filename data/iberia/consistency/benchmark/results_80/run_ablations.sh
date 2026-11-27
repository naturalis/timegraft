#!/bin/bash
# Noise and subsampling replicates for each ablation of the consistency layer (#80).
#
#   run_ablations.sh BENCH REPO OTOL [JOBS]
#
# BENCH is a gate folder after run_experiments.sh (README.md of the benchmark): it needs exp/ablate/*,
# exp/union_cache.json, exp/sub/s*/units, the units of the full run in REPO/results/iberia/units and its cache
# in REPO/resources/api_cache.json. The builds start from that cache merged with the benchmark cache and the
# caches of exp/ablate/*. Builds go to BENCH/exp80/<ablation>/{clean_norequery,noise/s*,sub/s*}.
# A build with seconds.txt is skipped, so the script can be run again after an interruption or a failure;
# it exits non-zero if any build of this attempt failed (listed in exp80/failed.txt).
set -u
BENCH=$1; REPO=$2; OTOL=$3; JOBS=${4:-2}
H=$(cd "$(dirname "$0")/.." && pwd); E=$BENCH/exp; O=$BENCH/exp80; U=$REPO/results/iberia/units
export PYTHONPATH=$H:$REPO/src
mkdir -p "$O"
rm -f "$O/failed.txt"  # failures of an earlier attempt are retried below, as in run_experiments.sh
START=$O/start_cache.json
# with the caches of the gate's clean ablation builds, which hold the questions only those settings ask
[ -s "$START" ] || python "$H/merge_caches.py" "$START" "$H/timetree_cache.json.gz" "$REPO/resources/api_cache.json" "$E"/ablate/*/cache.json
declare -A SET=(
  [precedence_ages]="--set backbone_build.ages.method=precedence"
  [mrca_anchoring]="--set backbone_build.anchoring.rule=mrca --set backbone_build.anchoring.prefer_sampled=false"
  [identifiers_from_ott]="--set backbone_build.col.join=ott"
  [no_reanchoring]="--set backbone_build.reanchor.mode=report"
  [global_crowns]="--set backbone_build.crown_query.basis=global"
  [raw_timetree_time]="--set backbone_build.ages.timetree_time.use=raw"
)
bb() {  # name units cache [args...]
  local name=$1 units=$2 cache=$3; shift 3
  [ -s "$O/$name/seconds.txt" ] && return 0
  rm -rf "$O/$name"; mkdir -p "$O/$name"; cp "$cache" "$O/$name/cache.json"
  local t0=$(date +%s)
  if python "$H/backbone.py" --repo "$REPO" --units "$units" --otol "$OTOL" --cache "$O/$name/cache.json" \
      --out "$O/$name" "$@" > "$O/$name.log" 2>&1; then
    mv "$O/$name.log" "$O/$name/stdout.log"; echo $(( $(date +%s) - t0 )) > "$O/$name/seconds.txt"
  else echo "$name" >> "$O/failed.txt"; fi
}
export -f bb; export O H REPO OTOL

# 1. the clean build of the raw_timetree_time ablation, which the gate does not make
bb raw_timetree_time/clean "$U" "$START" ${SET[raw_timetree_time]}

# 2. the references of the noise runs: each ablation without re-queries (as default_norequery)
for a in "${!SET[@]}"; do
  bb "$a/clean_norequery" "$U" "$START" ${SET[$a]} --set backbone_build.ages.requery.mode=off
done

# 3. noise caches: every query of the earlier union and of the ablations, perturbed per seed. The factor of a
#    query depends only on the seed and its key, so queries in the gate's noise caches keep their factors.
python - "$O/union_cache.json" "$E/union_cache.json" "$E"/ablate/*/cache.json "$O"/*/clean_norequery/cache.json \
    "$O/raw_timetree_time/clean/cache.json" <<'PY'
import json, sys
out, *ins = sys.argv[1:]
u = {}
for p in ins:
    for k, v in json.load(open(p)).items():
        u.setdefault(k, v)
json.dump(u, open(out, "w")); print(len(u), "answers in the union")
PY
for s in 1 2 3 4 5; do
  [ -s "$O/noise_s$s.json" ] || python "$H/perturb_cache.py" "$O/union_cache.json" "$s" "$O/noise_s$s.json"
done

# 4. replicates: noise offline from the perturbed caches, subsamples from the start cache (may ask TimeTree)
for a in "${!SET[@]}"; do for s in 1 2 3 4 5; do
  echo "bb $a/noise/s$s $U $O/noise_s$s.json ${SET[$a]} --set backbone_build.ages.requery.mode=off --offline"
done; done | xargs -P "$JOBS" -I{} bash -c '{}'
for a in "${!SET[@]}"; do for s in 1 2 3 4 5; do
  echo "bb $a/sub/s$s $E/sub/s$s/units $START ${SET[$a]}"
done; done | xargs -P "$JOBS" -I{} bash -c '{}'
[ -s "$O/failed.txt" ] && { echo "failed:"; cat "$O/failed.txt"; exit 1; }
echo "all builds done"
