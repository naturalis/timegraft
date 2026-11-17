#!/usr/bin/env bash
# The downstream experiments of gate version 2 (#81): the genus stage, grafting and placement run with
# downstream.py on the clean backbones and on replicates 1 and 2 of the injection, the age perturbation and the
# subsampling, for both settings. Run after run_experiments.sh, on the same BENCH. Writes exp/downstream/ as
# GATE.md lays it out; evaluate.py reads it. Runs one at a time, each with JOBS shard processes. A finished run
# is not repeated, so rerunning the script resumes; exits non-zero if any run failed (exp/downstream/failed.txt).
#
#   run_downstream.sh BENCH REPO [JOBS]
#
# Needs MAFFT and VSEARCH on the PATH. The clean runs go first: they fill the cache of build_genus() results
# (exp/downstream/memo) that the experiments draw on wherever a genus's inputs are unchanged.
set -euo pipefail
BENCH=$1; REPO=$2; JOBS=${3:-2}
H=$(cd "$(dirname "$0")" && pwd)
E=$BENCH/exp; D=$E/downstream; U=$BENCH/branch/results/iberia/units
mkdir -p "$D"; rm -f "$D/failed.txt"
export PYTHONPATH=$H

ds() {  # ds NAME BACKBONE ARM [downstream.py options]: one downstream run
  local name=$1 backbone=$2 arm=$3; shift 3
  local units=$U overlay=()
  if [ "$arm" = legacy ]; then units=$E/legacy_units; overlay=(--overlay "$H/legacy.yaml"); fi
  [ -s "$D/$name/seconds.txt" ] && return 0
  mkdir -p "$(dirname "$D/$name")"
  echo "$(date '+%H:%M') $name"
  if ! python "$H/downstream.py" --repo "$REPO" --backbone "$backbone" --units "$units" --out "$D/$name" \
      --memo "$D/memo" --jobs "$JOBS" --threads "$JOBS" "${overlay[@]}" "$@" > "$D/${name//\//_}.log" 2>&1; then
    echo "FAILED $name (log: $D/${name//\//_}.log)" >> "$D/failed.txt"
  fi
}

for arm in default legacy; do ds clean/$arm "$E/${arm}_clean" $arm; done
for s in 1 2; do
  for arm in default legacy; do
    ds inject/s$s/$arm "$E/inject/s$s/$arm" $arm
    ds sub/s$s/$arm "$E/sub/s$s/$arm" $arm --genera "$E/sub/s$s/units/genera.txt"
    ds noise/s$s/$arm "$E/noise/s$s/$arm" $arm
  done
done
if [ -e "$D/failed.txt" ]; then cat "$D/failed.txt" >&2; exit 1; fi
