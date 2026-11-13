"""Run the screening outside Snakemake, with the unchanged workflow script (#56).

    python screen.py --repo REPO --units DIR --out DIR [--overlay FILE ...] [--threads N]

``--units`` holds ``representatives.tsv``, ``unassigned.tsv`` and ``filtered.tsv``; writes ``screened.tsv``,
``screened_representatives.tsv``, ``screened_unassigned.tsv`` and ``screening_calibration.tsv``.
"""

import argparse
import runpy
import time
import types
from pathlib import Path

from backbone import NS, settings_for

ap = argparse.ArgumentParser()
ap.add_argument("--repo", required=True)
ap.add_argument("--units", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--overlay", action="append", default=[])
ap.add_argument("--threads", type=int, default=2)
args = ap.parse_args()
cfg = settings_for(args.repo, args.overlay)
marker = cfg["markers"][cfg["marker"]]
out, units = Path(args.out), Path(args.units)
out.mkdir(parents=True, exist_ok=True)
sm = types.SimpleNamespace(
    params=NS(src=f"{args.repo}/src", marker=cfg["marker"], settings=marker.get("screening", {}),
              filter=marker["filter"]),
    input=NS(reps=str(units / "representatives.tsv"), unassigned=str(units / "unassigned.tsv"),
             records=str(units / "filtered.tsv")),
    output=NS(screened=str(out / "screened.tsv"), reps=str(out / "screened_representatives.tsv"),
              unassigned=str(out / "screened_unassigned.tsv"), calibration=str(out / "screening_calibration.tsv")),
    threads=args.threads, log=[str(out / "screen_units.log")],
)
start = time.time()
runpy.run_path(f"{args.repo}/workflow/scripts/screen_units.py", init_globals={"snakemake": sm}, run_name="__main__")
(out / "seconds.txt").write_text(f"{time.time() - start:.1f}\n")
