"""Build the backbone outside Snakemake, with the unchanged workflow script (#56).

    python backbone.py --repo REPO --units DIR --otol DIR --cache FILE --out DIR
                       [--overlay FILE ...] [--taxonomy FILE] [--set key.path=value ...] [--offline]

Settings are ``backbone_build`` of ``config/config.yaml`` and ``config/iberia.yaml`` of the repository,
with the overlays on top (e.g. ``legacy.yaml``) and then the ``--set`` values. ``--units`` holds
``representatives.tsv``, ``filtered.tsv`` and ``screened.tsv``. ``--taxonomy`` replaces the OTT taxonomy
(for injected errors). ``--offline`` points the proxies at a closed port, so that any request fails.
"""

import argparse
import os
import runpy
import sys
import time
import types
from pathlib import Path

import yaml


def merge(a, b):
    for k, v in b.items():
        a[k] = merge(a.get(k) or {}, v) if isinstance(v, dict) and isinstance(a.get(k) or {}, dict) else v
    return a


class NS(types.SimpleNamespace):
    def __getitem__(self, i):
        return list(vars(self).values())[i]


def settings_for(repo, overlays=(), sets=()):
    cfg = {}
    for path in [f"{repo}/config/config.yaml", f"{repo}/config/iberia.yaml", *overlays]:
        merge(cfg, yaml.safe_load(open(path)) or {})
    for item in sets:
        key, value = item.split("=", 1)
        node = cfg
        *parents, leaf = key.split(".")
        for p in parents:
            node = node.setdefault(p, {})
        node[leaf] = yaml.safe_load(value)
    return cfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--units", required=True)
    ap.add_argument("--otol", required=True)
    ap.add_argument("--cache", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--overlay", action="append", default=[])
    ap.add_argument("--taxonomy")
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()
    if args.offline:
        for name in ("HTTPS_PROXY", "HTTP_PROXY", "https_proxy", "http_proxy", "ALL_PROXY"):
            os.environ[name] = "http://127.0.0.1:9"
        os.environ.pop("NO_PROXY", None)
        os.environ.pop("no_proxy", None)
    cfg = settings_for(args.repo, args.overlay, args.set)
    settings = dict(cfg["backbone_build"])
    settings["cache"] = args.cache
    out, units = Path(args.out), Path(args.units)
    out.mkdir(parents=True, exist_ok=True)
    synthesis, taxonomy = settings.get("synthesis", "opentree16.1"), settings.get("taxonomy", "ott3.7.3")
    overrides = (settings.get("reanchor") or {}).get("overrides") or []
    annotations = f"{args.otol}/{synthesis}/annotations.json" if settings.get("synthesis_support", True) else []
    sm = types.SimpleNamespace(
        params=NS(src=f"{args.repo}/src", settings=settings),
        input=NS(reps=str(units / "representatives.tsv"), records=str(units / "filtered.tsv"),
                 screened=str(units / "screened.tsv"), synth=f"{args.otol}/{synthesis}/labelled_supertree.tre",
                 taxonomy=args.taxonomy or f"{args.otol}/{taxonomy}/taxonomy.tsv", annotations=annotations,
                 overrides=overrides),
        output=NS(tree=str(out / "backbone.tre"), genera=str(out / "genera.tsv"), nodes=str(out / "node_ages.tsv"),
                  crowns=str(out / "crown_ages.tsv"), metadata=str(out / "metadata.json"),
                  col_species=str(out / "col_species.tsv"), concordance=str(out / "concordance.tsv")),
        log=[str(out / "build_backbone.log")],
    )
    start = time.time()
    runpy.run_path(f"{args.repo}/workflow/scripts/build_backbone.py", init_globals={"snakemake": sm},
                   run_name="__main__")
    (out / "seconds.txt").write_text(f"{time.time() - start:.1f}\n")


if __name__ == "__main__":
    sys.exit(main())
