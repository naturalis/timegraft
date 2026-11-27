"""Which part of the consistency layer gives its robustness (#80): metrics 5a to 5d for each ablation.

    python attribute.py BENCH OUTDIR

BENCH is a gate folder after run_experiments.sh and run_ablations.sh. For every arm (legacy, the default and each
ablation) the perturbed backbones are compared with the clean backbone of the same arm, as evaluate.py does for
the gate:

- noise (5a, 5b): exp/noise/s*/<arm> against the arm's clean build without re-queries;
- subsampling (5c, 5d): exp/sub/s*/<arm> against the arm's clean build.

Writes noise.tsv and subsample.tsv to OUTDIR (one row per arm, pooled over the five seeds, and one per arm and
seed) and prints both tables.
"""
import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchlib import Backbone, factor, quantile  # noqa: E402

SEEDS = (1, 2, 3, 4, 5)
ABLATIONS = ("precedence_ages", "mrca_anchoring", "identifiers_from_ott", "no_reanchoring", "global_crowns",
             "raw_timetree_time")


def arms(bench):
    e, o = bench / "exp", bench / "exp80"
    yield ("legacy", e / "legacy_clean", e / "legacy_clean",
           [e / f"noise/s{s}/legacy" for s in SEEDS], [e / f"sub/s{s}/legacy" for s in SEEDS])
    yield ("default", e / "default_norequery", e / "default_clean",
           [e / f"noise/s{s}/default" for s in SEEDS], [e / f"sub/s{s}/default" for s in SEEDS])
    for a in ABLATIONS:
        clean = o / a / "clean" if a == "raw_timetree_time" else e / "ablate" / a
        yield (a, o / a / "clean_norequery", clean,
               [o / a / f"noise/s{s}" for s in SEEDS], [o / a / f"sub/s{s}" for s in SEEDS])


def noise_moves(ref, folder):
    """Factors of node ages and of genus stems of one noise build against its reference (5a, 5b)."""
    bb = Backbone(folder / "backbone.tre")
    ages, clean = bb.node_ages(), ref.node_ages()
    nodes = [factor(ages[k], clean[k]) for k in ages if k in clean and clean[k] > 0 and ages[k] > 0]
    stems = [factor(bb.stem_age(g), ref.stem_age(g)) for g in bb.leaf if g in ref.leaf]
    return nodes, stems


def sub_moves(ref, folder):
    """Factors of clade ages and of genus stems of one subsampled build against its reference (5c, 5d)."""
    bb = Backbone(folder / "backbone.tre")
    clades, stems = [], []
    for tips, age in bb.node_ages().items():
        if tips <= set(ref.leaf) and age > 0:
            clades.append(factor(age, ref.mrca_age(tips)))
    for genus in bb.leaf:
        node = bb.parent[bb.leaf[genus]]
        if bb.tips[node] <= set(ref.leaf) and bb.age[node] > 0:
            stems.append(factor(bb.age[node], ref.mrca_age(bb.tips[node])))
    return clades, stems


def summary(values):
    return (f"{statistics.median(values):.3f}", f"{quantile(values, 0.9):.3f}",
            f"{sum(x > 1.5 for x in values) / len(values):.3f}")


def main(bench, out):
    out.mkdir(parents=True, exist_ok=True)
    noise, sub = [], []
    for arm, noise_ref, sub_ref, noise_dirs, sub_dirs in arms(bench):
        nref, sref = Backbone(noise_ref / "backbone.tre"), Backbone(sub_ref / "backbone.tre")
        pooled = {"n": [], "ns": [], "c": [], "cs": []}
        for seed, nd, sd in zip(SEEDS, noise_dirs, sub_dirs):
            n, ns = noise_moves(nref, nd)
            c, cs = sub_moves(sref, sd)
            for k, v in zip(("n", "ns", "c", "cs"), (n, ns, c, cs)):
                pooled[k] += v
            noise.append({"arm": arm, "seed": seed, "nodes": len(n), "median": summary(n)[0],
                          "p90 (5a)": summary(n)[1], "share_over_1.5 (5b)": summary(n)[2],
                          "stem_p90": summary(ns)[1]})
            sub.append({"arm": arm, "seed": seed, "clades": len(c), "median": summary(c)[0],
                        "p90 (5c)": summary(c)[1], "share_over_1.5": summary(c)[2], "stem_p90 (5d)": summary(cs)[1]})
        n, ns, c, cs = (pooled[k] for k in ("n", "ns", "c", "cs"))
        noise.append({"arm": arm, "seed": "all", "nodes": len(n), "median": summary(n)[0], "p90 (5a)": summary(n)[1],
                      "share_over_1.5 (5b)": summary(n)[2], "stem_p90": summary(ns)[1]})
        sub.append({"arm": arm, "seed": "all", "clades": len(c), "median": summary(c)[0], "p90 (5c)": summary(c)[1],
                    "share_over_1.5": summary(c)[2], "stem_p90 (5d)": summary(cs)[1]})
    for name, rows in (("noise.tsv", noise), ("subsample.tsv", sub)):
        with open(out / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(f"== {name}")
        for r in rows:
            if r["seed"] == "all":
                print("\t".join(str(v) for v in r.values()))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
