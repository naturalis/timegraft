"""Metric 5c and 5d of the gate for one setting: clade and stem ages of the subsampled backbones against the clean one.

Computes what evaluate.py computes for 5c and 5d, for any clean backbone and its five subsampled builds, so that
the settings of #108 can be compared one step at a time (README.md, "The gate with the node check rejecting").

    python subsample_by_step.py CLEAN SUB1 SUB2 SUB3 SUB4 SUB5

Each argument is a folder with backbone.tre. Prints the 90% quantile of the factor by which clade ages move (5c),
the share of clades that move by more than 1.5, and the 90% quantile for genus stems (5d).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from benchlib import Backbone, factor, quantile  # noqa: E402


def subsample_metrics(clean, subs):
    c = Backbone(Path(clean) / "backbone.tre")
    clades, stems = [], []
    for folder in subs:
        bb = Backbone(Path(folder) / "backbone.tre")
        clades += [factor(age, c.mrca_age(tips)) for tips, age in bb.node_ages().items()
                   if tips <= set(c.leaf) and age > 0]
        for genus in bb.leaf:
            node = bb.parent[bb.leaf[genus]]
            if bb.tips[node] <= set(c.leaf) and bb.age[node] > 0:
                stems.append(factor(bb.age[node], c.mrca_age(bb.tips[node])))
    return quantile(clades, 0.9), sum(x > 1.5 for x in clades) / len(clades), quantile(stems, 0.9)


if __name__ == "__main__":
    q, share, stem = subsample_metrics(sys.argv[1], sys.argv[2:])
    print(f"5c {q:.3f}\tshare over 1.5 {share:.3f}\t5d {stem:.3f}")
