"""Phylogenetic diversity per sample and between samples."""


rule diversity:
    input:
        tree=rules.graft.output.tree,
        samples=config.get("samples") or [],
    output:
        alpha=f"{OUT}/diversity/alpha.tsv",
        beta=f"{OUT}/diversity/beta.tsv",
        missing=f"{OUT}/diversity/units_not_in_tree.tsv",
    params:
        src=SRC,
    log:
        f"{OUT}/logs/diversity.log",
    conda:
        "../envs/timegraft.yaml"
    script:
        "../scripts/diversity.py"
