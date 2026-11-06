"""Dated genus-level backbone from Open Tree of Life and TimeTree."""


rule build_backbone:
    input:
        reps=rules.representatives.output.reps,
    output:
        tree=f"{OUT}/backbone/backbone.tre",
        genera=f"{OUT}/backbone/genera.tsv",
        nodes=f"{OUT}/backbone/node_ages.tsv",
    params:
        src=SRC,
        settings=config.get("backbone_build", {}),
    log:
        f"{OUT}/logs/build_backbone.log",
    conda:
        "../envs/timegraft.yaml"
    script:
        "../scripts/build_backbone.py"
