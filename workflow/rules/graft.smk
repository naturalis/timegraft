"""Grafting the dated genus subtrees onto the dated backbone."""


rule graft:
    input:
        backbone=config["backbone"],
        subtrees=expand(f"{OUT}/subtrees/shard_{{i}}.tsv", i=SHARDS),
        singles=rules.shard_genera.output.singles,
    output:
        tree=f"{OUT}/timegraft.tre",
        provenance=f"{OUT}/provenance.tsv",
    params:
        src=SRC,
        prune_empty=config["prune_empty"],
    log:
        f"{OUT}/logs/graft.log",
    conda:
        "../envs/timegraft.yaml"
    script:
        "../scripts/graft.py"
