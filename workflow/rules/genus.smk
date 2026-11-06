"""Per-genus subtrees, processed in shards so that genera run in parallel."""


rule shard_genera:
    input:
        reps=rules.representatives.output.reps,
        backbone=config["backbone"],
    output:
        shards=expand(f"{OUT}/genera/shard_{{i}}.tsv", i=SHARDS),
        singles=f"{OUT}/genera/singles.tsv",
        missing=f"{OUT}/genera/missing_from_backbone.tsv",
    params:
        src=SRC,
        n_shards=len(SHARDS),
        n_outgroup=config["n_outgroup"],
    log:
        f"{OUT}/logs/shard_genera.log",
    conda:
        "../envs/timegraft.yaml"
    script:
        "../scripts/shard_genera.py"


rule build_shard:
    input:
        shard=f"{OUT}/genera/shard_{{i}}.tsv",
        crown_ages=config["crown_ages"] if config.get("crown_ages") else [],
    output:
        f"{OUT}/subtrees/shard_{{i}}.tsv",
    params:
        src=SRC,
        marker_cfg=config["markers"][config["marker"]],
    log:
        f"{OUT}/logs/build_shard_{{i}}.log",
    conda:
        "../envs/timegraft.yaml"
    script:
        "../scripts/build_shard.py"
