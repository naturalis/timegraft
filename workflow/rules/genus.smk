"""Per-genus subtrees, processed in shards so that genera run in parallel."""


rule shard_genera:
    input:
        reps=rules.screen_units.output.reps,
        records=rules.filter_records.output.kept,
        backbone=BACKBONE,
    output:
        shards=expand(f"{OUT}/genera/shard_{{i}}.tsv", i=SHARDS),
        singles=f"{OUT}/genera/singles.tsv",
        missing=f"{OUT}/genera/missing_from_backbone.tsv",
        taxonomy=f"{OUT}/genera/taxonomy.tsv",
    params:
        src=SRC,
        n_shards=len(SHARDS),
        rooting=ROOTING,
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
        timetree_crowns=rules.build_backbone.output.crowns if BUILD_BACKBONE else [],
        taxonomy=rules.shard_genera.output.taxonomy,
        rate_table=MARKER_CFG["crown_age"].get("rate_table") or [],
        lsd2=LSD2_BIN if USE_LSD2 else [],
    output:
        f"{OUT}/subtrees/shard_{{i}}.tsv",
    params:
        src=SRC,
        marker_cfg=MARKER_CFG,
        rooting=ROOTING,
    log:
        f"{OUT}/logs/build_shard_{{i}}.log",
    conda:
        "../envs/chronos.yaml" if USE_CHRONOS else "../envs/timegraft.yaml"
    script:
        "../scripts/build_shard.py"
