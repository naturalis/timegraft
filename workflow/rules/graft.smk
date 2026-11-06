"""Grafting the dated genus subtrees onto the dated backbone."""


rule graft:
    input:
        backbone=BACKBONE,
        subtrees=expand(f"{OUT}/subtrees/shard_{{i}}.tsv", i=SHARDS),
        singles=rules.shard_genera.output.singles,
    output:
        tree=f"{OUT}/grafted.tre",
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


rule place_units:
    input:
        tree=rules.graft.output.tree,
        unassigned=rules.representatives.output.unassigned,
        records=rules.filter_records.output.kept,
        reps=rules.representatives.output.reps,
    output:
        tree=f"{OUT}/timegraft.tre",
        placements=f"{OUT}/units_placed.tsv",
    params:
        src=SRC,
        enabled=config.get("place_unassigned_units", True),
        settings=config["markers"][config["marker"]]["placement"],
        rate=config["markers"][config["marker"]]["crown_age"].get("rate"),
    threads: 4
    log:
        f"{OUT}/logs/place_units.log",
    conda:
        "../envs/timegraft.yaml"
    script:
        "../scripts/place_units.py"
