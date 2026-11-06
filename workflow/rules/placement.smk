"""Matching metabarcoding sequences to units, and placing the unmatched ones."""


rule place_asvs:
    input:
        tree=rules.place_units.output.tree,
        asvs=config.get("asvs") or [],
        asv_table=config.get("asv_table") or [],
        records=rules.filter_records.output.kept,
        reps=rules.representatives.output.reps,
    output:
        tree=f"{OUT}/placement/timegraft_placed.tre",
        placements=f"{OUT}/placement/placements.tsv",
        samples=f"{OUT}/placement/sample_tips.tsv",
    params:
        src=SRC,
        settings=config["markers"][config["marker"]]["placement"],
        rate=config["markers"][config["marker"]]["crown_age"].get("rate"),
    threads: 4
    log:
        f"{OUT}/logs/place_asvs.log",
    conda:
        "../envs/timegraft.yaml"
    script:
        "../scripts/place_asvs.py"
