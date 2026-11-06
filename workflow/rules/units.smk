"""From raw records to one representative per reference unit."""


rule filter_records:
    input:
        config["records"],
    output:
        kept=f"{OUT}/units/filtered.tsv",
        rejected=f"{OUT}/units/rejected.tsv",
    params:
        src=SRC,
        marker=config["marker"],
        options=config["markers"][config["marker"]]["filter"],
    log:
        f"{OUT}/logs/filter_records.log",
    conda:
        "../envs/timegraft.yaml"
    script:
        "../scripts/filter_records.py"


rule representatives:
    input:
        rules.filter_records.output.kept,
    output:
        reps=f"{OUT}/units/representatives.tsv",
        unassigned=f"{OUT}/units/unassigned.tsv",
    params:
        src=SRC,
        max_records=config["max_records_per_unit"],
        genus_agreement=config["genus_agreement"],
    log:
        f"{OUT}/logs/representatives.log",
    conda:
        "../envs/timegraft.yaml"
    script:
        "../scripts/representatives.py"
