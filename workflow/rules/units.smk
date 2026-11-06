"""From raw records to one representative per reference unit."""


rule bcdm_records:
    input:
        BCDM.get("path") or [],
    output:
        f"{OUT}/units/records.tsv",
    params:
        src=SRC,
        marker=BCDM.get("marker", "COI-5P"),
        taxa=BCDM.get("taxa", {}),
        countries=BCDM.get("countries", []),
    log:
        f"{OUT}/logs/bcdm_records.log",
    conda:
        "../envs/timegraft.yaml"
    script:
        "../scripts/bcdm_records.py"


rule filter_records:
    input:
        RECORDS,
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
