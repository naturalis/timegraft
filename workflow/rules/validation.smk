"""Validation: how far each barcode distance tracks a trusted reference tree."""


rule validate:
    input:
        reference=config.get("validation", {}).get("reference_tree") or [],
        reps=rules.representatives.output.reps,
        records=rules.filter_records.output.kept,
    output:
        pairs=f"{OUT}/validation/pairs.tsv",
        summary=f"{OUT}/validation/summary.tsv",
        bins=f"{OUT}/validation/binned_medians.tsv",
        plot=f"{OUT}/validation/saturation.png",
    params:
        src=SRC,
        settings=config.get("validation", {}),
    log:
        f"{OUT}/logs/validate.log",
    conda:
        "../envs/timegraft.yaml"
    script:
        "../scripts/validate.py"
