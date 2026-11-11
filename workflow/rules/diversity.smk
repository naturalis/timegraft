"""Phylogenetic diversity per sample and between samples."""


rule diversity:
    input:
        tree=rules.place_asvs.output.tree if USE_ASVS else rules.place_units.output.tree,
        samples=rules.place_asvs.output.samples if USE_ASVS else (config.get("samples") or []),
        provenance=rules.graft.output.provenance,
        subtrees=expand(f"{OUT}/subtrees/shard_{{i}}.tsv", i=SHARDS),
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
