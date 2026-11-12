"""Dated genus-level backbone from Open Tree of Life and TimeTree."""

BB = config.get("backbone_build", {})
STATIC = BB.get("source", "static") == "static"
OTOL_DIR = BB.get("otol_dir", "resources/otol")
SYNTH = f"{OTOL_DIR}/{BB.get('synthesis', 'opentree16.1')}/labelled_supertree.tre"
TAXONOMY = f"{OTOL_DIR}/{BB.get('taxonomy', 'ott3.7.3')}/taxonomy.tsv"
# support of each synthesis node by source trees, recorded in node_ages.tsv (#47)
SUPPORT = STATIC and BB.get("synthesis_support", True)
ANNOTATIONS = f"{OTOL_DIR}/{BB.get('synthesis', 'opentree16.1')}/annotations.json"


rule download_otol:
    output:
        synth=SYNTH,
        taxonomy=TAXONOMY,
    params:
        synthesis=BB.get("synthesis", "opentree16.1"),
        taxonomy=BB.get("taxonomy", "ott3.7.3"),
        directory=lambda w, output: str(Path(output.taxonomy).parent.parent),
    log:
        f"{OUT}/logs/download_otol.log",
    conda:
        "../envs/download.yaml"
    shell:
        """
        (
        curl -sSfL -o {output.synth} \
          https://files.opentreeoflife.org/synthesis/{params.synthesis}/{params.synthesis}_tree/labelled_supertree/labelled_supertree.tre
        curl -sSfL https://files.opentreeoflife.org/ott/{params.taxonomy}/{params.taxonomy}.tgz \
          | tar xzf - -C {params.directory} {params.taxonomy}/taxonomy.tsv
        ) > {log} 2>&1
        """


rule download_otol_annotations:
    output:
        ANNOTATIONS,
    params:
        synthesis=BB.get("synthesis", "opentree16.1"),
    log:
        f"{OUT}/logs/download_otol_annotations.log",
    conda:
        "../envs/download.yaml"
    shell:
        """
        curl -sSfL -o {output} \
          https://files.opentreeoflife.org/synthesis/{params.synthesis}/{params.synthesis}_tree/annotations.json \
          > {log} 2>&1
        """


rule build_backbone:
    input:
        reps=rules.representatives.output.reps,
        records=rules.filter_records.output.kept,
        synth=SYNTH if STATIC else [],
        taxonomy=TAXONOMY if STATIC else [],
        annotations=ANNOTATIONS if SUPPORT else [],
    output:
        tree=f"{OUT}/backbone/backbone.tre",
        genera=f"{OUT}/backbone/genera.tsv",
        nodes=f"{OUT}/backbone/node_ages.tsv",
        crowns=f"{OUT}/backbone/crown_ages.tsv",
        metadata=f"{OUT}/backbone/metadata.json",
        col_species=f"{OUT}/backbone/col_species.tsv",
    params:
        src=SRC,
        settings=BB,
    log:
        f"{OUT}/logs/build_backbone.log",
    conda:
        "../envs/timegraft.yaml"
    script:
        "../scripts/build_backbone.py"
