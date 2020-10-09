include: "rules/common.smk"

rule all:
    input:
        expand(OUT_DIR+"{subset}_gistic_results/{subset}.txt",subset=SUBSETS)

include: "rules/gistic.smk"
