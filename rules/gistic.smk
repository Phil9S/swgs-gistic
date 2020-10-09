rule gistic:
    input:
        list=get_list
    output:
        OUT_DIR+"{subset}_gistic_results/{subset}.txt"
    script:
        "../scripts/focal_CNA_analysis.R"
