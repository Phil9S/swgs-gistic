include: "rules/common.smk"

rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.

include: "rules/gistic.smk"
