from snakemake.utils import validate
import pandas as pd
from pathlib import Path
from snakemake.utils import min_version
import os.path

min_version("5.10.0")

configfile: "config/config.yaml"

##### load config and sample sheets #####
samplesheet = pd.read_table(config["samplesheet"]).set_index(["subset"], drop=False)
SUBSETS = list(samplesheet['subset'])

def get_list(wildcards):
    files = list(samplesheet.loc[(wildcards.subset), ["file"]])
    return files

OUT_DIR=config["output_dir"]
