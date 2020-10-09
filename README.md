# Snakemake workflow: swgs-gistic

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg)](https://snakemake.bitbucket.io)  

## Authors

* Philip Smith (@phil9s)

## Usage

Generates GISTIC 2.0 results from shallow WGS-derived absolute copy number data using various sample subsets.

### Step 1: Obtain a copy of this workflow

2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Install GISTIC 2.0

Set up directory for GISTIC 2.0 in the repository

```
cd swgs-gistic/
```

Download and set up GISTIC 2.0

```
wget -c ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTIC_2_0_23.tar.gz
tar zxf GISTIC_2_0_23.tar.gz
cd MCR_Installer/
unzip MCRInstaller.zip
./install -mode silent -agreeToLicense yes -destinationFolder {FULL_PATH}/swgs-gistic/MATLAB_Compiler_Runtime/
cd ../
```
_note: {FULL_PATH} should be replaced with the absolute path for the swgs-gistic pipeline repository_

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

