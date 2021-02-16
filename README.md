# COVGAP
#### COVID-19 Genome Analysis Pipeline, for whole genome sequencing and annotation of SARS-CoV-2
## Introduction:
SARS
Basic version: new versions will be uploaded within short time.

## Installation:
### Requirements:
- pyhton >= 3.7
- snakemake >= 5.26
- git >=1.8.3

You can install snakemake by creating a conda environment, as illustrated [here] (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html):

### Installation:
Portability expansion is on its way, for the moment, we reccommend to clone the repository in a directory of your choice. The procedure was tested successfully in Linux and Mac operative systems.
```
$ git clone https://github.com/appliedmicrobiologyresearch/covgap/ /path/to/covgap
```

## Using covgap:
This procedure works with raw demultiplexed short reads (Illumina) produced by an amplicon sequencing approach targeting SARS-CoV-2 (e.g. ARTIC)
Covgap will quality filter your reads and map them against the reference SARS-CoV-s genome, providing the variants in vcf format, annotate them, and include them in the consensus genome.




