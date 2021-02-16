# COVGAP
#### COVID-19 Genome Analysis Pipeline, for whole genome sequencing and annotation of SARS-CoV-2
## Introduction:
COVGAP will process the raw demultiplexed short reads (Illumina) produced by an amplicon sequencing approach targeting SARS-CoV-2 (e.g. ARTIC)
COVGAP will quality filter your reads and map them against the reference SARS-CoV-s genome, providing the variants in vcf format, and include them in a draft consensus genome. Only the variants showing high coverage and high alternative frequency (primary alleles variants) are included in the consensus. 
Alleles below the coverage or frequency thresholds (secondary alleles variants), and those alleles occurring as alternative to primary allele variants (alternative allele variants) are stored and annotated separately for the user consultation, but will NOT be part of the consensus.
Loci showing not enough coverage to allow a confident variant call are masked with N, instead of calling the reference. The consensi are labelled according to a user defined threshold of ambiguous base count (N) into LowCov (Ns above threshold), HighCov (Ns below threshold).

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
The initial reads need to be demultiplexed, paired end and tagged with `_R1.fastq.gz` and `_R2.fastq.gz`. To run the command with default options simply type:
```
snakemake -s path/to/covgap/covgap.smk --use-conda --cores 4
```
parameters are customisable by adding the flag `--config` followed by the parameter to be changed. A full list of the parameters follows:

| Parameter  | Description  | Default |
| :--------------- |:---------------------------| :-------|
| read_dir      | Directory containing the demultiplexed raw reads | reads/ |
| QC_sliding_windows      | Sliding window used to scan the 5‟ end. It clips the read once the average quality within the window falls below a threshold         |   4 |
| QC_phred_score | Minimum quality threshold per base (phred score) used to trim the read (see above sliding window) | 20 |
| adapters | The sequence of adapters to be trimmed (fasta format) | Nextera Flex adapters |
| primers | The sequences of forward and reverse primers used for the amplicon sequence (fasta format) | ARTIC primers V3 |
| ref_genome | The full SARS-CoV-2 genome used as reference | Wuhan-Hu-1 (NC_045512.2) |
| mapr | Mapping criteria to consider a read as mapping (json format) | map, mate_mapped |
| unmapr | Unmapping criteria to consider a read as NOT mapping (json format)   | unmapped, mate_unmapped |
| uptrim_threshold | Upperbound threshold to prune the read pileup if above the threshold. NOTE: this trimming is used for plot display only. It is NOT affecting the variant call | 1000 |
| variant_freq | Minimum alternative frequency in the read pileup to consider a variant as primary. We reccommend to choose values higher than 0.7 NOTE: they will be considered primary only those variants showing a variant_freq AND a variant_depth above threshold | 0.7 |
| variant_depth | Minimum locus depth in the read pile up to consider a variant as primary. NOTE: they will be considered primary only those variants showing a variant_freq AND a variant_depth above threshold | 50 |
| n_threshold | Ambiguous base call percentage overall the consensus to consider the genome HighCov or LowCov | 10 |


