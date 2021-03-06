# COVGAP
#### COVID-19 Genome Analysis Pipeline, for whole genome sequencing and annotation of SARS-CoV-2 {basic version, more on the way}
## Introduction:
COVGAP will process the raw demultiplexed short reads (Illumina) produced by an amplicon sequencing approach targeting SARS-CoV-2 (e.g. ARTIC)
COVGAP will quality filter your reads and map them against the reference SARS-CoV-2 genome, providing the variants in vcf format, and include them in a draft consensus genome. Only the variants showing high coverage and high alternative frequency (primary alleles variants) are included in the consensus. 
Alleles below the coverage or frequency thresholds (secondary alleles variants), and those alleles occurring as alternative to primary allele variants (alternative allele variants) are stored and annotated separately for the user consultation, but will NOT be part of the consensus.
Loci showing not enough coverage to allow a confident variant call are masked with N, instead of calling the reference. The consensi are labelled according to a user defined threshold of ambiguous base count (N) into LowCov (Ns above threshold), HighCov (Ns below threshold).

## Installation:
### Requirements:
- git >= 1.8.3
- python >= 3.7
- conda >= 4.9 (for mac only)
- snakemake >= 5.26

In case you don't have snakemake already, we recommend to install it via mamba:

```
$ conda install -c conda-forge mamba
```
then to install the full package run:

```
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
more options can be found at the [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) page. 
Finally source the environment:

```
$ conda activate snakemake
```

### Installation:
Portability expansion is on its way, for the moment, we reccommend to clone the repository in a directory of choice:
```
$ git clone https://github.com/appliedmicrobiologyresearch/covgap/ path/to/covgap
```

## Using covgap: 
The initial reads need to be already demultiplexed, paired end and tagged with `_R1.fastq.gz` and `_R2.fastq.gz`. Anything prior this tagging will be interpreted as the sample name. For the moment COVGAP can run only in linux or macOS. Choose the appropriate version while running snakemake. To run the command with default options simply type:
#### For linux
```
$ snakemake -s path/to/covgap/repo/covgap_linux.smk --use-conda --cores [number of cores reserved]
```
#### For Mac
```
$ snakemake -s path/to/covgap/repo/covgap_mac.smk --use-conda --cores [number of cores reserved]
```
## Parametrisation and options
Parameters are customisable by adding the flag `--config` followed by the parameter to be changed. For example, to change the directory containing reads from the default to `../my_reads/` , the command will be:
```
$ snakemake -s path/to/covgap/repo/covgap_mac.smk --use-conda --cores 4 --config read_dir=../my_reads/
```

A full list of the parameters follows:

| Parameter  | Description  | Default |
| :--------------- |:---------------------------| :-------|
| read_dir      | [path] Directory containing the demultiplexed raw reads | reads/ |
| QC_sliding_windows | [numeric] Sliding window used to scan the 5‟ end. It clips the read once the average quality within the window falls below a threshold         |   4 |
| QC_phred_score | [numeric] Minimum quality threshold per base (phred score) used to trim the read (see above sliding window) | 20 |
| adapters | [file.fasta] The sequence of adapters to be trimmed | Nextera Flex adapters |
| primers | [file.fasta] The sequences of forward and reverse primers used for the amplicon sequence | ARTIC primers V3 |
| ref_genome | [file.fasta] The full SARS-CoV-2 genome used as reference, NOTE: if using a custom reference, make sure it is indexed before the run | Wuhan-Hu-1 (NC_045512.2) |
| mapr | [file.json] Mapping criteria to consider a read as mapping | map, mate_mapped |
| unmapr | [file.json] Unmapping criteria to consider a read as NOT mapping | unmapped, mate_unmapped |
| uptrim_threshold | [numeric] Upperbound threshold to prune the read pileup if above the threshold. NOTE: this trimming is used for plot display only. It is NOT affecting the variant call | 1000 |
| variant_freq | [numeric] Minimum alternative frequency in the read pileup to consider a variant as primary. We reccommend to choose values higher than 0.7 NOTE: they will be considered primary only those variants showing a variant_freq AND a variant_depth above threshold | 0.7 |
| variant_depth | [numeric] Minimum locus depth in the read pile up to consider a variant as primary. NOTE: they will be considered primary only those variants showing a variant_freq AND a variant_depth above threshold | 50 |
| n_threshold | [numeric] Ambiguous base call percentage overall the consensus to consider the genome HighCov or LowCov | 10 |

## Output files:

| File | Description |
| :----------------- |:---------------------------| 
| result/[sample]/Mapping/[sample].alignment.removed.duplicates.mapped.reads.only.sorted.[High/Low]Cov.bam | final alignment containing only reads mapping the reference |
| result/[sample]/Mapping/[sample].alignment.stats.tab | stats containing the amount of mapped vs unmapped reads |
| result/[sample]/Mapping/[sample].average.coverage.tab | stats containing the average coverage for the sample and its standard deviation |
| result/[sample]/variantcall/[sample].final.vcf | the final vcf including all the primary alleles being found |
| result/[sample]/variantcall/[sample].minority_alleles.vcf | the final vcf including all the minority alleles being found |
| result/[sample]/variantcall/[sample].minority_alleles_report.tech | additional report on minority alleles, featuring the denomination, the frequency score, as well as depth |
| result/[sample]/variantcall/[sample].consensus.[High/Low]Cov.fasta | the consensus genome featuring all the primary allele variants |
| result/[sample]/variantcall/[sample].classification.tab | a statement file including the classification of the sample and its percentage of ambiguous basecalls |


## Contact:

Alfredo Mari: alfredo.mari@unibas.ch

## Cite COVGAP:

Please refer to the following publications to properly cite COVGAP:

1. Stange M, Mari A, Roloff T, Seth-Smith HM, Schweitzer M, Brunner M, et al. (2021) SARS-CoV-2 outbreak in a tri-national urban area is dominated by a B.1 lineage variant linked to a mass gathering event. PLoS Pathog 17(3): e1009374. https://doi.org/10.1371/journal.ppat.1009374
2. Mari A, Roloff T, Stange M, Søgaard KK, Asllanaj E, Tauriello G, Alexander LT, Schweitzer M, Leuzinger K, Gensch A, Martinez AE, Bielicki J, Pargger H, Siegemund M, Nickel CH, Bingisser R, Osthoff M, Bassetti S, Sendi P, Battegay M, Marzolini C, Seth-Smith HMB, Schwede T, Hirsch HH, Egli A. Global Genomic Analysis of SARS-CoV-2 RNA Dependent RNA Polymerase Evolution and Antiviral Drug Resistance. Microorganisms. 2021; 9(5):1094. https://doi.org/10.3390/microorganisms9051094





