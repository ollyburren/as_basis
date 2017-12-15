# December Notes

Decided to rewrite the pipeline so that we align to a reference panel.

## Step 1

This involves reprocessing raw summary statistics so that they have the format:

| pid | a1 | a2 | or | p.value |
|:---:|:--:|:--:|:--:|:-------:|
|1:123353|A|G|0.99|1e-3|

scripts to do this are
* ./align_OR/process_raw_GWAS.R - This script processes non BioBank summary summary_statistics
* ./align_OR/process_bb_summary_stats.R - This script processes BioBank summary statistics

The latter script has a dependency in that it only processes variants that will be included in the basis.

This produces a set of tab delim files with the above format in `/home/ob219/scratch/as_basis/gwas_stats/processed`

## Step 4

QC summary statistics. This involves:

1. Check summary statistics vs ImmunoBase annotations and which is the effect allele.

* ./align_OR/check_processed_vs_immunobase.R

2. Check for isolated associations and remove. Perhaps compute posterior probability of inclusion - generate snp set and flag those where this contains just one SNP. For these flagged  SNPs use LD to pull SNPs that are in LD (0.8) and examine their associations. If there are no snps manually examine.




## Step 3

1. First prepare reference. In this case due to input samples we use EUR 1KG reference dataset and prepare a file with the following

| pid	| a1 | a2 |	a1.af | ld.block |
|:---:|:--:|:--:|:-----:|:--------:|
|1:10583 |G |	A |	0.7929 | 1983  |

This involves:

- selecting variants that are common between basis traits

- removing variants (such as those in MHC)

- Adding ld block information

* ./align_OR/create_as_basis_snp_support.R

Here we select variants that are found in all the cohorts and remove MHC (chr6:20-40Mb)

4. aligning OR with a reference data set.

* ./align_OR/align_with_reference.R. This script does the alignment note that to create a shell script based on manifest file there is some code at the end. Note that this process is iterative. We do the alignment and then may notice additional missing SNPs we can feed these back into `create_as_basis_snp_support.R` and then do the realignment. When finished should have a directory of aligned,filtered variants for basis creation.

## Analysis

* ./alignOR/bb_plot.R This script creates the basis and then projects biobank summary statistics on.

* ./alignOR/jia_plot.R This script creates the basis and then projects jia subtype summary statistics on.

## Projecting individuals

* ./align_OR/jia_project_individual.R - This code takes `annotSnpStats` objects (it's much quicker if these have been reduced by prefiltering against the basis reference using `./jia_ss/cut_down_annoSnpStats_jia.R`) and then using reference we have previously generated it assigns log(OR) on an individual basis.

* ./align_OR/split_individual_lor.R To make computation easier (i.e. so we don't have a huge matrix) we split files down. This file also contains code to plot scree and compute one sample t-tests for individual principle components vs control. 
