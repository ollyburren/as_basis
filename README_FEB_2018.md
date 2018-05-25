# Instructions and notes Feb 2018

# Process raw summary stats
  - `process_astle_raw_gwas.R` Script computes dichotomous OR from cont trait.
  - `process_bb_summary_stats.R` Script processes BioBank summary stats
  - `process_other_gwas.R` Script copies old data see `process_raw_GWAS.R` for further information

# What SNPs  in common ?
  - `perl ../perl/variants_in_common.pl | bash` parses all files included in a manifest and works out what SNPs are in common across all studies. TODO - make output file (currently ~/tmp/intersect.txt) and input manifest file parameters.
  - `create__as_basis_snp_support.R` using 1K genomes reference dataset and LD block assignments create a support file. (Perhaps use LD blocks from [`Berisa,T. and Pickrell,J.K. (2016) Approximately independent linkage disequilibrium blocks in human populations. Bioinformatics, 32, 283–285.`] - will require some callibration).
  - `filter_gwas.R` takes a snp_support_file generated above and filters a GWAS to obtain source files for basis analysis.
  - `align_with_reference.R` **RUN TWICE** ! contains code to create a shell script to create neccessary files for using in the basis creation (TODO perhaps link `filter_gwas.R` functionality.). Note because although initial filtering finds SNPs in common sometimes the alleles don't match with the refefence. This script will remove those where the alleles are different. Until we have run the whole manifest we don't know what this exclusion list will be so we need to run this twice.

  - UPDATE basis support ?

# Create Basis.

If we leave out MS get 315,072 shared SNPs for the basis if we include we loose 9% and therefore get 286,915. How does the loss affect the basis ?
