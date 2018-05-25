#Notes on bootstrap for T1D

We used imputed wtccc data `/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-imputed`.


This was filtered so that we only had SNPs that were in the basis I used the program [gtools](http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html) to filter see `run_gtools.R` for more details.

I then used snpStats to filter and convert these into snpStats objects see `convert_t1d_imputed.R` for more details.

`prep_bootstrap_files.R`  takes the filtered snpStas generated above and creates support bootstrap files. Some code need to be run once (1-93) - this creates basis and shrinkage files and checks allele alignment and case support files. The latter lines need to be run for each bootstrap case control configuration - fair bit of hardcoding here but OK as will only be run a few times.

`run_boostrap.R` this runs on the Q carries out GWAS and projects the results onto the basis with correct shrinkage.

Simple shell script for running the above

```
IN_DIR=/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/sims/wtccc_1000_1000_500/
RSCRIPT="Rscript /home/ob219/git/as_basis/R/align_OR/t1d_bootstrap/run_bootstrap.R"

for i in `\ls ${IN_DIR}bs_*`;do
echo "$RSCRIPT -m $i"
done
```
