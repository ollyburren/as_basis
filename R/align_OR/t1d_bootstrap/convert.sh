IN_DIR='/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-imputed/tmp/'
OUT_DIR='/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-imputed/tmp/snpStats'
SAMPLE_FILE_DIR='/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-imputed/'
RSCRIPT='Rscript /home/ob219/git/as_basis/R/align_OR/t1d_bootstrap/convert_t1d_imputed.R' 

for i in `\ls ${IN_DIR}*.gz`; do
	echo "$RSCRIPT -i $i -o $OUT_DIR -s $SAMPLE_FILE_DIR" 
done
