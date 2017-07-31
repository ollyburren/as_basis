# as_basis Using covariance to leverage rare disease GWAS

## Outline

Some diseases are likely polygenic but rare. This makes convential GWAS approaches underpowered, here we investigate PCA based methods to analyse the relationships between larger GWAS in related traits focussing on autoimmunity to see if we can use these to better understand the biology of rarer phenotypes.

## R :: This contains most R scripts used to develop the method

- computeHaplotypeFreq1KG.R - Use Mary Fortunes simGWAS to simulate GWAS where we know what the causal variant relatedness is between studies.

- process_raw_summ_stats.R - Code used to process GWAS summary stats into format to build basis

- wakefield.R -  Code used to compute posterior probabilities to weight SNPs as to whether they are causal or not

- central_functions.R - Functions that are used in other scripts

- empirical_prior.R - Estimate q for each SNP which is taken from the posterior probability of inclusion for any trait.

- computeBetaCovariates.R - This is to speed up simulations for examining how sample size effects things. We can compute beta covariates and then rescale to simulate different sample sizes in the projections.

- simulate_AI_weightings.R - This computes simulations for overall weightings. This is a considerably more complicated simulation as weightings have to be recalculated for each basis and projection simulation pair. This should be run with

```bash
OUTDIR=
RSCRIPT=./simulate_AI_weightings.R
NO_SIMS=100
NO_SIMS_PER_JOB=10
TOTAL=$((NO_SIMS / NO_SIMS_PER_JOB))
COUNTER=1

while true; do
	if [[ "$COUNTER" -gt "$TOTAL" ]]; then
		exit 1
	fi
	outfile=${OUTDIR}$COUNTER.RDS
	if [[ ! -e "$outfile" ]]; then
		echo "Rscript --vanilla $RSCRIPT --sim_size $NO_SIMS_PER_JOB --out ${OUTDIR}$COUNTER.RDS"
	fi
	COUNTER=$((COUNTER+1))
done
```

- simulate_AI_weightings.R - This is parallel to *simulate_AI_weightings.R* but we simulate from aff.t1d and project these onto basis constructed using simulations from ill.t1d. Use similar BASH script as *simulate_AI_weightings.R* to run.


## R/figures

This dir contains code for creating figures and misc analysis

- analysis1_w_QT.R - Early analysis to look at the effect of different scaling metrics on PCA eucledian distances between traits. We create a basis using the scaling and then project on t1d.affy - we look to see what metrics minimise the distance across all PC to t1d.ill. This was early work and proof of principal that included quantitative traits.

- analysis1.R - This is the same as above but with QT removed and extra category for non scaled projection which appears to improve things. This stores metrics that are used for all further figures. **NOTE: If the basis changes then so should this**.

- analysis1_w_sim.R - This is the same as above but includes a simulation mechanism to examine the effect of 100 random basis generated with ill.t1d and projection of 100 randomised GWAS (simulated from ill.t1d dataset so related). This allows us to calibrate Eucledian distances and see how well different scaling metrics perform.

- aff_ill_w_sim.R - Here we simulate 100 aff.t1d and then project these onto  100 bases created with ill.t1d (simulated 100 times).

- analysis2_w_sim.R - This is the same as above but we simulate under no sharing - we hope that for good metrics discovered above there is no overlap between no sharing and sharing confidence intervals.

- analysis3_w_sim.R - Here we examine what happens as we alter sample size for projections - this code runs on the hpc it uses betaCovariates computed above.

`Rscript --vanilla analysis3_w_sim.R -o /tmp/sample_size/3_3.RData -c 3000  -t 3000`
