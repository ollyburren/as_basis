# Systematic comparison of JIA using cupcake

## Summary statistics

`./R/analysis/distance_simulations.R` - This script needs to be adapted to run on the HPC and do a pairwise comparison between all disease subtypes.

Note that we set all OR to 1 so in fact we are simulating under the null with no shrinkage. Thus the plots generated make sense as ERA has smaller sample size so will have greater spread.

Note also that cc and ERA do share cases therefore these metrics are likely conservative.

From Chris

It so happens, that their correlation depends only on the sample size, and how much of that is shared.  if they share no samples, they are independent, and you can just simulate each, project, calculate difference, and repeat. if they share only controls, *and* the control size is large relative to the number of cases, they are still nearly independent, and we can do as above in fact, if they shared lots of samples and were not independent, this would make them *closer* than two independent samples (to see this, think about what happens as they share more and more cases and controls until the two datasets are identical and therefore have distance 0)

I guess we should implement the covariance-variance matrix to take into account sharing of cases. See [slack](https://wallacegroup.slack.com/archives/C3QDEP0V7/p1508760637000020). But results seem pretty convincing if we assume that both beta's are centred on zero (we only care about relative differences) then you never see such a large distance in 200 simulations.

I reimplemented the above so that it scales across all JIA subtypes. We can run a set of simulations and project onto a cached basis using the following:

`Rscript --vanilla /home/ob219/git/as_basis/R/analysis/null_by_sample_size.R -o /home/ob219/scratch/as_basis/jia_summ_analysis/ -c /home/ob219/tmp/as_basis_cache.RDS -t jia_sys`


* Write up simulation method.
* Add parameters for different subsets of data.
* Split final matrix creation into chunks to make more performant.

## Individual level data


`./R/jia_project_individuals_aligned_to_basis.R` - This script projects takes JIA genotype data that has been aligned to the basis and then computes posterior log(OR).

* Fold functions into cupcake. This will likely involve making `snpStats` or annotSnpStats a dependency.
* Scree plot to work out what is the optimum number of PC's to include in the t-test.
