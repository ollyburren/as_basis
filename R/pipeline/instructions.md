#1_merge_all.GWAS.R

This R script reads in a directory of GWAS results

1. works out which variants are in all studies, prunes accordingly
2. adds allele frequencies from a reference dataset (e.g. 1K genome)
3. assigns variants to LD blocks based on 1cM blocks derived from HapMap project.
4. compute posterior probabilities for inclusion in causal set (using [Wakefields](http://faculty.washington.edu/jonno/papers/GE08.pdf) approximate Bayes Factors)

## Inputs

- directory containing a list of GWAS summary statistics in the following format saved as an RData [data.table object](https://cran.r-project.org/web/packages/data.table/). The filename should reference the CSV file below. For example, mystudy.RData will have a row in mystudy in the study configuration file below

| id          | chr | position | risk.allele | other.allele |   or | p.val |
|-------------|----:|---------:|-------------|--------------|-----:|------:|
| chr1:751343 |   1 |   751343 | T           | A            | 1.18 |  0.01 |
| chr1:751756 |   1 |   751756 | T           | C            | 1.17 |  0.01 |
| rs3094315   |   1 |   752566 | A           | G            | 1.14 |  0.01 |
| rs3131972   |   1 |   752721 | G           | A            | 1.14 |  0.01 |
| rs3131971   |   1 |   752894 | C           | T            | 1.15 |  0.01 |
| chr1:753405 |   1 |   753405 | A           | C            | 1.17 |  0.01 |

Note we expect all the OR to be with respect to the risk allele (i.e > 1)

- CSV file containing study configuration options (file)
"disease","cases","controls","pmid"
"mystudy",1234,4567,"1234"

- Reference allele frequency file (mine comes from EUR 1KG) as an RData data.table object.

| name        | chr | position | a1 | a2 | a2.f |
|-------------|-----|---------:|----|----|-----:|
| rs58108140  | 1   |    10583 | G  | A  | 0.21 |
| rs189107123 | 1   |    10611 | C  | G  | 0.02 |
| rs180734498 | 1   |    13302 | C  | T  | 0.14 |
| rs144762171 | 1   |    13327 | G  | C  | 0.04 |
| rs151276478 | 1   |    13980 | T  | C  | 0.02 |
| rs140337953 | 1   |    30923 | G  | T  | 0.73 |

- A file of tab delimited LD blocks.

| chr| start-end     |
|---:|---------------------|
| 18 | 74080457-74331870   |
|  3 | 132059281-133226149 |
| 16 | 1997615-2394258     |
| 14 | 51263446-51769865   |
| 20 | 37012011-38436593   |


## Outputs

- An intermediate .RData file containing a data.table as follows for all shared SNPs. This file is useful for adding further traits. By default ends up in `OUTDIR/shared_support_file_with_AF.RData`.

| chr | position | name      | minor\_allele |  maf |
|-----|---------:|-----------|---------------|-----:|
| 1   |  2033256 | rs908742  | A             | 0.35 |
| 1   |  2040936 | rs4648808 | T             | 0.11 |
| 1   |  2058023 | rs3128291 | A             | 0.11 |
| 1   |  2068906 | rs3128296 | G             | 0.11 |
| 1   |  2071340 | rs424079  | C             | 0.38 |
| 1   |  2082489 | rs262669  | A             | 0.38 |

- A rather large .RData output file containing all the information required for basis creation and projection.

## To Do

This currently works on all traits at the same time which means that it does not scale to a HPC context. Need to rewrite so that given the intermediate file we can parallelise by trait. One idea is to have a specific function that works out which SNPs are shared across traits as this is one of the main drivers for processing all traits at once. If we also added REF and ALT alleles to the intermediate file then we would be assured aligned OR.
