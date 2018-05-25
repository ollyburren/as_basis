# read in a VCF file of genotypes and create a matrix for projection on to a basis.


bcftools_bin<-'~/bin/bcftools-1.4/bcftools'

parseVCF <- function(bcftools_bin,vcf,pos_file,samp_file){
  header_cmd<-sprintf("%s view --header-only -S %s -R %s  %s",bcftools_bin,samp_file,region_file,vf)
  my.pipe<-pipe(header_cmd)
  header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
  close(my.pipe)
  cnames<-unlist(strsplit(header,"\t"))
}
