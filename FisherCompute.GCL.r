FisherCompute.GCL <- function(freq, loci,prec=4){
	
  nloci <- length(loci)

  nalleles=LocusControl$nalleles[loci]

  pval=rep(NA,nloci)
  names(pval)=loci

  for(locus in loci){
    pval[locus]=fisher.test(freq[,locus,1:nalleles[locus]],workspace=8000000,hybrid=TRUE)$p.value
  }
  pval[pval==0]=min(pval[pval!=0])
  return(data.frame(Locus=c(loci,"Overall"),pval=round(c(as.numeric(pval),pchisq(q=-2*sum(log(pval)), df=2*nloci,lower.tail=F)),prec)))  
}



