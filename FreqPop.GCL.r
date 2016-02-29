FreqPop.GCL=function(sillyvec,loci){

##############################################################################################################################################################################
#
#  This function creates an array of allele frequencies
#
# Written by JJ 6/3/2011
##############################################################################################################################################################################

  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}

  nsilly=length(sillyvec)

  nloci=length(loci)

  alleles=LocusControl$alleles[loci]

  nalleles=LocusControl$nalleles[loci]

  ploidy=LocusControl$ploidy[loci]

  Freq=array(NA,c(nsilly,nloci,max(nalleles)),dimnames=list(sillyvec,loci,paste("Allele ",1:max(nalleles),sep="")))
  
  for(silly in sillyvec){
    my.gcl=get(paste(silly,".gcl",sep=""),pos=1)
    scores=my.gcl$scores
    for(locus in loci){
      myscores=scores[,locus,1:ploidy[locus]][!is.na(scores[,locus,1:ploidy[locus]])]
      allelevec=paste("Allele ",1:nalleles[locus],sep="")
      for(allele in 1:nalleles[locus]){
        Freq[silly,locus,allelevec[allele]]=sum(myscores==alleles[[locus]][allele])
      }
    }
  }
  return(Freq)
}
