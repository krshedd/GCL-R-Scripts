SampSizeByLocus.GCL=function(sillyvec,loci){
##############################################################################################################################################################################
#
#  This function creates a data.frame of sample size by locus with nrow=length(sillyvec) and ncol=length(loci). 
#
#  "sillyvec" is a vector of silly names (without the ".gcl") you wish to know sample size by locus.
# 
#  "loci" is a vecor of locus names you want sample sizes for.
#  
#  Written by JJ 4/24/2011
##############################################################################################################################################################################

  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}

  nsilly=length(sillyvec)

  nloci=length(loci)

  nalleles=LocusControl$nalleles[loci]

  ploidy=LocusControl$ploidy[loci]

  results=array(NA,c(nsilly,nloci),dimnames=list(sillyvec,loci) )

  for(silly in sillyvec){
    counts=get(paste(silly,".gcl",sep=""),pos=1)$counts[, loci,, drop = FALSE]
    n=nrow(counts)
    if(nloci>1){   
      results[silly,loci]=sapply(loci,function(locus){sum(counts[!is.na(counts[1:n,locus,1]),locus,1:nalleles[locus]])/ploidy[locus]})
    }
    if(nloci==1){   
      results[silly,loci]=sum(counts[!is.na(counts[1:n,1]),1:nalleles[loci]])/ploidy[loci]
    }    
  }
  return(data.frame(results))
}
