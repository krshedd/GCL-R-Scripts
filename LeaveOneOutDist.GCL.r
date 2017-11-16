LeaveOneOutDist.GCL=function(sillyvec,loci,groupvec,pgroup=rep(1/max(groupvec),max(groupvec))){
  ####################################################################################################################################################################################################################################################################
  #
  # This function calculates a genetic likelihood for each individual belonging
  # to each population, based on the probabilitie of each individuals genotype 
  # originating from population allele frequencies (leave one out). These
  # individual to population probabilities are also rolled up into reporting
  # groups. This genetic likelihood profile is intended to show how well
  # populations / reporting groups are differentiated by the markers at hand.
  # This function can aid in marker selection.
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # sillyvec = Kodiak49Pops
  #   ~ A character vector of sillys (pops make more sense, but collecitons could be used).
  # loci = loci90
  #   ~ A character vector of loci. Can include haploid markers and/or combined loci.
  # groupvec = Kodiak49GroupVec9
  #   ~ A numberic vector of length(sillyvec) sepcifying group assignments.
  # pgroup = rep(1/max(groupvec),max(groupvec))
  #   ~ Default is for a regionally flat prior
  #     This means that each individual has an equal probability of coming from any group,
  #     otherwise, there is an implicit population flat prior, which would mean that groups
  #     with lots of populations would "pull" in more fish.
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Unnamed List of 2
  #  List of (number of RGs)
  #   List of (number of pops)
  #    numberic vector (length of number of individuals in pop) 
  #     probability of each individual originated from a pop in that group
  #  List of (number of pops)
  #   List of (number of pops)
  #    numberic vector (length of number of individuals in pop) 
  #     probability of each individual originated from a pop in that pop
  #
  # e.g. list[[1]][[1]][[1]][[1]] = probability of indv 1 in pop 1 originating from a pop in RG 1
  # e.g. list[[2]][[1]][[1]][[1]] = probability of indv 1 in pop 1 originating from  pop 1
  #
  # Created by Jim Jasper on unknown???
  # Modified by Kyle Shedd on unkown to include progress bars
  #
  ####################################################################################################################################################################################################################################################################
  
  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}
  
  cat("3 main function tasks: \n1) Determine qstar\n2) Determine likelihood profile for each population (longest task)\n3) Wrap up to Reporting Groups\n")
  if (.Platform$OS.type == "windows") flush.console()
  
  Y=FreqPop.GCL(sillyvec,loci)
  
  if(any(sapply(loci, function(locus) sum(colSums(Y[, locus, ]) > 0, na.rm = TRUE) ) == 1)) {
    
    fixed_loci <- which(sapply(loci, function(locus) sum(colSums(Y[, locus, ]) > 0, na.rm = TRUE) ) == 1)
    
    invisible(sapply(fixed_loci, function(i) {message(paste(loci[i], "is fixed for all 'sillyvec' and was removed"))}))
    
    loci <- loci[-fixed_loci]
    
  }  # Checking for fixed loci
  
  ORD=order(groupvec,decreasing=FALSE)
  
  sillyvec=sillyvec[ORD]
  
  groupvec=groupvec[ORD]
  
  C=as.vector(table(groupvec))
  
  cmC=c(0,cumsum(C))
  
  G=length(C)
  
  nsilly=length(sillyvec)
  
  p=setNames(rep(pgroup/C,C),sillyvec)
  
  nloci=length(loci)
  
  alleles=LocusControl$alleles[loci]
  
  nalleles=LocusControl$nalleles[loci]
  
  ploidy=LocusControl$ploidy[loci]
  
  my.gcl=lapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)})
  names(my.gcl)=sillyvec  
  
  n=sapply(my.gcl,function(gcl){gcl$n})
  names(n)=sillyvec
  
  counts=lapply(my.gcl,function(gcl){gcl$counts[,loci, , drop = FALSE]})
  names(counts)=sillyvec
  
  Y=FreqPop.GCL(sillyvec,loci)  # re-running 'Y' with potentially new 'loci' (fixed loci removed)
  
  N=sapply(loci,function(locus){sapply(sillyvec,function(silly){sum(Y[silly,locus,1:nalleles[locus]])})})
  dimnames(N)=list(sillyvec,loci)
  
  q = array(data = NA, dim = c(nsilly, nloci, max(nalleles)), dimnames = list(sillyvec, loci, paste("Allele ", 1:max(nalleles), sep = "")))
  for(silly in sillyvec){
    for(locus in loci){
      q[silly,locus,1:nalleles[locus]] = (Y[silly,locus,1:nalleles[locus]]+1/nalleles[locus]) / (N[silly,locus]+1)
    }
  } 
  
  cat("\nDetermine leave-one-out q (qstar) for each silly\n", sep = '')
  if (.Platform$OS.type == "windows") flush.console()
  pb <- txtProgressBar(min = 0, max = length(sillyvec), style = 3)
  
  qstar=lapply(sillyvec,function(silly){array(NA,c(n[silly],nloci,max(nalleles)),dimnames=list(seq(n[silly]),loci,paste("Allele ",1:max(nalleles),sep=""))) });names(qstar)=sillyvec 
  for(silly in sillyvec){
    setTxtProgressBar(pb = pb, value = which(sillyvec == silly))
    for(k in seq(n[silly])){
      for(locus in loci){
        for(allele in seq(nalleles[locus])){
          qstar[[silly]][k,locus,allele]=ifelse(is.na(counts[[silly]][k,locus,1]),NA,
                                                (Y[silly,locus,allele]-counts[[silly]][k,locus,allele]+1/nalleles[locus])/(N[silly,locus]-sum(counts[[silly]][k,locus,seq(nalleles[locus])])+1)
          )
        }
      }
    }
  }
  
  genofunc=function(x,theta,NAind){return(prod(sapply(loci[NAind],function(locus){2*dmultinom(x=x[locus,1:nalleles[locus]],size=ploidy[locus],prob=theta[locus,1:nalleles[locus]])} )))}
  
  
  T=lapply(sillyvec,function(silly1){TT=lapply(sillyvec,function(silly2){rep(NA,n[silly2])});names(TT)=sillyvec;TT});names(T)=sillyvec
  P=lapply(sillyvec,function(silly1){PP=lapply(sillyvec,function(silly2){rep(NA,n[silly2])});names(PP)=sillyvec;PP});names(P)=sillyvec
  
  cat("\nDetermine likelihood profile for each population\n", sep = '')
  if (.Platform$OS.type == "windows") flush.console()
  
  for(silly1 in sillyvec){
    setTxtProgressBar(pb = pb, value = which(sillyvec == silly1))
    for(m in seq(n[silly1])){
      x=array(data = counts[[silly1]][m,1:nloci,1:max(nalleles)], dim = c(nloci, max(nalleles)), dimnames = list(c(loci), paste("Allele", seq(max(nalleles)))))
      NAind=!is.na(x[1:nloci,1])
      theta=q[1:nsilly,1:nloci,1:max(nalleles), drop = FALSE]
      theta[silly1,1:nloci,1:max(nalleles)]=qstar[[silly1]][m,1:nloci,1:max(nalleles)]
      for(silly2 in sillyvec){
        P[[silly2]][[silly1]][m]=p[silly2]*genofunc(x=x,theta=array(data = theta[silly2,1:nloci,1:max(nalleles)], dim = c(nloci, max(nalleles)), dimnames = list(c(loci), paste("Allele", seq(max(nalleles))))), NAind)
      }  
      for(silly2 in sillyvec){
        T[[silly2]][[silly1]][m]=P[[silly2]][[silly1]][m]/sum(sapply(sillyvec,function(silly){P[[silly]][[silly1]][m]}))
      }  
    }
  }
  
  RT=lapply(seq(G),function(g){RTT=lapply(sillyvec,function(silly1){rep(NA,n[silly1])});names(RTT)=sillyvec;RTT});names(RT)=paste("Group" ,seq(G),sep="")
  
  cat("\nDetermine likelihood profile for each Reporting Group\n", sep = '')
  if (.Platform$OS.type == "windows") flush.console()
  
  for(silly1 in sillyvec){
    setTxtProgressBar(pb = pb, value = which(sillyvec == silly1))
    for(m in seq(n[silly1])){
      for(g in seq(G)){
        rt=0
        for(silly2 in sillyvec[(cmC[g]+1):cmC[g+1]]){
          rt=rt+T[[silly2]][[silly1]][m] 
        }
        RT[[g]][[silly1]][m]=rt
      }
    }
  }  
  
  return(list(RT,T))
  
}