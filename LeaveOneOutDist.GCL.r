LeaveOneOutDist.GCL=function(sillyvec,loci,groupvec,pgroup=rep(1/max(groupvec),max(groupvec))){

  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}

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

  counts=lapply(my.gcl,function(gcl){gcl$counts[,loci,]})
  names(counts)=sillyvec
 
  Y=FreqPop.GCL(sillyvec,loci)

  N=sapply(loci,function(locus){sapply(sillyvec,function(silly){sum(Y[silly,locus,1:nalleles[locus]])})})
  dimnames(N)=list(sillyvec,loci)

  q=array(NA,c(nsilly,nloci,max(nalleles)),dimnames=list(sillyvec,loci,paste("Allele ",1:max(nalleles),sep="")))
  for(silly in sillyvec){
    for(locus in loci){
      q[silly,locus,1:nalleles[locus]]=(Y[silly,locus,1:nalleles[locus]]+1/nalleles[locus])/(N[silly,locus]+1)
    }
  } 

  qstar=lapply(sillyvec,function(silly){array(NA,c(n[silly],nloci,max(nalleles)),dimnames=list(seq(n[silly]),loci,paste("Allele ",1:max(nalleles),sep=""))) });names(qstar)=sillyvec 
  for(silly in sillyvec){
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

  genofunc=function(x,theta,NAind){return(prod(sapply(loci[NAind],function(locus){2*dmultinom(x=x[locus,1:nalleles[locus]],size=ploidy[locus],prob=theta[locus,1:nalleles[locus]])})))}
  

  T=lapply(sillyvec,function(silly1){TT=lapply(sillyvec,function(silly2){rep(NA,n[silly2])});names(TT)=sillyvec;TT});names(T)=sillyvec
  P=lapply(sillyvec,function(silly1){PP=lapply(sillyvec,function(silly2){rep(NA,n[silly2])});names(PP)=sillyvec;PP});names(P)=sillyvec
  for(silly1 in sillyvec){
    for(m in seq(n[silly1])){
      x=counts[[silly1]][m,1:nloci,1:max(nalleles)]
      NAind=!is.na(x[1:nloci,1])
      theta=q[1:nsilly,1:nloci,1:max(nalleles)]
      theta[silly1,1:nloci,1:max(nalleles)]=qstar[[silly1]][m,1:nloci,1:max(nalleles)]
      for(silly2 in sillyvec){
        P[[silly2]][[silly1]][m]=p[silly2]*genofunc(x=x,theta=theta[silly2,1:nloci,1:max(nalleles)],NAind)
      }  
      for(silly2 in sillyvec){
        T[[silly2]][[silly1]][m]=P[[silly2]][[silly1]][m]/sum(sapply(sillyvec,function(silly){P[[silly]][[silly1]][m]}))
      }  
    }
  }

  RT=lapply(seq(G),function(g){RTT=lapply(sillyvec,function(silly1){rep(NA,n[silly1])});names(RTT)=sillyvec;RTT});names(RT)=paste("Group" ,seq(G),sep="")
  for(silly1 in sillyvec){
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
