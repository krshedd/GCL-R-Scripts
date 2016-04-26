LeaveOneOutLikeProfile.GCL=function(popvec,loci,groupvec,groupnames,groupcomps=NULL,ncores=2){
#
#  rm(list=ls(all=TRUE))
#
#  load("V:\\WORK\\Chinook\\Cook Inlet\\Kenai River\\2011 Baseline\\KenaiChinook2011.RData")
#
#  popvec=KenaiNames10;loci=loci40;groupvec=groupvec3;groupnames=groups;groupcomps=NULL;ncores=8
#
#  likeprof=LeaveOneOutLikeProfile.GCL(popvec=KenaiNames10,loci=loci40,groupvec=groupvec3,groupnames=groups,groupcomps=NULL,ncores=8)
#
#  rm(list=ls(all=TRUE))
#
#  load("V:\\DATA\\WORK\\CHINOOK\\SE Alaska\\GAPS3.0\\GAPSbaseline_Feb2013.RData")
#
#  popvec=SEAKPops356;loci=GAPSloci;groupvec=GroupVec25RG;groupnames=Groups25;groupcomps=NULL;ncores=8
#
#  likeprof=LeaveOneOutLikeProfile.GCL(popvec=SEAKPops356,loci=GAPSloci,groupvec=GroupVec25RG,groupnames=Groups25,groupcomps=NULL,ncores=8)
#    
#

  begin.time=Sys.time()

  if(sum(!loci%in%LocusControl$locusnames)){stop(paste0("'",loci[!loci%in%LocusControl$locusnames],"' from argument 'loci' not found in 'LocusControl' object!!!"))}

  while(!require(abind)){install.packages("abind")}

  while(!require(foreach)){install.packages("foreach")}

  while(!require(doParallel)){install.packages("doParallel")}

  cl=makePSOCKcluster(ncores)

  registerDoParallel(cl,cores=ncores)  

  G=max(groupvec)

  names(groupvec)=popvec

  groupnamevec=setNames(groupnames[groupvec],popvec)

  if(is.null(groupcomps)){groupcomps=rep(1,max(groupvec))}

  p=setNames((groupcomps/tabulate(groupvec,G))[groupvec],popvec)

  nalleles=LocusControl$nalleles[loci]

  gclvec=sapply(popvec,paste0,sep=".gcl")

  gcls=sapply(gclvec,get,pos=1,simplify=FALSE)

  vials=sapply(gcls,function(gcl){rownames(gcl$counts)})

  atts=Reduce(rbind,lapply(popvec,function(pop){data.frame(FromGroup=as.vector(groupnamevec[pop]),FromPop=pop,Pop_Vial=paste(pop,vials[[pop]],sep="_"),stringsAsFactors=FALSE)}))

  rownames(atts)=atts$Pop_Vial

  x=Reduce(f=function(...){abind(...,along=1)},lapply(gcls,function(gcl){gcl$counts[,loci,]}))

  x <- x[, , colSums(is.na(x), dims = 2) < prod(dim(x)[1:2])]  # added by KS on Thu Dec 31 13:11:51 2015 to remove alleles that do not exist (NA) for ANY loci, this prevents issues with non-conformable arrays later on because y reduces array to only alleles that exist for loci
  
  x[is.na(x)]=0

  rownames(x)=atts$Pop_Vial

  y=FreqPop.GCL(popvec,loci)

  beta=array(0,dim=dim(y),dimnames=dimnames(y))

  for(locus in loci){

   beta[,locus,seq(nalleles[locus])]=1/nalleles[locus] 

  }#locus

  lnqFunc=function(Y){

     Y0=Y

     Y0[is.na(Y)]=0

     N=array(apply(Y0,1:2,sum),dim=dim(y),dimnames=dimnames(y))

     q=(Y0+beta)/(N+1)

     lnq=log(q)

     lnq[is.na(Y)]=0

     return(lnq)

  }#lnqFunc 

  genefreq=foreach(vial=atts$Pop_Vial,.combine=rbind)%dopar%{

    ats=atts[vial,]

    xx=x[vial,,]

    yy=y

    yy[ats$FromPop,,]=y[ats$FromPop,,]-xx

    lnq=lnqFunc(yy)

    gg=p*exp(Reduce("+",lapply(seq(max(nalleles)),function(allele){xx[,allele]%*%t(lnq[,,allele])})))[1,]

    gg/sum(gg)
      
  }#vial

  stopCluster(cl)

  dimnames(genefreq)=list(atts$Pop_Vial,popvec)

  Rgenefreq=t(rowsum(t(genefreq),group=groupvec))

  colnames(Rgenefreq)=groupnames

  time=round(Sys.time()-begin.time,2)

  print(time)

  return(list(Attributes=atts,IndividualByGroup=Rgenefreq,IndividualByPop=genefreq,RunTime=time))
  
}
