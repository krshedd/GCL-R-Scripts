LinkageCorrelationJAGS.GCL=function(sillyvec,markerset,groupvec,groupnames){
#####################################################################################################################################################################################################################################################################
#
#
#
# Written by Jim Jasper 2/17/2014
#####################################################################################################################################################################################################################################################################

  if(sum(is.na(match(markerset,LocusControl$locusnames)))){stop(paste("'",markerset[is.na(match(markerset,LocusControl$locusnames))],"' from argument 'markerset' not found in 'LocusControl' object!!!",sep=""))}

  model="V:\\DATA\\R_GEN\\BUGS Scripts\\LinkageCorrelationCoefficientJAGS.txt"

  ORD=order(groupvec)

  Cg=as.vector(table(groupvec[ORD]))

  cmCg=c(0,cumsum(Cg))

  G=length(Cg)

  C=length(sillyvec) 

  locus=paste(markerset,collapse=".")  

  tempsillyvec=paste(sillyvec[ORD],".temp",sep="")

  IND=is.na(match(locus,LocusControl$locusnames))

  if(IND){

    CombineLoci.GCL(sillyvec,markerset,update=FALSE)

    counts=sapply(tempsillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)$counts[,locus,1:9]},simplify=FALSE)
  
    y=t(sapply(tempsillyvec,function(silly){apply(counts[[silly]],2,function(col){sum(col[!is.na(col)])})}))

    suppressWarnings(rm(list=paste(tempsillyvec,".gcl",sep="")))

  }
  if(!IND){

    counts=sapply(sillyvec[ORD],function(silly){get(paste(silly,".gcl",sep=""),pos=1)$counts[,locus,1:9]},simplify=FALSE)

    y=t(sapply(sillyvec[ORD],function(silly){apply(counts[[silly]],2,function(col){sum(col[!is.na(col)])})}))

  }
 
  require(R2jags)

  n=apply(y,1,sum)

  data=list(y=y,n=n,C=C,G=G,cmCg=cmCg)

  pars2save=c("rGroup","p1Gst","p2Gst","fGst","GstRatio")#  

  fit=jags(data=data,parameters.to.save=pars2save,model.file=model,n.chains=1,n.iter=10000,n.thin=1,DIC=FALSE)$BUGSoutput$summary
  
  detach("package:R2jags")

  return(fit)

}

