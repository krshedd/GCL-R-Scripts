RandomInits.GCL=function(groupvec,groupweights=rep(1/max(groupvec),max(groupvec)),nchains=5){
########################################################################################################################
#  This function generates a matrix of randomly generated initial stock proportions.
#
#
#
# Written by Jim Jasper 6/16/2011
########################################################################################################################

  chains=paste("Chain",1:nchains,sep="")

  Cg=table(groupvec)
  
  popweights=groupweights[groupvec]/Cg[groupvec]

  ans=sapply(chains,function(chain){sapply(popweights,function(popweight){rgamma(1,popweight,1)})})

  ans=apply(ans,2,function(x){x/sum(x)})

  return(ans)
}

