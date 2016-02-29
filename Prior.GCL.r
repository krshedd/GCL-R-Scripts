Prior.GCL=function(groupvec,groupweights,minval=0.01){

  Cg=table(groupvec)

  IND=groupweights<minval

  groupweights[!IND]=(1-minval*sum(IND))*groupweights[!IND]/sum(groupweights[!IND])  

  groupweights[IND]=minval

  popweights=groupweights[groupvec]/Cg[groupvec]

  return(popweights)
}