ConfusionMatrices.GCL=function(LeaveOneOutDist,groupnames,groupvec,sillyvec){

  #LeaveOneOutDist=list(RT,T)

  myORD=order(groupvec)

  groupvec=groupvec[myORD]

  sillyvec=sillyvec[myORD] 

  C=length(groupvec)

  G=max(groupvec)

  Cg=as.vector(table(groupvec))

  cmCg=c(0,cumsum(Cg)) 

  ConfusionMatrixPopByPop=t(matrix(rapply(LeaveOneOutDist[[2]],mean),nrow=C,byrow=TRUE,dimnames=list(paste(groupnames[groupvec],sillyvec,sep="."),paste(groupnames[groupvec],sillyvec,sep="."))))

  ConfusionMatrixGroupByPop=t(matrix(rapply(LeaveOneOutDist[[1]],mean),nrow=G,byrow=TRUE,dimnames=list(groupnames,paste(groupnames[groupvec],sillyvec,sep="."))))

  ConfusionMatrixGroupByGroup=apply(ConfusionMatrixGroupByPop,2,function(clm){sapply(seq(G),function(g){mean(clm[groupvec==g]) }) })
  dimnames(ConfusionMatrixGroupByGroup)[[1]]=groupnames
                               

  return(list(GroupByGroup=ConfusionMatrixGroupByGroup,GroupByPop=ConfusionMatrixGroupByPop,PopByPop=ConfusionMatrixPopByPop))

}