ConfusionMatrices.GCL=function(LeaveOneOutDist,groupnames,groupvec,sillyvec){
  ####################################################################################################################################################################################################################################################################
  #
  # This function takes the genetic likelihood profile from LeaveOneOutDist.GCL
  # and summarizes into matrices at the group, group by pop, and pop by pop level.
  # This genetic likelihood profile is intended to show how well populations / 
  # reporting groups are differentiated by the markers at hand. These probabilities
  # can give an idea of potential misallocation and where.
  # This function can aid in marker selection.
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # LeaveOneOutDist = 
  #   ~ Unmodified list output from LeaveOneOutDist.GCL.
  # groupnames = Groups15
  #   ~ A character vector groupnames.
  #   ~ This provided dimnames for matrices
  # groupvec = Kodiak49GroupVec9
  #   ~ A numberic vector of length(sillyvec) sepcifying group assignments.
  # sillyvec = Kodiak49Pops
  #   ~ A character vector of sillys (pops make more sense, but collecitons could be used).
  #   ~ This provided dimnames for matrices
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Named List of 3
  #  GroupByGroup matrix (number of RGs x number of RGs)
  #   average probability of each individual from a group originating from a pop in that group
  #  GroupByPop matrix (number of pops x number of RGs)
  #   average probability of each individual from a pop originating from a pop in that group
  #  PopByPop matrix (number of RGs x number of RGs)
  #   average probability of each individual from a pop originating from a pop
  #
  # Created by Jim Jasper on unknown???
  # Modified by Kyle Shedd on unkown to include description
  #
  ####################################################################################################################################################################################################################################################################
  
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