treeColor.GCL=function(tree,currentnames,treenames,groupvec,regioncol,regionpch=NULL){
##############################################################################################################################################################
#
#  "currentnames", "treenames", and "groupvec" all have the same order.
#
#  "currentnames" are the names on the tree (but in different order) that you want changed to "treenames".  
#
#  "groupvec" is a vector group membership.
#
#  "regioncol" is a vector the same length as max(groupvec) (i.e. the number of groups), each element is an R color number. 
#
#  "regionpch" has the same structure as "regioncol", each element is an R pch number. 
#
# Written By Jim Jasper 6/22/2011
##############################################################################################################################################################

  xx=tree
  pops1=currentnames
  pops2=treenames
  pops3=groupvec
  C=length(pops1)

  xx$tip.label=pops2[match(xx$tip.label,pops1)]

  regInd=pops3[match(xx$tip.label,pops2)]
  cols=regioncol[regInd]
  edgeind=match(xx$edge[,2],1:C)
  EdgeReg=regInd[edgeind]


  all=cbind(1:length(EdgeReg),xx$edge[,1],EdgeReg)

  edgereg=EdgeReg
  MM=length(EdgeReg)
  kk=1
  while(kk){
  NN=dim(all)[1]
  if(NN==0) break()
    ind=(apply(all[-NN,2:3]==all[-1,2:3],1,prod)==1)
    ind[is.na(ind)]=F
    ind[1]=F
    ind2=match(all[,1][ind],all[,1])
    all[ind2-1,3]=all[ind2,3]
    edgereg[all[,1][ind]-1]=all[ind2,3]
    all=all[-c(ind2,ind2+1),] 
  }

  EdgeReg=edgereg

  edgecol=regioncol[EdgeReg]

  edgecol[is.na(edgecol)]="black"

  if(!is.null(regionpch)){
    mypch=regionpch[regInd]
  }
  else{mypch=NULL}

  return(list(tree=xx,color=edgecol,pch=mypch))

}











