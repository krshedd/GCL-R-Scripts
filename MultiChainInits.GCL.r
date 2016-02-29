MultiChainInits.GCL<-function(npops,nchains,prop=.9){
############################################################################
#This function creates a matrix of initial start values, where the values for each 
#chain add to 1 and a higher start values are given to a different range of pops for each chain.
#For example, if you have 500 pops and need starting values for 5 chains and "prop"
#is set at .9, the first 100 pops will have starting values = .9/100 and the rest
#of the pops will have start values = .1/400.  Each chain will have a different 
#set of pops with this higher start value.  
############################################################################
#npops - number of populations or length(silyvec)
#nchains - number of BAYES chains that you need starting values for.
#prop - start value proportion you want to give npops/nchains, default is .9
############################################################################
#Written by Tyler Dann and Andy Barclay 11/22/2011
############################################################################

  if(nchains<2){stop("Need more than one chain, hoser!!!")}

  r=npops%%nchains

  avgpop=floor(npops/nchains)

  groupvec=c(rep(1:(nchains),each=avgpop),rep(nchains,r)) 

  GroupWeights=array((1-prop)/(nchains-1),c(nchains,nchains))

  diag(GroupWeights)=prop  

  initmat=apply(GroupWeights,1,function(groupweights){
            Prior.GCL(groupvec=groupvec,groupweights=groupweights,minval=0)
          })  

  
  dimnames(initmat)=list(1:npops,paste("Chain",1:nchains,sep=''))

  return(initmat)
}


