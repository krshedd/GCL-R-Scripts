IndividualAssignmentSummary.GCL=function(GroupNames,groupvec,mixnames,BAYESoutputDir,nchains,nreps,burn=1/2,thin=100){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function will summarize BAYES individual assignment .CLS files.
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # GroupNames = reporting group names
  #   ~ A character vector indicating group names.
  # groupvec = reporting group affiliation for each baseline populations
  #   ~ A numeric vector indicating reporting group affiliation for each
  #     baseline population. Numbers correspond  to GroupNames.
  # mixnames = mixtures names
  #   ~ A character vector indicating mixtures names (i.e. name of folder where
  #     BAYES output files exist.
  # BAYESoutputDir = filepath to where BAYES output folders exist.
  #   ~ A character vector of length one indicating where BAYES output folders
  #     exist.
  # nchians = number of BAYES chains
  #   ~ A numeric vector indicating number of BAYES chains. Can be a vector of
  #     length(mixnames) if there are a different number of chains per mixture.
  # nreps = number of MCMC interations
  #   ~ A numeric vector indicating number of MCMC reps.  Can be a vector of
  #     length(mixnames) if there are a different number of reps per mixture.
  # burn = proportion of MCMC iterations for burn-in
  #   ~ Typically 0.5
  # thin = factor BAYES MCMC reps were thinned by for Individual Assignment
  #   ~ Number from the 3rd argument for thin in CreateControlFile.GCL.
  #     Recommend thin=c(1,1,100) for IA, so this argument is set to 100.
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # List of length mixnames
  #   Matrix with nrow = number of individuals in mixture and ncol = number of
  #   reporting groups. Numbers indicate proportion of thinned MCMC iterations
  #   in which an individual was assigned to a specific reporting group.
  #   Individual names come from FK_FISH_ID if the mixture .gcl object exists,
  #   otherwise they are sequential.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # GroupNames=c("KarlukEarly","KarlukLate")
  # groupvec=c(1,1,2,2,1,1,1,2,1,2,1,2,2,1,1,2)
  # mixnames=c("KarlukEarly_0.1_1")
  # BAYESoutputDir="V:/Analysis/4_Westward/Sockeye/Karluk Smolt 2013-2015/Baseline/BAYES/BAYES.output/IndividualAssignment/Test"
  # nchains=5
  # nreps=40000
  # burn=0.5
  # thin=100
  #
  #  Written By JJ 6/28/13
  #  Modified by Kyle Shedd 2/25/15
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  names(mixnames)=mixnames
  
  if(length(nchains)==1){
    
    nchains=rep(nchains,length(mixnames))
    
    names(nchains)=mixnames
    
  }
  
  if(length(nreps)==1){
    
    nreps=rep(nreps,length(mixnames))
    
    names(nreps)=mixnames
    
  }
  
  myoutputdirs=sapply(mixnames,function(mixname){paste(BAYESoutputDir,mixname,sep="/")},simplify=FALSE)
  
  myfilenames=sapply(mixnames,function(mixname){paste(myoutputdirs[mixname],"/",mixnames[mixname],"Chain",seq(nchains[mixname]),"CLS.CLS",sep="")},simplify=FALSE)
  
  output=results=vector("list",length(mixnames))
  
  names(output)=names(results)=mixnames
  
  for(mixname in mixnames){
    
    testoutput=read.table(myfilenames[[mixname]][1],header=FALSE)
    
    mixsampsize=dim(testoutput)[1]/floor(nreps[mixname]/thin)
    
    chains=paste("Chain",seq(nchains[mixname]),sep="")
    
    skip=mixsampsize*floor(nreps[mixname]*burn/thin)
    
    output[[mixname]]=NULL
    
    if(exists(paste(mixname, ".gcl", sep = ''))) {
      
      my.gcl <- get(paste(mixname, ".gcl", sep = ''))
      
      indnames <- my.gcl$attributes$FK_FISH_ID
      
    } else {
      
      indnames <- seq(mixsampsize)
      
    }
    
    results[[mixname]]=array(NA,c(mixsampsize,max(groupvec)),list(indnames,GroupNames))
    
    for(chain in chains){
      
      output[[mixname]]=rbind(output[[mixname]],read.table(myfilenames[[mixname]][match(chain,chains)],header=FALSE,skip=skip))
      
    }#chain 
    nits=nrow(output[[mixname]])/mixsampsize
    
    C=ncol(output[[mixname]])
    
    z=apply(output[[mixname]],1,which.max)
    
    zsummry=t(sapply(tapply(z,rep(seq(mixsampsize),nits),function(zz){tabulate(zz,C)/nits}),function(zzz)zzz))
    
    for(g in 1:max(groupvec)){
      
      results[[mixname]][,g]=apply(cbind(zsummry[,groupvec%in%g]),1,sum)
      
    }#g
    
  }#mixname 
  
  return(results)
  
}
