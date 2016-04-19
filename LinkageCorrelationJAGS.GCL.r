LinkageCorrelationJAGS.GCL=function(sillyvec,markerset,groupvec,groupnames){
  #####################################################################################################################################################################################################################################################################
  #
  # This function calculates Gst values for the all combinations of the markerset
  # (i.e. each locus separately and also combined) for your sillyvec according to 
  # the groupvec that is given. Additioanlly, it provides regional and population
  # level correlation coefficients (r) between the markers. 
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # sillyvec = KMA473Pops
  #   ~ A character vector of sillys (pops make more sense, but collecitons could be used).
  # markerset = c("One_MHC2_190", "One_MHC2_251")
  #   ~ A character vector of two loci.
  # groupvec = KMA473PopsGroupVec15
  #   ~ A numberic vector of length(sillyvec) sepcifying group assignments.
  # groupnames = Groups15Short
  #   ~ A character vector of length unique(sillyvec) specifying group names.
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # List of 3
  #  Gst
  #   numeric matrix of Gst values
  #    number of rows = 4 (Gst Ratio = combined loci Gst / highest single locus Gst;
  #    Combined loci; Locus 1; Locus 2)
  #    number of columns = 7 (mean, sd, 2.5%, 25%, 50%, 75%, and 97.5% quartiles)
  #  Region
  #   numeric matrix of reporting group correlation coefficients
  #    number of rows = length(groupnames)
  #    number of columns = 7 (mean, sd, 2.5%, 25%, 50%, 75%, and 97.5% quartiles)
  #  Pop
  #   numeric matrix of population correlation coefficients
  #    number of rows = length(sillvec)
  #    number of columns = 7 (mean, sd, 2.5%, 25%, 50%, 75%, and 97.5% quartiles)
  #
  # Written by Jim Jasper 2/17/2014
  # Edited by Kyle Shedd 7/29/2015 to get Pop correlation coefficients
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
  
  pars2save=c("rGroup","p1Gst","p2Gst","fGst","GstRatio","r")#  
  
  fit=jags(data=data,parameters.to.save=pars2save,model.file=model,n.chains=1,n.iter=10000,n.thin=1,DIC=FALSE)$BUGSoutput$summary
  
  detach("package:R2jags")
  
  Gst = fit[c(1:4), ]
  rownames(Gst) = c("GstRatio", paste(paste(markerset, collapse = "."), "Gst", sep = "."), paste(markerset[1], "Gst", sep = "."), paste(markerset[2], "Gst", sep = "."))
  
  Region = fit[-c(1:(4 + C)), ]
  rownames(Region) = groupnames
  
  Pop = fit[5:(4+C), ]
  rownames(Pop) = sillyvec
  
  fit = list(Gst = Gst, Region = Region, Pop = Pop)
  
  return(fit)
  
}