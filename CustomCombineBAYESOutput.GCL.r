CustomCombineBAYESOutput.GCL=function(groupvec, groupnames, maindir, mixvec, prior="", ext="RGN", nchains=3, burn=0.5, alpha=0.1,PosteriorOutput=FALSE,threshhold=5e-7){

# groupvec <- c(1,1,1,1,1,1,1,1,2,3,5,5,4,5) 
# groupnames <- c("Alaska","Nass","Skeena","Fraser","Other") 
# maindir <- "V:\\WORK\\SARA\\Sockeye\\SEAK\\Analysis\\Mixtures\\0405Dists101104\\Bayes\\Feb2009\\Output\\2005dist104" 
# mixvec <- c("year2005dist104stat2829a","year2005dist104stat30a","year2005dist104stat31a","year2005dist104stat32a","year2005dist104stat33a","year2005dist104stat34a") 
# prior <- "Flat"

  require(coda)

  G <- max(groupvec)

  C <- length(groupvec)

  nummix <- length(mixvec)

  results <- setNames(vector("list",nummix),mixvec)

  Output <- setNames(vector("list",nummix),mixvec)
  
  for(mix in mixvec){

    filenames <- paste(maindir,"\\",mix,"\\",mix,prior,"Chain",1:nchains,ext,".",ext,sep="")

    files <- lapply(filenames,function(filename){mcmc(as.matrix(read.table(filename)[,-1]))})

    end <- sapply(files,nrow)

    if(length(unique(end))>1){stop("Chains must be the same length!!!")}

    end <- end[1]

    begin <- floor(burn*end)+1

    files4GR <- vector("list",nchains)

    for(chain in seq(nchains)){

      files4GR[[chain]] <- as.mcmc(t(rowsum(t(files[[chain]][begin:end,]),group=groupvec)))

    }#chain

    files4GR <- as.mcmc.list(files4GR)    

    if(nchains<2){

      GR <- list(cbind(rep(NA,G)))

    }else{    

      GR <- gelman.diag(files4GR,multivariate=FALSE,transform=TRUE) 

    }
    
    Output[[mix]] <- Reduce(rbind,lapply(files4GR,as.matrix))

    results[[mix]] <- array(NA,c(G,7),dimnames=list(groupnames,c("mean","sd","median",paste(round(alpha/2,3)*100,"%",sep=""),paste(round(1-alpha/2,3)*100,"%",sep=""),"P=0","GR")))

    results[[mix]][groupnames,1] <- apply(Output[[mix]],2,mean)

    results[[mix]][groupnames,2] <- apply(Output[[mix]],2,sd)

    results[[mix]][groupnames,3] <- apply(Output[[mix]],2,median)

    results[[mix]][groupnames,4] <- apply(Output[[mix]],2,quantile,probs=alpha/2)

    results[[mix]][groupnames,5] <- apply(Output[[mix]],2,quantile,probs=1-alpha/2)

    results[[mix]][groupnames,6] <- apply(Output[[mix]],2,function(clm){sum(clm<threshhold)/length(clm)})

    results[[mix]][groupnames,7] <- GR[[1]][,1]    

  }

  if(PosteriorOutput){

     ans <- list(Stats=results,Output=Output)
 
  }else{

    ans <- results

  }

  return(ans)

}

