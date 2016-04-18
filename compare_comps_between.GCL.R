compare_comps_between.GCL <- function(mixnames, groupnames, mixdir, d = 0.05, onesided = FALSE, nchains = 5, burn = 1/2){

########################################################################################################################################################################################################################################################################################
#
#  load("V:/WORK/Sockeye/Cook Inlet/2012 Baseline/Mixture/2013 UCIfisheryMixtures/2013UCIfisheryMixtureAnalysis.RData")
#
#  mixnames <- c("DriftExpCorr.Jul11","Drift.Jul8")
#
#  d <- 0.05    
#
#  groupnames <- groups  
#
#  mixdir <- "V:/Analysis/2_Central/Sockeye/Cook Inlet/2012 Baseline/Mixture/2013 UCIfisheryMixtures/BAYES/Output"
#
#  nchains <- 5
#
#  burn <- 1/2
#
########################################################################################################################################################################################################################################################################################

  dirs <- paste(dir, "/", mixnames, "/", mixnames, sep="")

  chains=paste("Chain",seq(nchains), "RGN.RGN",sep="")

  fnames <- matrix(apply(expand.grid(dirs, chains), 1, paste, collapse=""), nrow=nchains, byrow=TRUE)

  colnames(fnames) <- mixnames

  files <- sapply(mixnames, function(mix){lapply(fnames[,mix], read.table)}, simplify = FALSE)

  output <- sapply(mixnames, function(mix){Reduce(rbind, lapply(files[[mix]], function(file){file[seq(floor(nrow(file) * burn) + 1, nrow(file)), -1]}))}, simplify = FALSE)

  dsim <-Reduce("-", output)

  colnames(dsim) <- groupnames

  if(onesided){

    dsimBOOL  <- dsim < 0

    pvals <- c(apply(dsimBOOL, 2, mean))

    pvals <- pmin(pvals, 1 - pvals)

  }
   

  if(!onesided){

    absd <- abs(dsim)

    absdBOOL  <- absd < d

    pvals <- c(apply(absdBOOL,2,mean), Overall = mean(apply(absdBOOL,1,all)))

  }

  

  return(pvals)

}  



