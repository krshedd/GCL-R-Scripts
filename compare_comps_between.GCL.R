compare_comps_between.GCL <- function(mixnames, groupnames, mixdir, d = seq(10)/100, nchains = 5, burn = 1/2){

##########  Arguments  #####################################################################################################################################################################################################################################################################
#
#  mixnames -- characacter vector of length 2 containing the mixture names to compare.
#  
#  groupnames -- characacter vector giving the names of the stock groups resolved by the mixtures, with length equal to the number columns in the RGN files.
#  
#  mixdir -- atomic string giving the folder path where the mixture sub-folders are located.
#
#  d -- atomic numeric giving a biologically significant difference between stock compositions. Ignored if 'onesided' is set to TRUE. 
#
#  onesided -- should a one-sided test be cunducted?
#
##########  Return Value  ###################################################################################################################################################################################################################################################################
# 
#   if 'onesided' is set to TRUE, retuns a vector of p-values, same length as 'groupnames', with each element corresponding to the proportion of 
#   MCMC relizations where a group's composition in one mixture exceeds that of the other mixture.
#  
#   if 'onesided' is set to FALSE, 'd' must be specified, and a numeric vector of p-values is returned the same length as 1 + length(groupnames). 
#   Each element corresponds the proportion of MCMC relizations where a group's composition in one mixture exceeds a distance 'd' of the other mixture. 
#   The last element provides the overall proportion of realizations where all groups simultaneously differ.   
#
##########  Example  ###################################################################################################################################################################################################################################################################
#  
#  load("V:/Analysis/2_Central/Sockeye/Cook Inlet/2012 Baseline/Mixture/2013 UCIfisheryMixtures/2013UCIfisheryMixtureAnalysis.RData")
#
#  mixnames <- c("DriftExpCorr.Jul11","Drift.Jul8")
#
#  d <- seq(10)/100 
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

  dirs <- paste(mixdir, "/", mixnames, "/", mixnames, sep="")

  chains=paste("Chain",seq(nchains), "RGN.RGN",sep="")

  fnames <- matrix(apply(expand.grid(dirs, chains), 1, paste, collapse=""), nrow=nchains, byrow=TRUE)

  colnames(fnames) <- mixnames

  files <- sapply(mixnames, function(mix){lapply(fnames[,mix], read.table)}, simplify = FALSE)

  output <- sapply(mixnames, function(mix){Reduce(rbind, lapply(files[[mix]], function(file){file[seq(floor(nrow(file) * burn) + 1, nrow(file)), -1]}))}, simplify = FALSE)

  dsim <-Reduce("-", output)

  colnames(dsim) <- groupnames

  pvals1side <- c(apply(dsim < 0, 2, mean), Overall = NA)

  pvals1side <- pmin(pvals1side, 1 - pvals1side)

  pvals <- sapply(setNames(d,d), function(dd){c(apply(abs(dsim) < dd,2,mean), Overall = mean(apply(abs(dsim) < dd,1,all)))})

  return(list(pvals = cbind(one.sided = pvals1side, diff.mean = c(apply(dsim,2,mean), NA), pvals), differnce.output = dsim))

}  



