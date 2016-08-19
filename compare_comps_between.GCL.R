compare_comps_between.GCL <- function(mixnames, groupnames, mixdir, diffs = seq(10)/100, nchains = 5, burn = 1/2){

##########  Arguments  #####################################################################################################################################################################################################################################################################
#
#  mixnames -- characacter vector of length 2 containing the mixture names to compare.
#  
#  groupnames -- characacter vector giving the names of the stock groups resolved by the mixtures, with length equal to the number columns in the RGN files.
#  
#  mixdir -- atomic string giving the folder path where the mixture sub-folders are located.
#
#
##########  Return Value  ###################################################################################################################################################################################################################################################################
# 
#   retuns a vector of p-values, same length as 'groupnames', with each element corresponding to the proportion of 
#   MCMC relizations where a group's composition in one mixture exceeds that of the other mixture.
#  
##########  Example  ###################################################################################################################################################################################################################################################################
#  
#  load("V:/Analysis/2_Central/Sockeye/Cook Inlet/2012 Baseline/Mixture/2013 UCIfisheryMixtures/2013UCIfisheryMixtureAnalysis.RData")
#
#  mixnames <- c("DriftExpCorr.Jul11","Drift.Jul8")
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

  dsim0 <- Reduce("-", output)

  colnames(dsim0) <- groupnames

  layout(mat = matrix(seq.int(length(groupnames) + length(groupnames) %% 2), byrow = TRUE, ncol = 2))

  xy <- lapply(lapply(lapply(dsim0, hist, plot = FALSE), FUN = "[", value = c("density", "mids")), FUN =  setNames, nm = c("y", "x"))

  xlim <- range(unlist(lapply(xy, "[[", "x")))

  ylim <- range(unlist(lapply(xy, "[[", "y")))

  par(mar = c(2, 0, 2, 0) + 0.01)

  plots <- lapply(groupnames, function(nm){plot(xy[[nm]], type = "l", yaxt = "n", xlim = xlim, ylim = ylim, lwd = 2) ; abline(v = 0, lty = "dotted", lwd = 2) ; text(xlim[2], ylim[2], labels = nm, pos = 2)})

  dsim <- t(diag((-1) ^ (apply(dsim0, 2, mean) < 0)) %*% t(dsim0))

  one.sided <- cbind(one.sided = c(setNames(apply(dsim < 0, 2, mean), groupnames), Overall = mean(apply(dsim < 0, 1, all))), diff.mean = c(apply(dsim0, 2, mean), NA))

  two.sided <- sapply(setNames(diffs, diffs), function(dd){c(apply(abs(dsim0) < dd, 2, mean), Overall = mean(apply(abs(dsim0) < dd, 1, all)))})  

  return(list(one.sided = one.sided, two.sided = two.sided))

}  



