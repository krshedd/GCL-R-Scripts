FailureRate.GCL <- function(sillyvec) {
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  This function calculates failure rates by SILLY, loci, and project
#
#  It is IMPORTANT to run this function after "RemoveNAindv.GCL" but before any
#  other steps in the QC (i.e. removing duplicate fish, alternate fish, 
#  missing fish, etc.). This function does NOT connect to LOKI. It merely
#  calculates failure rates (0's / total fish run) from the silly.gcl objects.
#
#  Argument(s):  
#  sillyvec <- character vector of SILLYs in the project
#
#  Output:
#  List of failure rate by SILLY, locus, and project
#
#  Written by Kyle Shedd 10/16/15  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  numran <- t(sapply(sillyvec, function(silly) {
    silly.gcl <- get(paste(silly, ".gcl", sep = ""))
    apply(silly.gcl$scores, 2, function(loci) {sum(!is.na(loci[, "Dose1"]))} )
  }))
  
  require(Matrix)
  numscored <- t(sapply(sillyvec, function(silly) {
    silly.gcl <- get(paste(silly, ".gcl", sep = ""))
    apply(silly.gcl$scores, 2, function(loci) {nnzero(loci[, "Dose1"])} )
  }))
  
  
  if(length(sillyvec) == 1){
    silly.failure <- 1 - (apply(numscored, 1, sum) / sum(numran))
    locus.failure <- 1 - (apply(numscored, 2, sum) / as.vector(numran))
    overall.failure <- silly.failure
  }
  
  if(length(sillyvec) > 1){
    silly.failure <- 1 - (apply(numscored, 1, sum) / apply(numran, 1, sum))
    locus.failure <- 1 - (apply(numscored, 2, sum) / apply(numran, 2, sum))
    overall.failure <- sum(apply(numran, 1, sum) * silly.failure) / sum(apply(numran, 1, sum))
  }
  
  output <- list(Silly_Failure_Rate = silly.failure, Locus_Failure_Rate = locus.failure, Overall_Failure_Rate = overall.failure)

  return(output)
}