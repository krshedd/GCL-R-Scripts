RemoveNAindv.GCL <- function(sillyvec) {
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  This function removes any individuals from a SILLY that were not genotyped 
#  for all markers in markersuite.
#
#  LOKI2R.GCL pulls all individuals from a SILLY genotyped for any loci 
#  existing in LocusControl$locusnames. Sometimes we have individuals that were
#  only genotyped for a subset of loci in our markersuite. This function
#  removes those individuals from SILLYs
#
#  Argument(s):  
#  sillyvec <- character vector of SILLYs
#
#  Output:
#  1) this function will modify silly.gcls in your environment to remove NA fish
#  2) output is a list of FK_FISH_ID of fish removed for each SILLY
#
#  Written by Kyle Shedd 10/16/15  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  na.individuals.to.remove <- sapply(sillyvec, function(silly) {
    my.gcl <- get(paste(silly, ".gcl", sep = ""))
    my.gcl$attributes$FK_FISH_ID[apply(my.gcl$scores, 1, function(ind) {
      sum(is.na(ind[, 1])) > 0 
    })]
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  sapply(sillyvec, function(silly) {
    if(length(na.individuals.to.remove[[silly]]) > 0){
      my.gcl <- get(paste(silly, ".gcl", sep = ""))
      counts <- my.gcl$counts
      scores <- my.gcl$scores
      n <- my.gcl$n
      attributes <- my.gcl$attributes
      
      ind.to.keep <- !my.gcl$attributes$FK_FISH_ID %in% na.individuals.to.remove[[silly]]
      
      counts <- counts[ind.to.keep, , ,drop=FALSE]
      scores <- scores[ind.to.keep, , ,drop=FALSE]
      n <- nrow(scores)
      attributes <- attributes[ind.to.keep, ,drop=FALSE]
      
      assign(paste(silly, ".gcl", sep = ""), list(counts = counts, scores = scores, n = n, attributes = attributes), pos = 1)
    }} )
  
  return(na.individuals.to.remove)
}