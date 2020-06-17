RemoveDups.GCL <- function(dupcheck){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function removes duplicated "IDs" with the highest number of missing loci. 
  #  If both IDs have the same number of missing loci, ID1 is removed.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   dupcheck - an object created by the function "CheckDupWithinSilly.GCL". 
  #
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #  This function returns a tibble of IDs removed for each silly in dupcheck.
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  password = "************"
  #  CreateLocusControl.GCL(markersuite = "Sockeye2011_96SNPs", username = "awbarclay", password = password)
  #  sillyvec = c("SMCDO03", "SNEVA13")
  #  LOKI2R.GCL(sillyvec = sillyvec, username = "awbarclay", password = password)
  #  RemoveIndMissLoci.GCL(sillyvec = sillyvec)
  #
  #  dupcheck <- CheckDupWithinSilly.GCL(sillyvec = sillyvec, loci = LocusControl$locusnames, quantile = NULL, minproportion = 0.95, ncores = 8)
  #  removed_dups <- RemoveDups.GCL(dupcheck)
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   When quantile is set to NULL this function utilizes rubias::close_matching_samples() to perform the duplicate check and it much faster than when you set a quantile.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(is_empty(dupcheck)){stop("Nothing removed. There are no duplicates to remove in dupcheck.", call. = FALSE)}
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse)  # Install packages, if not in library and then load them.

  if(!tibble::is_tibble(dupcheck)){dupcheck <- dupcheck$report}
  
  dupcheck_names <-  c("silly", "ID1", "ID2", "Missing1", "Missing2", "proportion")
  
  if(length(setdiff(dupcheck_names, names(dupcheck))) > 0){
    
    stop(paste0("Nothing removed. Dupcheck must contain the folling variables: ", paste0(dupcheck_names, collapse = ", ")))
    
    }
  
  sillys <- dupcheck$silly %>% 
    unique()
  
  to_remove <- dupcheck %>%
    mutate(remove = case_when(Missing1 > Missing2~ID1,
                              Missing2 > Missing1~ID2,
                              Missing1 == Missing2~ID1))
  
  results <- lapply(sillys, function(silly){
    
    remove <- to_remove %>% 
      filter(silly==!!silly) %>% 
      select(silly, removed_IDs = remove)
    
    RemoveIDs.GCL(silly = silly, IDs = remove$removed_IDs)
    
    remove
    
  }) %>% bind_rows()

  return(results)
}