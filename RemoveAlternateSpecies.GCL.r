RemoveAlternateSpecies.GCL <- function(AlternateSpeciesReport, AlternateCutOff = 0.5, FailedCutOff = 0.5, NonmissingCutOff){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function removes fish that have been indentified, via genetic markers, as the wrong species.
  #   It relies on output from FindAlternateSpecies.GCL. 
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   AlternateSpeciesReport - The object created by FindAlternateSpecies.GCL
  #   AlternateCutOff - The minimum proportion of alternate alleles an individual must have in order to be considered the 'wrong species'; default 0.5
  #   FailedCutOff - The minimum proportion of failed loci an individual must have in order to be considered the 'wrong species'; default 0.5
  #   NonmissingCutOff - The minimum number of nonmissing alternate loci an individual must have in order to be considered the 'wrong species'
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Tibble showing which fish were removed and their stats (alternate, failed markers)
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   source("/R/Functions.GCL.R")
  #   source(file = "Examples/QCExample.R")
  #   wrong_spp <- FindAlternateSpecies.GCL(sillyvec = sillys, species = "sockeye")
  #   rm_spp <- RemoveAlternateSpecies.GCL(AlternateSpeciesReport = wrong_spp , AlternateCutOff = 0.2, FailedCutOff = 0.2, NonmissingCutOff = 3)
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function requires an object from FindAlternateSpecies.GCL.
  #   
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse)  # Install packages, if not in library and then load them.
  
  RemovedSpp <- AlternateSpeciesReport %>% # This is output from FindAlternateSpecies.GCL
    tidyr::separate(col = "silly_fish", into = c("silly", "ID"), sep = "\\_(?=[^\\_]+$)") %>% # Split out silly and fish
    dplyr::mutate(
      fate =
        dplyr::case_when(
          alternate >= AlternateCutOff & failure >= FailedCutOff & non_missing_alt >= NonmissingCutOff ~ "remove",
          TRUE ~ "keep"
        ) # If alternates AND failures AND # of non_missing_alts are greater than cutoffs, mark as remove
    ) %>%
    dplyr::filter(fate == "remove") # Just get fish marked to remove
  
  output <- lapply(RemovedSpp$silly %>% unique(), function(silly){
    
    IDsToRemove <-  RemovedSpp %>% 
      dplyr::filter(silly == !!silly) %>% 
      dplyr::pull(ID)
    
    if(!purrr::is_empty(IDsToRemove)) {
      
      RemoveIDs.GCL(silly = silly, IDs = IDsToRemove)
      
    }else(data.frame())

  }) %>%
    dplyr::bind_rows()
  
  if(is_empty(output)){
    
    warning("No fish were identified as the wrong species.")
    
  }
  
  return(output)
  
}

