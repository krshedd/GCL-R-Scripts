RemoveAlternateSpecies.GCL=function(AlternateSpeciesReport, AlternateCutOff = 0.5, FailedCutOff = 0.5, NonmissingCutOff){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function removes fish that have been indentified, via genetic markers, as the wrong species.
  #   It relies on output from FindAlternateSpecies.GCL. 
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   AlternateSpeciesReport - The object created by FindAlternateSpecies.GCL
  #   AlternateCutOff - The percent similarity that you want to call a fish as 'wrong species'; default 0.5
  #   FailedCutOff - The percent similarity that you want to call a fish as 'wrong species'; default 0.5
  #   NonmissingCutOff - The number of nonmissing alternate loci
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
    tidyr::separate(col = "silly_fish", into = c("silly", "ID"), sep = "_") %>% # Split out silly and fish
    dplyr::mutate(
      fate =
        dplyr::case_when(
          alternate >= AlternateCutOff & failure >= FailedCutOff & non_missing_alt >= NonmissingCutOff ~ "remove",
          TRUE ~ "keep"
        ) # If alternates AND failures AND # of non_missing_alts are greater than cutoffs, mark as remove
    ) %>%
    dplyr::filter(fate == "remove") # Just get fish marked to remove
  
  #RemoveIDs.GCL(silly = RemovedSpp$silly, IDs = RemovedSpp$ID) # Call on RemoveIDs to take care of removal.
  
  output <- lapply(
    RemovedSpp$silly %>% unique(), function(silly){
    
    IDsToRemove <-  RemovedSpp %>% 
      dplyr::filter(silly == !!silly) %>% 
      dplyr::pull(ID)
    
    if(!purrr::is_empty(IDsToRemove)) {
      
      RemoveIDs.GCL(silly = silly, IDs = IDsToRemove)
      
    }

  }) %>%
    dplyr::bind_rows()
  
  return(output)
}

