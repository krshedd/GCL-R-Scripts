RemoveIDs.GCL <- function(silly, IDs){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function removes "IDs" (FK_FISH_ID) from the "*.gcl" object associated with "silly".
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   silly - a silly code without the ".gcl" extention.
  #
  #   IDs - a numeric or character vector of FK_FISH_ID's to remove from "silly".
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  #  Assigns "*.gcl" objects with removed IDs to the current workspace. 
  #
  #  The function returns a tibble with 3 variables: 
  #                  SILLY_CODE <chr> = the silly with IDs removed 
  #                  IDs <list> = the IDs removed
  #                  is_empty <lgl> = were all IDs removed?
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  password <- "************"
  #
  #  CreateLocusControl.GCL(markersuite = "Sockeye2011_96SNPs", username = "awbarclay", password = password)
  #  sillyvec <- c("SMCDO03", "SNEVA13")
  #  LOKI2R.GCL(sillyvec = sillyvec, username = "awbarclay", password = password)
  #  RemoveIndMissLoci.GCL(sillyvec = sillyvec)
  #
  #  RemoveIDs.GCL(silly = "SMCDO03", IDs = 1:10)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse)  # Install packages, if not in library and then load them.
  
  if(purrr::is_empty(silly)){
    
    stop( paste0("Silly is empty - is this expected?")) # check to see if silly contains anything
    
  }
  
  my.gcl <- get(paste0(silly, ".gcl"), pos = 1)
  
  if(purrr::is_character(IDs)) {IDs <- as.numeric(IDs)}
  
  if(!all(IDs %in% my.gcl$FK_FISH_ID)){
    
    stop(paste0("No IDs were removed. Some of the IDs do not exist for ", silly, "."))
    
  }

  assign(paste(silly, ".gcl", sep = ""),
         my.gcl %>% dplyr::filter(!FK_FISH_ID %in% IDs),
         pos = 1)
  
  message(paste0(length(IDs), " IDs were removed from ", silly, ".gcl"))
  
  if(dim(get(paste0(silly, ".gcl"), pos = 1))[1]==0){
    
    warning(paste0("All IDs were removed from ", silly, ".gcl"))
    
  }
  
  return(tibble::tibble(SILLY_CODE = silly, IDs_Removed= list(IDs) %>% purrr::set_names(silly), is_empty = dim(get(paste0(silly, ".gcl"), pos = 1))[1]==0))
  
}