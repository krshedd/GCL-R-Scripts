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
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  CreateLocusControl.GCL(markersuite = "Sockeye2011_96SNPs", username = "awbarclay", password = password)
  #  sillyvec = c("SMCDO03", "SNEVA13")
  #  password = "************"
  #  LOKI2R.GCL(sillyvec = sillyvec, username = "awbarclay", password = password)
  #  RemoveIndMissLoci.GCL(sillyvec = sillyvec)
  #
  #  RemoveIDs.GCL(silly = "SMCDO03", IDs = 1:10)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse)  # Install packages, if not in library and then load them.
  
  my.gcl <- get(paste0(silly, ".gcl"), pos = 1)
  
  if(purrr::is_character(IDs)) {IDs <- as.numeric(IDs)}
  
  if(!all(IDs %in% my.gcl$FK_FISH_ID)){
    
    stop(paste0("No IDs were removed. Some of the IDs do not exist for ", silly, "."))
    
  }

  assign(paste(silly, ".gcl", sep = ""),
         my.gcl %>% dplyr::filter(!FK_FISH_ID %in% IDs),
         pos = 1)
  
  message(paste0(length(IDs), " IDs were removed from ", silly, ".gcl"))
  
}