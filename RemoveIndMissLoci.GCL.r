RemoveIndMissLoci.GCL <- function(sillyvec, loci = NULL, proportion = 0.8){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function removes individuals from "*.gcl" objects that have fewer non-missing loci than that specified by "proportion".
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   proportion - the cut-off proportion of the number of non-missing loci.
  #
  #   loci - optional vector of locus names to be considered when removing individuals for missing data; default is LocusControl$locusnames
  #          This argument is useful if you have loci in your markersuite that do not perform well in the lab and they will be dropped from the final dataset. 
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Returns a tibble of indiduals removed by silly.
  #    Assigns the ".gcl" objects back to workspace after removing individuals with missing loci.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  # 
  #   missloci_ind <- RemoveIndMissLoci.GCL(sillyvec = sillyvec157, proportion = 0.8)
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # This function only removes individuals for failed scores (zero scores), not unanalyzed loci (NAs).
  # Therefore, individuals with a lot of NAs will not be removed. Use RemoveNAindv.GCL to remove individuals with NAs for all loci.
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!exists("LocusControl")){
    
    stop("'LocusControl' not yet built.")
    
  }
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse) #Install packages, if not in library and then load them.
  
  if(is.null(loci)){
    
    loci <- LocusControl$locusnames #Use all loci if loci argument set to NULL.
    
  }
  
  if(!length(loci)==sum(loci%in%LocusControl$locusnames)){
    
   stop("The following loci are not in LocusControl-", paste0(loci[- na.omit(match(LocusControl$locusnames, loci))], collapse = " ,")) #Making sure loci are spelled correctly.
    
  }
  
  output <- lapply(sillyvec, function(silly){
    
    my.gcl <- get(paste0(silly, ".gcl"), pos = 1)
    
    tmp <- my.gcl %>% 
      dplyr::select(tidyselect::all_of(loci)) 
    
    IDsToRemove <- my.gcl %>% 
      dplyr::mutate(nloci = rowSums(!is.na(tmp)), nmissing = rowSums(tmp == "0")) %>% 
      dplyr::mutate(prop_loci = 1-(nmissing/nloci)) %>% 
      dplyr::select(prop_loci, tidyselect::everything()) %>% 
      dplyr::filter(prop_loci <= proportion) %>% 
      dplyr::pull(FK_FISH_ID)
    
    assign(x = paste0(silly, ".gcl"), value = my.gcl %>% dplyr::filter(!FK_FISH_ID%in%IDsToRemove), pos = 1, envir = .GlobalEnv )
    
    tibble::tibble(SILLY_CODE = silly, IDs_Removed = IDsToRemove)
    
  }) %>% 
    dplyr::bind_rows()
  
 if(dim(output)[1]==0){
   
  message("No individuals were removed")
   
 } else{
   
   message(paste0("A total of ", dim(output)[1]), " individuals were removed from sillys in sillyvec.")
   
 }

  return(output)

}
