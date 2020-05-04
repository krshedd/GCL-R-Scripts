old2new_LocCtrl.GCL <- function (LocCtrl = LocusControl, overwrite = FALSE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function converts the old style LocusControl object from a nested list to a tibble.
  # The old style LocusControl object can be saved as LocusControl_old before it is overwritten.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   LocCtrl - the locus control object
  #   overwrite - logical; whether you want the old LocusControl object overwritten without saving 
  #              (TRUE) or assign LocusControl to LocusControl_old before overwriting (FALSE) 
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function assigns a converted LocusControl object to the current workspace.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # old2new_LocCtrl.GCL(LocCtrl = LocusControl, overwrite = TRUE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); p_load(tidyverse) #Install packages, if not in library and then load them.
  
  if(tibble::is_tibble(LocCtrl)){
    
    stop("The LocusControl object is already a tibble.")
    
    }
  
  if(!overwrite){
    
    assign(paste0("LocusControl_old"), value = LocCtrl)
    
    } #Saving old gcl
  
  locus_info <- with(LocCtrl, tibble(MarkerSuite, locusnames, Publishedlocusnames, nalleles, ploidy))
  
  alleles <- tibble::tibble(alleles = lapply(LocusControl$locusnames, function(locus){
    
    a <- LocusControl$alleles[[locus]]
    
    tibble::tibble(allele = seq(length(a)), call = a)
    
  }) %>% purr::set_names(LocusControl$locusnames)
  
  ) 
  
  LocCtrl_tidy = dplyr::bind_cols(locus_info, alleles)
  
  assign("LocusControl", value = LocCtrl_tidy)
  
  message("LocusControl converted from a list to a tibble")
  
  return(NULL)
  
}