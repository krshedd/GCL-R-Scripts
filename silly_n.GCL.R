silly_n.GCL <- function(sillyvec){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function gets the sample size of .gcl objects.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   sillyvec - A character vector of silly codes with associated *.gcl objects that you want sample sizes for.
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function returns a two-column tibble of sillys (silly) and sample sizes (n)
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  #
  # SampleSize.GCL(sillyvec = sillyvec157)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); p_load(tidyverse) #Install packages, if not in library and then load them.
  
  lapply(sillyvec, function(silly){
    
    my.gcl <- get(paste0(silly, ".gcl"))
    
    if(!tibble::is_tibble(my.gcl)){
      
      stop("Make sure all of the sillys in sillyvec are in the new-style tibble format")
      
    }
    
    tibble::tibble(silly = silly, n = dim(my.gcl)[[1]])
    
  }) %>% 
    bind_rows()
  
}