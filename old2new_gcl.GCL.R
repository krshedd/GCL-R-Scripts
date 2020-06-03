old2new_gcl.GCL <- function (sillyvec, save_old = FALSE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function converts old style "*.gcl" objects from nested arrays to tibbles of scores and attributes. 
  # The old style objects can be saved as *.gcl_old before they are overwritten.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   sillyvec - A character vector of silly codes with associated *.gcl objects that need to be converted.
  #   overwrite - logical; whether you want the old objects overwritten without saving 
  #              (TRUE) or assign the new objects to *.gcl_old before overwriting (FALSE) 
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function assigns a converted *.gcl object to the current workspace.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  # 
  # old2new_gcl.GCL(sillyvec = sillyvec157, save_old = TRUE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse) #Install packages, if not in library and then load them.
  
  if(!all(sillyvec %in% stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))) {  # Do all sillys exist in the environment?
    
    missing_sillys <- setdiff(sillyvec, stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))
    
    stop(paste0("The following sillys are not in your environment:\n", paste0(missing_sillys, collapse = ", ")))
    
  }
  
  sillyvec0 <- sillyvec[sapply(sillyvec, function(silly){ #Excluding objects that are already in tidy format.
    
    !tibble::is_tibble(get(paste0(silly, ".gcl")))
    
  })]
  
  sapply(sillyvec0, function(silly){
    
    my.gcl <- get(paste0(silly, ".gcl"))
    
    if(save_old){assign(paste0(silly, ".gcl_old"), value = my.gcl, pos = -1, envir = .GlobalEnv)} #Saving old gcl
    
    attr <- my.gcl$attributes %>% 
      dplyr::mutate(SILLY_CODE = silly) #Some older attributes did not include SILLY_CODE
    
    scores <- lapply(seq(dim(my.gcl$scores)[[3]]), function(dim){
      
      s = my.gcl$scores[ , , dim, drop = FALSE] %>%
        tibble::as_tibble()
      
      colnames(s) <- dimnames(my.gcl$scores)[[2]] 
      
      if(dim > 1){
        
        colnames(s) <- paste(colnames(s), dim-1, sep =".")
        
      }
      
      s
      
    }) %>% 
      dplyr::bind_cols() %>% 
      na_if("0")
    
    tidy.gcl <- scores %>%
      dplyr::bind_cols(attr) %>% 
      dplyr::select(colnames(attr), sort(colnames(scores)))
    
    assign(paste0(silly, ".gcl"), value = tidy.gcl, pos = -1, envir = .GlobalEnv)
    
  })#End silly
  
  message(paste0("The following *.gcl objects have been converted to tibbles:\n", paste0(sillyvec, collapse = ", ")))
  
}