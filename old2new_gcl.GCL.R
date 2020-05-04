old2new_gcl.GCL <- function (sillyvec, save_old = FALSE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function converts the old style "*.gcl" objects from nested arrays to a tibbles.
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
  # 
  # old2new_gcl.GCL(sillyvec = sillyvec, overwrite = TRUE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); p_load(tidyverse) #Install packages, if not in library and then load them.
  
  sillyvec0 <- sillyvec[sapply(sillyvec, function(silly){ #Excluding objects that are already in tidy format.
    
    !tibble::is_tibble(get(paste0(silly, ".gcl")))
    
  })]
  
  sapply(sillyvec0, function(silly){
    
    my.gcl <- get(paste0(silly, ".gcl"))
    
    if(save_old){assign(paste0(silly, ".gcl_old"), value = my.gcl, pos = -1, envir = .GlobalEnv)} #Saving old gcl
    
    attr <- my.gcl$attributes
    
    scores <- lapply(seq(dim(my.gcl$scores)[[3]]), function(dim){
      
      s = my.gcl$scores[ , , dim, drop = FALSE] %>% 
        tibble::as_tibble()  
      
      if(dim > 1){
        
        colnames(s) <- paste(colnames(s), dim-1, sep =".")
        
      }
      
      s
      
    }) %>% 
      dplyr::bind_cols() 
    
    tidy.gcl <- scores %>%
      dplyr::bind_cols(attr) %>% 
      dplyr::select(colnames(attr), sort(colnames(scores)))
    
    assign(paste0(silly, ".gcl"), value = tidy.gcl, pos = -1, envir = .GlobalEnv)
    
  })#End silly
  
  message(paste0("The following *.gcl objects have been converted to tibbles: ", paste0(sillyvec, collapse = " ,")))
  
}