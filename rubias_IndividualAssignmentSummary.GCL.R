rubias_IndividualAssignmentSummary.GCL <- function(rubias_output = NULL, mixnames = NULL, path = "rubias/output"){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function will summarize rubias individual assignments.
  # 
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # rubias_output - output list object from `run_rubias_mixture` or `infer_mixture`
  #
  # mixnames - a character vector of mixture names to include in the individual assignment summary. 
  #
  # path - character vector of where to find output from each mixture as a .csv (created by `run_rubias_mixture`)
  #
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   A tibble with the following variables:
  #     mixture_collection - the name of the mixture;
  #     indiv - the individual (aka SillySource);
  #     and a variable for each group in group_names containing the proportion of MCMC iterations
  #     in which an individual was assigned to a specific reporting group.
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # path <- "C:/Users/awbarclay/Documents/Analysis/Sockeye/UCI_sockeye_2020_postseason/rubias/output"
  # 
  # rubias_output <- readRDS(file = "C:/Users/awbarclay/Documents/Analysis/Sockeye/UCI_sockeye_2020_postseason/output/rubias_output.rds")
  # 
  # mixnames <- c("DriftDW_20", "DriftCorr_20", "UpperSub_20", "Kasilof600ft_20", "WestKalgin_20", "Eastern_20", "General_North_20", "General_South_20")
  # 
  # rubias_IndividualAssignmentSummary.GCL(rubias_output = rubias_output, mixnames = mixnames)
  #  
  # rubias_IndividualAssignmentSummary.GCL(rubias_output = NULL, mixnames = mixnames, path = path)
  #     
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse) #Install packages, if not in library and then load them.
  
  if(is.null(rubias_output)){
    
    if(is.null(mixnames)){
      
      stop("If rubias_output is NULL, then a vector of mixnames must be supplied for reading in the individual probabilities for each mixture.")
    
      }else{
        
      file_check <- sapply(mixnames, function(mix){
        
        !file.exists(paste0(path, "/", mix, "_indiv_posteriors.csv"))
        
      })
      
      if(sum(file_check) > 0){
        
        missing <- names(file_check[file_check])
        
        stop("No individual posterior files exist in the supplied path for the following mixtures: ", paste0(missing, collapse = ", "))
        
      }
      
     indiv_post <- lapply(mixnames, function(mix){
        
        read_csv(file = paste0(path, "/", mix, "_indiv_posteriors.csv"), col_types = cols(
          mixture_collection = col_character(),
          indiv = col_character(),
          repunit = col_character(),
          collection = col_character(),
          PofZ = col_double(),
          log_likelihood = col_double(),
          z_score = col_double(),
          n_non_miss_loci = col_double(),
          n_miss_loci = col_double()
        ))
        
      }) %>% bind_rows()
      
      results <- indiv_post%>% 
        mutate(repunit = factor(repunit, levels = unique(repunit))) %>% 
        group_by(mixture_collection, indiv, repunit) %>% 
        summarize(prob = sum(PofZ), .groups = "drop") %>% 
        pivot_wider(names_from = repunit, values_from = prob)
      
    }
    
  } else{
    
    indiv_post <- rubias_output$indiv_posteriors
    
    mix_check <- !mixnames == indiv_post$mixture_collection %>% unique() %>% set_names(mixnames)
    
   if(sum(mix_check) > 0){
     
     stop("The folling mixnames are not included in the supplied rubias output object: ", paste0(names(mix_check[mix_check]), collapse = ", "))
     
   }
    
    results <- rubias_output$indiv_posteriors %>% 
      filter(mixture_collection %in% mixnames) %>% 
      mutate(repunit = factor(repunit, levels = unique(repunit))) %>% 
      group_by(mixture_collection, indiv, repunit) %>% 
      summarize(prob = sum(PofZ), .groups = "drop") %>% 
      pivot_wider(names_from = repunit, values_from = prob)
    
  }
  
  return(results)
  
}
  
 