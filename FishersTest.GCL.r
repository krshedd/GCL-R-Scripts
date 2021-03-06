FishersTest.GCL <- function(freq, loci, tests){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function is a wrapper for FisherCompute.GCL() to run Fisher's Exact Tests for homogeneity of allele frequencies for multiple sets of collections.
  #   Null hypothesis: there are no significant differences in allele freqencies.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   freq - a tibble of allele freqencies produced by FreqPop.GCL() for the collections (sillys) you want to test.
  #
  #   loci - a character vector of locus names
  #   
  #   tests - a list of character vectors. The collections is each element of the list are tested separately. 
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Produces a tibble of overall pvalues produced by combining individual locus pvalues 
  #    using Fisher's method and a nested tibble of pvalues by locus for each test.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  #   old2new_LocCtrl.GCL()
  #   old2new_gcl.GCL(sillyvec = c("KKILL05","KKILL06", "KFUNN05", "KFUNN06"), save_old = FALSE)
  #
  #   freq <- FreqPop.GCL(sillyvec = c("KKILL05","KKILL06", "KFUNN05", "KFUNN06"), loci = loci443)
  # 
  #   FishersTest.GCL(freq = freq, loci = loci443, tests = list(c("KKILL05","KKILL06"), c("KFUNN05", "KFUNN06")))
  # 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse) #Install packages, if not in library and then load them.

  if(sum(is.na(match(loci, LocusControl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci,LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
    }
  
  testnames <- sapply(tests, function(test){
    
    paste(test, collapse = ".")}
    
    )
  
  #Looping through tests with lapply - may want to use foreach %dopar% if this is too slow for people.
  output <- lapply(tests, function(test){
    
    my.test.freq <- freq %>% 
      dplyr::filter(silly%in%test)
    
    FisherCompute.GCL(freq = my.test.freq, loci = loci)
    
    }) %>% 
    dplyr::bind_rows() %>%
    dplyr::mutate(bylocus = bylocus %>% purrr::set_names(testnames))
    
  return(output)
  
}