FreqPop.GCL <- function(sillyvec, loci, ncores = 4 ){
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function gets the allele frequency for each locus in loci for each silly in sillyvec.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   loci - vector of locus names; if set to NULL all loci in the ".gcl" obejects will be used.
  #
  #   ncores - the number of cores to use in a foreach %dopar% loop. If the nubmer of core exceeds the number on your device, then ncores defaults to detectCores()
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  A tibble with the following variables: silly, locus, allele_no (from LocusControl), freq (allele frequency), and proportion
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  #  old2new_LocCtrl.GCL()  
  #  old2new_gcl.GCL(sillyvec67)
  #
  #  Freq <- FreqPop.GCL(sillyvec = sillyvec67, loci = loci413)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  start.time <- Sys.time() 
  
  if(sum(is.na(match(loci, LocusControl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci,LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  if(ncores > detectCores()){
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
  
    }
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse, doParallel, parallel, foreach)  # Install packages, if not in library and then load them

  all.gcl <- lapply(sillyvec, function(silly){get(paste(silly, ".gcl", sep = ""), pos = 1)}) %>% 
    purrr::set_names(sillyvec)
  
  scores_cols <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
    as.vector() 
  
  alleles <- LocusControl$alleles[loci] %>% 
    dplyr::bind_rows(.id = "locus") %>% 
    dplyr::rename(allele_no = allele)
  
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  

  #Start parallel loop
  freqs <- foreach::foreach(silly = sillyvec, .packages = c("tidyverse")) %dopar% {
    
    my.gcl <- all.gcl[[silly]]
    
    my.gcl %>% 
      dplyr::select(all_of(scores_cols)) %>% 
      tidyr::gather(key = "locus", value = "allele") %>% 
      dplyr::filter(locus%in%loci) %>%
      dplyr::full_join(alleles, by = c("locus"="locus", "allele"="call"), keep = TRUE) %>% 
      dplyr::mutate(locus = gsub("\\.1$", "", locus.y)) %>% 
      dplyr::group_by(locus, allele_no, call) %>% 
      dplyr::summarize(freq = sum(!is.na(allele)), .groups = "drop_last") %>% 
      dplyr::mutate(silly = !!silly) %>% 
      dplyr::filter(!is.na(allele_no))
   
  } %>% 
    dplyr::bind_rows() %>% 
    dplyr::relocate(silly, before = locus)
  
  parallel::stopCluster(cl) #End parallel loop
  
  output <- freqs %>% 
    dplyr::rename(allele = call) %>% 
    dplyr::group_by(silly, locus) %>% 
    dplyr::mutate(total = sum(freq)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(proportion = freq/total) %>% 
    dplyr::select(-total)
  
  print(Sys.time()-start.time)
  
  return(output)
  
}
