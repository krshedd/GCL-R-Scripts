create_bayes_mixture.GCL <- function(mixvec, loci, dir, ncores = 4){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function creates a BAYES mixure (.mix) file for each silly in mixvec.
  #   **Note: If you want to analyze more than one silly as a mixture, use PoolCollections.GCL to combine them into a new silly.gcl**
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   mixvec - character vector of ".gcl" objects you want to produce mixture files for.
  #   loci - character vector of the loci you wish to include
  #   dir - character vector of where to save the ".bse" file
  #            
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns the fortran format of the mixture file - this object is needed for create_bayes_control.GCL
  #   Saves each mixture in mixvec as a .mix file
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  # LOKI2R.GCL(sillyvec = c("KSUSC18FW", "KSUSCN18", "KSUSC19FW", "KSUSCN19"), username = "awbarclay", password = password)
  # CombineLoci.GCL(sillyvec = c("KSUSC18FW", "KSUSCN18", "KSUSC19FW", "KSUSCN19"), markerset = c("Ots_MHC1", "Ots_MHC2"))
  # PoolCollections.GCL(collections = c("KSUSC18FW", "KSUSCN18"), newname = "Susitna2018")
  # PoolCollections.GCL(collections = c("KSUSC19FW", "KSUSCN19"), newname = "Susitna2019")
  # loci <- c(loci82, "Ots_MHC1.Ots_MHC2")
  # mix_fortran <- create_bayes_mixture.GCL(mixvec = c("Susitna2018", "Susitna2019"), loci = loci, dir = getwd(), ncores = 8)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse, doParallel, parallel, foreach) #Install packages, if not in library and then load them.
  
  if(sum(is.na(match(loci, LocusControl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci, LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  start_time <- Sys.time()
  
  my.gcl <- sapply(mixvec, function(silly){
    
    get(paste(silly, ".gcl", sep = ""), pos = 1)
    
  }, simplify = FALSE)
   
  for(silly in mixvec){ #Start silly loop
    
    new.gcl <- my.gcl[[silly]]
  
    scores_names <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
      as.vector() 
    
    scores0 <- new.gcl %>% 
      dplyr::select(tidyselect::all_of(scores_names)) %>% 
      as.data.frame(stringsAsFactors = FALSE)
    
    id <- as.numeric(dimnames(scores0)[[1]])
    
    cl <- parallel::makePSOCKcluster(ncores)
    
    doParallel::registerDoParallel(cl, cores = ncores)  
    
    #Had to suppress an anoying message which only printed to console when running the foreach loop.
    suppressMessages(
      
      mix_scores <-  foreach::foreach(loc = loci, .packages = c("tidyverse", "tidyselect", "janitor"), .export = c("LocusControl"), .combine = dplyr::bind_cols) %dopar% { #Start multicore loop: loci
      
      variables <- c(loc, paste(loc, 1, sep = "."))
        
      my.alleles <- LocusControl$alleles[[loc]]
        
      comb_alleles0 <- scores0 %>%
        dplyr::select(tidyselect::all_of(variables)) %>% 
        dplyr::mutate(dplyr::across(dplyr::everything(), .fns = factor, levels = my.alleles$call)) %>% 
        dplyr::mutate(dplyr::across(dplyr::everything(), .fns = as.numeric)) %>% 
        tibble::as_tibble() %>% 
        dplyr::mutate(id = id) %>% 
        tidyr::pivot_longer(-id) %>% 
        dplyr::mutate(value = factor(value, levels = my.alleles$allele)) %>% 
        janitor::tabyl(id, value, show_missing_levels = TRUE, show_na = TRUE)
      
      comb_alleles <-  comb_alleles0 %>% 
        dplyr::select(dplyr::all_of(my.alleles$allele %>% as.character())) %>% 
        tidyr::unite(col = "comb", as.character(my.alleles$allele), sep = "") %>%
        dplyr::pull(comb) %>% 
        stringr::str_pad(width = length(my.alleles$allele)+1, pad = " ", side = "left") %>% 
        tibble::as_tibble()
      
      comb_alleles
      
      } %>% 
      tidyr::unite(col = "mix", dplyr::everything(), sep = "") #End multicore loop: loci
    
      ) #Suppress
    
    write.table(mix_scores, file = paste0(dir, "/", silly, ".mix"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    parallel::stopCluster(cl)
      
    } #End silly loop
  
  a_rle <- sapply(loci, function(loc){LocusControl$nalleles[[loc]]}) %>% 
    rle()
  
  mix_fortran <- paste0("(", paste0(sapply(seq(length(a_rle$values)), function(locus){
    
    paste(a_rle$lengths[locus], "(1X,", a_rle$values[locus], "I1)")
    
    }), collapse = ","), ")")
  
  Sys.time()-start_time
  
  return(mix_fortran)
  
}