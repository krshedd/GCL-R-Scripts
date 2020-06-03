gcl2Genepop.GCL <- function(sillyvec, loci, path, VialNums = TRUE, usat = FALSE, ncores = 4){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function create a GENEPOP input file from "*.gcl" objects.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - character vector of SILLYs - example: c(KCRES10", "KSTRA10", "KNIKOL12", "KNIKOL13")
  #
  #   loci - character vector of locus names as they are spelled in LOKI - example: c("GTH2B-550", "NOD1", "Ots_100884-287")
  #
  #   path - full file path to write out the GENEPOP file with "\\" or "/" separator between folders 
  #           example: "V:\\Analysis\\2_Central\\Chinook\\Cook Inlet\\2019\\2019_UCI_Chinook_baseline_hap_data\\GENEPOP\\genepopfile.gen"
  #                   or "V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/GENEPOP/genepopfile.gen"
  # 
  #   VialNums - logical; if TRUE (default), vial numbers will be included for each individual next to their silly code separated by an underscore (e.g. KCRESC10_1)
  #                       if FALSE, only the silly code will be included for each individual.
  #
  #   usat - logical; whether the data are from microsatellites (TRUE) or not (FALSE). 
  #          This is included because GENEPOP only accepts numberic alleles; SNP alleles have to be converted from character to numeric and microsatellite alleles are already numeric.
  #
  #   ncores - the number of cores for mulitcoring using doParallel and foreach. 
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Writes out a GENEPOP file.
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #  
  #  ### SNP
  #  attach("V:/Analysis/2_Central/Coho/Cook Inlet/2019/2019_Cook_Inlet_coho_baseline/2019_Cook_Inlet_coho_baseline.RData")
  #  old2new_LocCtrl.GCL()
  #  sapply(sillyvec104, function(silly){assign(paste0(silly, ".gcl"), get(paste0(silly, ".gcl")), pos = -1, envir = .GlobalEnv)})
  #  sillyvec <- sillyvec104
  #  loci <- loci81
  #  detach()
  #  old2new_gcl.GCL(sillyvec)
  #  gcl2Genepop.GCL(sillyvec = sillyvec, loci = loci, path = "V:/Analysis/2_Central/Coho/Cook Inlet/2019/2019_Cook_Inlet_coho_baseline/GENEPOP/genepopfile.gen", VialNums = TRUE, usat = FALSE, ncores = 4)
  #
  #  ### uSAT
  #  attach("V:/Analysis/2_Central/Chinook/Susitna-Watana/2017/Upper and Middle River Analysis/SuWuUpperMiddleRiverAnalysis_tidy.RData")
  #  sapply(Sillys, function(silly){assign(paste0(silly, ".gcl"), get(paste0(silly, ".gcl")), pos = -1, envir = .GlobalEnv)})
  #  sillyvec <- Sillys
  #  loci <- loci13
  #  LocusControl <- LocusControl
  #  detach()
  #  gcl2Genepop.GCL(sillyvec = sillyvec, loci = loci, path = "V:/Analysis/2_Central/Chinook/Susitna-Watana/2017/Upper and Middle River Analysis/Genepop/genepopfile.gen", VialNums = TRUE, usat = TRUE, ncores = 4)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  start_time <- Sys.time()
  
  if(!exists("LocusControl")){
    
    stop("'LocusControl' not yet built.")
    
  }
  
  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse, doParallel, foreach) #Install packages, if not in library and then load them.
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  alleles <- LocusControl$alleles[loci] %>% 
    bind_rows(.id = "locus")
  
  ploidy <- LocusControl %>% 
    filter(locusnames%in%loci) %>% 
    pull(ploidy) %>% 
    set_names(loci)
  
  if(sum(ploidy == 1)){
    
    stop("One or more of the loci supplied are haploid. GENEPOP only works with diploid data (ploidy = 2).")
    
  }
  
  nalleles <- LocusControl %>% 
    filter(locusnames%in%loci) %>% 
    pull(nalleles) %>% 
    set_names(loci)
  
  my.gcl <- lapply(sillyvec, function(silly){
    
    get(paste(silly, ".gcl", sep=""), pos = 1)
    
  }) %>% purrr::set_names(sillyvec)
  
  file <- "GENEPOP input format"
  
  file <- rbind(file, cbind(loci))
  
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  # multicore loop pop scores
  scores_all <- foreach::foreach(silly = sillyvec, .packages = c("tidyverse", "tidyselect")) %dopar% {
    
    new.gcl <- my.gcl[[silly]]
    
    IDs <- paste(new.gcl$SILLY_CODE, new.gcl$FK_FISH_ID, sep = "_") 
    
    if(!VialNums){ 
      
      vials <- new.gcl$SILLY_CODE
      
    } else{vials <- IDs}
    
    scores_names <- c(loci, paste(loci, "1", sep =".")) %>% 
      sort()
    
    scores <- new.gcl %>% 
      dplyr::select(tidyselect::all_of(scores_names)) %>% 
      as.data.frame(stringsAsFactors = FALSE)
    
    dimnames(scores)[[1]] = IDs 
    
    if(usat) {
      
      maxchar <- max(nchar(alleles$call))
      
      pop_scores <- lapply(loci, function(loc){
        
        variables <- c(loc, paste(loc, 1, sep = "."))
        
        scores %>%
          dplyr::select(tidyselect::all_of(variables)) %>%
          tidyr::replace_na(replace = list(0, 0) %>% 
                              purrr::set_names(variables)) %>%
          dplyr::mutate_all(.funs = stringr::str_pad, width = maxchar, pad = "0", side = "left") %>% 
          tidyr::unite(col = loc, c(loc, paste(loc, 1, sep = ".")), sep = "")
        
      }) %>% 
        dplyr::bind_cols() %>% 
        purrr::set_names(loci) %>% 
        dplyr::mutate(ID = vials, comma = ",") %>% 
        tidyr::unite("ID_comma", c("ID", "comma"), sep = " ") %>% 
        tidyr::unite("comb_scores", tidyselect::all_of(loci), sep = " ") %>% 
        tidyr::unite("ID_comma_scores", c("ID_comma","comb_scores"), sep = "  ") %>% 
        dplyr::pull(ID_comma_scores) 
      
      rbind("Pop", cbind(pop_scores))
      
    } else {
      
      pop_scores <- lapply(loci, function(loc){
        
        variables <- c(loc, paste(loc, 1, sep = "."))
        
        my.alleles <- alleles %>% 
          dplyr::filter(locus==loc)
        
        scores %>%
          dplyr::select(tidyselect::all_of(variables)) %>% 
          dplyr::mutate_all(.funs = factor, levels = my.alleles$call) %>% 
          dplyr::mutate_all(.funs = as.numeric) %>% 
          dplyr::mutate_all(.funs = as.character) %>%
          tidyr::replace_na(replace = list(0, 0) %>% 
                              set_names(variables)) %>%
          dplyr::mutate_all(.funs = stringr::str_pad, width = 2, pad = "0", side = "left") %>% 
          tidyr::unite(col = loc, c(loc, paste(loc, 1, sep = ".")), sep = "")
        
      }) %>% 
        dplyr::bind_cols() %>% 
        purrr::set_names(loci) %>% 
        dplyr::mutate(ID = vials, comma = ",") %>% 
        tidyr::unite("ID_comma", c("ID", "comma"), sep = " ") %>% 
        tidyr::unite("comb_scores", tidyselect::all_of(loci), sep = " ") %>% 
        tidyr::unite("ID_comma_scores", c("ID_comma","comb_scores"), sep = "  ") %>% 
        dplyr::pull(ID_comma_scores) 
      
      rbind("Pop", cbind(pop_scores))
      
    }
    
  }#End multicore loop
  
  parallel::stopCluster(cl)
  
  file <- rbind(file, scores_all %>% unlist() %>% cbind())
  
  utils::write.table(file, path, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  print(Sys.time() - start_time)
  
  return(NULL) 
  
}