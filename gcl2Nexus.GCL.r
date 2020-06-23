gcl2Nexus.GCL <- function(sillyvec, loci, path, VialNums = TRUE, PopNames = NULL, ncores = 4){

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function creates a Nexus file in GDA input format from "*.gcl" objects.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - character vector of SILLYs - example: c(KCRES10", "KSTRA10", "KNIKOL12", "KNIKOL13")
  #
  #   loci - character vector of locus names as they are spelled in LOKI - example: c("GTH2B-550", "NOD1", "Ots_100884-287")
  #
  #   path - full file path to write out the GENEPOP file with "\\" or "/" separator between folders 
  #           example: "V:\\Analysis\\2_Central\\Chinook\\Cook Inlet\\2019\\2019_UCI_Chinook_baseline_hap_data\\output\\nexusfile.nex"
  #                   or "V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/output/nexusfile.nex"
  # 
  #   VialNums - logical; if TRUE (default), vial numbers will be included for each individual next to their silly code separated by an underscore (e.g. KCRESC10_1)
  #                       if FALSE, only the silly code will be included for each individual.
  #
  #   PopNames - a character vector the same length as sillyvec to give populations new names. If no vector is given, PopNames defaults to "sillyvec". 
  #
  #   ncores - the number of cores for mulitcoring using doParallel and foreach. 
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Writes out a nexus file.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # 
  #   load("V:/Analysis/2_Central/Sockeye/Cook Inlet/2012 Baseline/Baseline/CI2012Baseline.RData")
  #   CreateLocusControl.GCL(markersuite = "Sockeye2011_96SNPs", locusnames = NULL, username = "awbarclay", password = password) 
  #   old2new_gcl.GCL(sillyvec = PooledNames71)
  #  
  #   gcl2Nexus.GCL(sillyvec = PooledNames71, loci = loci96, path = "V:/Analysis/2_Central/Sockeye/Cook Inlet/2012 Baseline/Baseline/Output/nexusfile.nex", VialNums = TRUE, PopNames = CINames71, ncores = 4)
  #
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
  
  ploidy <- LocusControl %>% 
    filter(locusnames%in%loci) %>% 
    pull(ploidy) %>% 
    set_names(loci) %>% 
    sort(decreasing = TRUE)
  
  loci <- names(ploidy)
  
  alleles <- LocusControl$alleles[loci] %>% 
    bind_rows(.id = "locus")
  
  nalleles <- LocusControl %>% 
    filter(locusnames%in%loci) %>% 
    pull(nalleles) %>% 
    set_names(loci)
  
  nloci <- length(loci)
  
  my.gcl <- lapply(sillyvec, function(silly){
    
    get(paste(silly, ".gcl", sep=""), pos = 1)
    
  }) %>% purrr::set_names(sillyvec)
  
  if(is.null(PopNames)){PopNames <- sillyvec}
  
  if(!length(PopNames)==length(sillyvec)){stop("PopNames is not the same length as sillyvec.")}
  
  PopNames <- PopNames %>% purrr::set_names(sillyvec)
  
  file <- "#nexus"
  
  file <- rbind(file,"")
  
  file <- rbind(file, "begin gdadata; [!GDA input format]")
  
  file <- rbind(file, paste0("dimensions npops=",length(sillyvec), "  nloci=", nloci, ";"))
  
  file <- rbind(file, paste0("format missing=? separator=/;"))
  
  if(min(ploidy)==1){
    
    hapset0 <- grep("1", ploidy)
    
    if(length(hapset0)>1){hapset <- paste(range(hapset0), collapse="-")}
    
    file <- rbind(file, paste0("hapset ", hapset, ";"))
  }
  
  file <- rbind(file, "locusallelelabels")
  
  max.char <- max(nchar(seq(nloci)))
  
  locuslabels <- cbind(c(paste0("   ", format(seq(nloci-1), width = max.char), "  ", loci[seq(nloci-1)], ","),
                         paste0("   ", format(nloci, width = max.char), "  ", loci[nloci]))
                       )
  
  file <- rbind(file, locuslabels)
  
  file <- rbind(file, cbind(c(";", "matrix")))
  
  my.gcl <- lapply(sillyvec, function(silly){
    
    get(paste(silly, ".gcl", sep=""), pos = 1)
    
  }) %>% purrr::set_names(sillyvec)
  
  if(exists("hapset")){
    
    scores_names <- sapply(loci[-hapset0], function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
      as.vector() 
    
  } else{
    
    scores_names <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
      as.vector() 
    
  }
  
  # Create scores tables
  
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  scores_all <- foreach::foreach(silly = sillyvec, .packages = c("tidyverse", "tidyselect")) %dopar% {
    
    new.gcl <- my.gcl[[silly]]
    
    IDs <- paste(new.gcl$SILLY_CODE, new.gcl$FK_FISH_ID, sep = "_") 
    
    if(!VialNums){ 
      
      vials <- new.gcl$SILLY_CODE
      
    } else{vials <- IDs}
 
    scores <- new.gcl %>% 
      dplyr::select(tidyselect::all_of(scores_names)) %>% 
      as.data.frame(stringsAsFactors = FALSE)
    
    dimnames(scores)[[1]] = IDs 
    
    scores[scores==0] <- "?"
      
    scores[is.na(scores)] <- "?"
    
    pop_scores <- lapply(loci[-hapset0], function(loc){
      
      variables <- c(loc, paste(loc, 1, sep = "."))
      
      scores %>%
        dplyr::select(tidyselect::all_of(variables))%>% 
        tidyr::unite(col = loc, variables, sep = "/  ")
      
    }) %>% 
      dplyr::bind_cols() %>% 
      dplyr::bind_cols(scores %>% select(all_of(loci[hapset0]))) %>% 
      purrr::set_names(loci) %>% 
      dplyr::mutate(ID = paste0(vials)) %>% 
      tidyr::unite("comb_scores", tidyselect::all_of(loci), sep = "   ") %>% 
      tidyr::unite("ID_scores", c("ID","comb_scores"), sep = "  ") %>% 
      dplyr::pull(ID_scores) 
    
    if(silly == dplyr::last(sillyvec)){
      
      rbind(paste0(PopNames[silly],":"), 
            cbind(pop_scores), 
            t(cbind(";", "end;")))
      
    } else{
      
      rbind(paste0(PopNames[silly],":"), 
            cbind(pop_scores), 
            t(cbind("", ",")))
      
      }
    
  } #End multicore loop
  
  parallel::stopCluster(cl)
  
  file <- rbind(file, scores_all %>% unlist() %>% cbind())
  
  write.table(file, path, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  Sys.time()-start_time
  
  return(NULL) 
 
}