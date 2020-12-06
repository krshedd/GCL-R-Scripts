new2old.GCL <- function(sillyvec, save_new = FALSE, ncores = 4){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function converts new style "*.gcl" objects from tibbles to a list of arrays. 
  # The new style objects can be saved as *.gcl_new before they are overwritten.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   sillyvec - A character vector of silly codes with associated *.gcl objects that need to be converted.
  #   save_new - logical; whether you want the new objects overwritten without saving 
  #              (TRUE) or assign the new objects to *.gcl_new before overwriting (FALSE) 
  #   ncores - a numeric vector of length one indicating the number of cores to use
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function assigns a converted *.gcl object to the current workspace.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  #
  # new2old.GCL(sillyvec = sillyvec31, save_new = TRUE, ncores = 8)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse, janitor, foreach, doParallel, parallel)  # Install packages, if not in library and then load them.
  
  # Do all sillys exist in the environment?
  if(!all(sillyvec %in% stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))) {  
    
    missing_sillys <- setdiff(sillyvec, stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))
    
    stop(paste0("The following sillys are not in your environment:\n", paste0(missing_sillys, collapse = ", ")))
    
  }
  
  # Is there a LocusControl object in the environment?
  if(!exists(x = "LocusControl", where = 1)){
    
    stop("There is no LocusControl in the current environment.")
    
  }
  
  # Is LocusControl in the "old" list form?
  if(!tibble::is_tibble(LocusControl)) {
    
    stop("LocusControl is the 'old' style, list-form. This function requires the 'new' style, tibble.")
    
  }
  
  # Is uSAT?
  if(max(LocusControl$nalleles) > 2) {
    
    stop("Sorry, this function can only handle up to 2 alleles at the moment...")
    
  }
  
  # Checking to make sure all loci in .gcl object are included in the LocusControl object.
  locuscheck <-  sapply(sillyvec, function(silly){
    
    locvars <- names(get(paste0(silly, ".gcl")))[-c(1:19)]
    
    loc <- locvars[-grep(pattern = "\\.1$", x = locvars)]  # Need to run this regular expression by Chase to make sure it will always work.
    
    setdiff(loc, LocusControl$locusnames)
    
  }) %>% lapply(FUN = purrr::is_empty) %>% unlist() == FALSE
  
  if(sum(locuscheck) > 0){
    
    stop(paste0("The following sillys contain loci that are not included in LocusControl: ", paste0(sillyvec[locuscheck], collapse = ", ")))
    
  }
  
  # Excluding objects that are already in old format.
  sillyvec_ <- sillyvec[sapply(sillyvec, function(silly){  
    
    tibble::is_tibble(get(paste0(silly, ".gcl")))
    
  })]
  
  if(is_empty(sillyvec_)){
    
    stop("All of the sillys in sillyvec are already in the old-style format")
    
  }
  
  # Create list of gcl objects for parallel loop
  all.gcl <- sapply(sillyvec_, function(silly){get(paste0(silly, ".gcl"), pos = 1)}, simplify = FALSE)
  
  alleles <- LocusControl$alleles %>% 
    bind_rows(.id = "locus")
  
  ploidy <- LocusControl %>% 
    select(locusnames, ploidy) %>% 
    rename(locus = locusnames)
  
  all_ploid <- left_join(alleles, ploidy, by = "locus")
  
  # Start parallel loop
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  gcls <- foreach::foreach(silly = sillyvec_, .packages = c("tidyverse", "janitor"), .export = "LocusControl") %dopar% {
    
    my.gcl <- all.gcl[[silly]]
    
    locvars <- names(my.gcl[-c(1:19)])
    
    loci <- locvars[-grep(pattern = "\\.1$", x = locvars)]# Need to run this regular expression by Chase to make sure it will always work.
  
    ids <- my.gcl$FK_FISH_ID
    
    # Scores array
    dose1 <- my.gcl %>% 
      dplyr::select(dplyr::all_of(loci)) %>% 
      as.matrix()
    
    dimnames(dose1)[[1]] <- my.gcl$FK_FISH_ID
    
    dose2 <- my.gcl %>% 
      dplyr::select(dplyr::all_of(paste0(loci, ".1"))) %>% 
      as.matrix()
    
    dimnames(dose2) <- list(my.gcl$FK_FISH_ID, loci)
    
    s_dims <- c(lapply(dimnames(dose1), FUN = length), 2) %>% 
      unlist()
    
    scores <- array(c(dose1, dose2), dim = s_dims, dimnames = list(ids, dimnames(dose1)[[2]], c("Dose1", "Dose2")))
    
    # Counts array
    max_allele <- LocusControl$alleles %>% 
      dplyr::bind_rows() %>% 
      dplyr::pull(allele) %>% 
      max()
    
    c_dims <- c(lapply(dimnames(dose1), FUN = length), max_allele) %>% 
      unlist()
    
    counts <- array(data = NA, dim = c_dims, dimnames = list(ids, dimnames(dose1)[[2]], paste0("Allele ", 1:max_allele)) )
    
    tab <- dplyr::bind_rows(dplyr::bind_cols(id = ids, dose1), dplyr::bind_cols(id = as.numeric(dimnames(dose1)[[1]]), dose2), .id = "dose") %>% 
      tidyr::pivot_longer(cols = c(-dose, -id), names_to = "locus", values_to = "call") %>% 
      dplyr::left_join(all_ploid, by = c("locus", "call")) %>% 
      janitor::tabyl(id, locus, allele, show_na = TRUE)
    
    for(a in 1:max_allele){
      
      nms <- names(tab[[a]][ ,-1])
      
      tmp <- tab[[a]][ ,-1] %>% 
        as.matrix()
      
      dimnames(tmp)[[1]] <- ids
      
      mt <- all_ploid %>% filter(ploidy == 1)
      
      mtloci <- mt$locus %>% 
        unique()
      
      dploci <- setdiff(loci, mtloci)
      
      if(!is.null(tab[["NA_"]])){
        
        if(min(ploidy$ploidy) == 1){
        
          mt_na <- cbind(sapply(dploci, function(loc){rep(FALSE, length(ids))}), is.na(dose1[, mtloci]))[ , loci]
          
          tmp[mt_na] <- NA
          
        } 
        
        if(max(ploidy$ploidy) == 2){
          
          dp_na <- cbind(sapply(mtloci, function(loc){rep(FALSE, length(ids))}), tab$NA_[ , -1] > 0)[ , loci]
          
          tmp[dp_na] <- NA
          
        }
        
      }
      
      counts[,,a] <- tmp
      
    }
    
    # Attributes
    attributes <- my.gcl[ , 1:19] %>% 
      mutate(CAPTURE_DATE = CAPTURE_DATE %>% as.character() %>%  as.POSIXct(format = "%Y-%m-%d"), END_CAPTURE_DATE = END_CAPTURE_DATE %>% as.character() %>%  as.POSIXct(format = "%Y-%m-%d")) %>% 
      as.data.frame()
    
   list(counts = counts, scores = scores, n = length(ids), attributes = attributes)
  
  } %>% set_names(sillyvec_)
  
  parallel::stopCluster(cl)  # End parallel loop
  
  # If save_new == TRUE, save new style objects with extension
  if(save_new == TRUE){
    
    for(silly in sillyvec_){
      
      assign(x = paste0(silly, ".gcl_new"), value = get(paste0(silly, ".gcl")), pos = 1)
      
    }
    
  }
  
  # Assign old style objects to workspace
  for(silly in sillyvec_){
    
    assign(x = paste0(silly, ".gcl"), value = gcls[[silly]], pos = 1)
    
  }
  
  message(paste0("The following *.gcl objects have been converted to old-style lists:\n", paste0(sillyvec_, collapse = ", ")))
  
}
