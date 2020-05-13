CheckDupWithinSilly.GCL <- function(sillyvec, loci, quantile = 0.99, minproportion = 0.95, ncores = 8){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function checks for duplicate individuals within each silly in "sillyvec".
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   loci - vector of locus names; if set to NULL all loci in the ".gcl" obejects will be used.
  #
  #   quantile and minproportion - these arguments are used together to determine the cutt-off proportion at which a pair of duplicates 
  #                                is defined: i.e. proportion = max(quantile(duplication, quantile), minproportion. 
  #Setting "quantile" equal to NULL will skip the calculation of the duplication distribution and will run much faster--REQUIRES the package "reshape" ##### Need to use something other than reshape
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    When quantile is set to NULL, returns a tibble of duplicate pairs of individuals by silly.
  #    When quantile is a number, a list containing a tibble of duplicate pairs of individuals by silly and tibble the porportion of duplication for each pair of individuals.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  # 
  #   dupcheck <- CheckDupWithinSilly.GCL(sillyvec = sillyvec157, loci = loci557, quantile = 0.99, minproportion = 0.95, ncores = 4)
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Sillys without duplicates will be included in the output report but will have NAs for IDs1, IDs2, Missing1, Missing2, and proportion
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  start.time <- Sys.time() 
  
  if(sum(is.na(match(loci, LocusControl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci, LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep=""))
    
  }
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse, reshape, doParallel, foreach) #Install packages, if not in library and then load them.
  
  
  if(ncores > detectCores()) {
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
  }

  nsilly <- length(sillyvec)

  nloci <- length(loci)

  scores_cols <- c(loci, paste0(loci, ".1")) %>% 
    sort()
  
  my.gcl <- lapply(sillyvec, function(silly){
    
    get(paste(silly, ".gcl", sep=""), pos = 1)
    
    }) %>% purrr::set_names(sillyvec)
  
  #Start if NULL quantile
  if(is.null(quantile)){
    
    cl <- doParallel::makePSOCKcluster(ncores)
    
    doParallel::registerDoParallel(cl, cores=ncores)  
    
    # multicore loop
    
    dupcheck <- foreach::foreach(silly = sillyvec) %dopar% {
      
      new.gcl <- my.gcl[[silly]]
        
      IDs <- new.gcl$FK_FISH_ID
      
      n <- silly_n.GCL("new")$n
      
      nloci <- length(loci)
      
      #Combine dose 1 and 2 into single columns for each locus separated by a period.
      scores0 <- new.gcl[ ,c("FK_FISH_ID", scores_cols)] 
    
      if(n < 2){
        
        report <- tibble::tibble(silly = silly, ID1 = NA, ID2 = NA, Missing1 = NA, Missing2 = NA, proportion = NA) #If only one individual the results will be a tibble of NAs. Can't do next() using applys
        
      } else{
        
        scores1 <- sapply(IDs, function(id){
          
          ind <- scores0 %>% 
            filter(FK_FISH_ID==id)
          
          sapply(loci,function(locus){
            
            ifelse(is.na(ind[ ,locus]), NA, paste(ind[ , locus], ind[ , paste0(locus,".1")], sep = "."))
            
            })
          
          }) %>% 
          t() %>% 
          data.frame()
        
        dimnames(scores1)[[1]] <- IDs
        
        sort.scores.df <- reshape::sort_df(scores1)
        
        sortIDs <- dimnames(sort.scores.df)[[1]]
        
     #Calculate the porpotion of duplicated loci for each silly
      duplication <- t(sapply(1:(n-1), function(id){ 
            
        vec <- sapply(1:nloci, function(locus){ 
          
          sort.scores.df[id, locus]==sort.scores.df[id+1, locus]
          
        });
        
        ifelse(sum(!is.na(vec)), sum(vec[!is.na(vec)])/sum(!is.na(vec)), 0 )
        
      }))
      
      dupIND <- as.vector(duplication > minproportion)
      
      if(sum(dupIND)){
        
        dups <- sapply((1:(n-1))[dupIND], function(id){c(ID1 = sortIDs[id], ID2 = sortIDs[id+1])}) %>% t()
        
        missing <- t(sapply(1:nrow(dups), function(dup){
          
          vec <- match(dups[dup,], sortIDs);
          
          c("Missing1" = sum(is.na(sort.scores.df[vec[1], ])), "Missing2" = sum(is.na(sort.scores.df[vec[2], ])))
          
        }, simplify = TRUE))
        
        report <- dplyr::bind_cols(dplyr::as_tibble(dups), 
                                  dplyr::as_tibble(missing), 
                                  proportion = duplication[dupIND]
        ) %>% dplyr::mutate(silly = silly) %>% 
          dplyr::select(silly, tidyr::everything())
        
      }
      
      if(!sum(dupIND)){
        
       report <- tibble::tibble(silly = silly, ID1 = NA, ID2 = NA, Missing1 = NA, Missing2 = NA, proportion = NA) 
        
      }
      
      }
      
      report
      
    } #End multicore loop
    
    doParallel::stopCluster(cl)
    
    output <- dupcheck %>% 
      dplyr::bind_rows()
    
    Sys.time()-start.time
    
    return(output)
    
  }#End if NULL quantile
  
  # Start if quantile 
  if(!is.null(quantile)){
    
    cl <- doParallel::makePSOCKcluster(ncores)
    
    doParallel::registerDoParallel(cl, cores = ncores)  
    
    # multicore loop
    dupcheck <- foreach::foreach(silly = sillyvec) %dopar% {
      
      new.gcl <- my.gcl[[silly]]
      
      IDs <- new.gcl$FK_FISH_ID
      
      n <- silly_n.GCL("new")$n
      
      nloci <- length(loci)
      
      #Combine dose 1 and 2 into single columns for each locus separated by a period.
      scores0 <- new.gcl[ ,c("FK_FISH_ID", scores_cols)] 
      
      if(n < 2){
        
        report <- tibble::tibble(silly = silly, ID1 = NA, ID2 = NA, Missing1 = NA, Missing2 = NA, proportion = NA) #If only one individual the results will be a tibble of NAs. Can't do next() using applys.
        
      } else{

      ncombs <- choose(n, 2)        

      mycombs <- combn(IDs, 2)

      scores1 <- sapply(IDs, function(id){
        
        ind <- scores0 %>% 
          filter(FK_FISH_ID==id)
        
        sapply(loci, function(locus){
          
          ifelse(is.na(ind[ ,locus]), NA, paste(ind[ , locus], ind[ , paste0(locus,".1")], sep = "."))
          
        })
        
      }) %>% t() %>% 
        data.frame () 
      
      dimnames(scores1)[[1]] <- IDs
      
      duplication <- sapply(1:ncombs, function(comb){
        
        compair <- scores1[as.character(mycombs[1, comb]), ] == scores1[as.character(mycombs[2, comb]), ];
        
        sum(compair[!is.na(compair)])/sum(!is.na(compair))
        
        }) %>% purrr::set_names(sapply(1:ncombs, function(comb){
          
          paste(mycombs[,comb], collapse=".")
          
          }))
 
      proportion <- max(quantile(duplication, quantile), minproportion)
 
      dupIND <- duplication > proportion
    
      if(sum(dupIND)){

        dups <- data.frame(ID1 = mycombs[1, dupIND], ID2 = mycombs[2, dupIND])
        
        missing <- t(sapply(1:nrow(dups),function(dup){
          
          vec <- match(as_vector(dups[dup, ]), IDs);
          
          c("Missing1" = sum(is.na(scores1[vec[1], ])), "Missing2" = sum(is.na(scores1[vec[2], ])))
          
        }, simplify = TRUE))
        
        report <- dplyr::bind_cols(dplyr::as_tibble(dups), 
                                   dplyr::as_tibble(missing), 
                                   proportion = duplication[dupIND])%>% mutate(silly = silly) %>% 
          select(silly, tidyr::everything())
                                   
      }
      if(!sum(dupIND)){
  
        report <- tibble::tibble(silly = silly, ID1 = NA, ID2 = NA, Missing1 = NA, Missing2 = NA, proportion = NA) 

      }
      
      }
      
      DupDist <- tibble::tibble(silly = silly, IDs = names(duplication), duplication) %>% 
        separate(IDs, into = c("ID1", "ID2"), sep = "[:punct:]")

      list(report = report, DupDist = DupDist)  
       
    }#End multicore loop
    
    parallel::stopCluster(cl)
  
    report <- lapply(1:length(sillyvec), function(i){
      
      dupcheck[[i]][["report"]]
      
      }) %>% bind_rows()
    
    DupDist <- lapply(1:length(sillyvec), function(i){
      
      dupcheck[[i]][["DupDist"]]
      
    }) %>% bind_rows()
    
    output <- list(report = report, DupDist = DupDist)
    
    Sys.time()-start.time
    
    return(output)

  }#end if quantile
    
}