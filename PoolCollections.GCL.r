PoolCollections.GCL <- function(collections, loci = LocusControl$locusnames, IDs = NULL, newname = paste(collections,collapse=".")){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function combines "*.gcl" objects into a new one called "newname.gcl".
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   collections - a character vector of silly codes without the ".gcl" extention (e.g. collections <- c("KQUART06","KQUART08","KQUART10")). 
  #                 Collections can be a single silly if you want create a new ".gcl" with only fish supplied in IDs.
  #
  #   loci - a character vector of locus names
  # 
  #   IDs - a named list of fish ID vectors (either character or numeric), each vector is associated with and named after a member of "collections".
  #         These will be used to subset each collection before pooling. If no IDs are supplied all individuals from each collection are used.
  #
  #   newname - is the name of the new "*.gcl" created. Do not provide ".gcl" extention. If no name supplied then the newname defaults to
  #             the collection names collapsed with a period between each name (e.g. "KQUART06.KQUART08.KQUART09")
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Assigns a new "pooled collection" to your workspace
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  # 
  #   removedInd <- (collections = c("KQUART06","KQUART08","KQUART10"), loci = loci, IDs = list(KQUART06 = 3:12, KQUART08 = 1:10, KQUART10 = 1:4), newname = "QuartzCr")
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function is also useful for producing "pooled mixture" objects for mixed stock analysis. 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!exists("LocusControl")){
    
    stop("'LocusControl' not yet built.")
    
  }
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse) #Install packages, if not in library and then load them.

  if(nchar(newname)>200){
    
    newname <- substr(newname, start = 1, stop = 200)
    
    }

  if(sum(is.na(match(loci,LocusControl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci,LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
    }

  ncollections <- length(collections)

  if(is.null(IDs)){
    
    IDs <- lapply(collections, function(collection){
      
      get(paste0(collection, ".gcl"), pos = 1)$FK_FISH_ID
      
      }) 
    
    }

  if(!is.list(IDs)){
    
    stop("'IDs' must be a list")
    
    }

  if(ncollections!=length(IDs)){
    
    stop("'IDs' must be same length as 'collections'")
    
    }

  IDs <- purrr::set_names(IDs, collections) #Making sure IDs has names
  
  SubsetLoci <- c(loci, paste0(loci, ".1")) %>% sort()#These are the locus score headers for subsetting by loci.
  
  output <- lapply(collections, function(collection){
    
    my.gcl <- get(paste0(collection, ".gcl"), pos = 1)
    
    attr <- my.gcl[ , 1:19] %>% 
      names() #The attribute names
    
    my.gcl %>% 
      dplyr::filter(FK_FISH_ID%in%IDs[[collection]]) %>% 
      dplyr::select(attr, SubsetLoci)
    
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(FK_FISH_ID = seq(length(unlist(IDs))))
    
  assign(paste(newname,".gcl",sep=""), output, pos=1, envir = .GlobalEnv)  

  return(NULL)
  
}
