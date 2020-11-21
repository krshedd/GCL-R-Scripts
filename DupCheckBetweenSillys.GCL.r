DupCheckBetweenSillys.GCL <- function(KeySillys, KeySillyIDs = NULL, BetweenSillys, loci, minnonmissing = 0.6, minproportion = 0.9, ncores = 4){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function check for duplicate individuals between silly ".gcl" objects. This function is mainly used for QC purposes. 
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   KeySillys - a vector of silly codes without the ".gcl" extention (e.g. KeySillys <-  c("KQUART06","KQUART08","KQUART10")). 
  #
  #   KeySillyIDs - a named list of character vector FK_FISH_IDs for each KeySilly to check against all sillys in BetweenSillys.(e.g. KeySillys <-  list(KQUART06 = c("25", "30") , KQUART08 =  ,KQUART10 =))
  #                 If NULL (default), all FK_FISH_IDs for each silly in KeySillys are checked against all sillys in BetweenSillys
  #
  #   BetweenSillys - a vector of silly codes without the ".gcl" extention (e.g. BetweenSillys <-  c("KQUART06","KQUART08","KQUART10"))
  #         
  #   loci - vector of locus names; if set to NULL all loci in the ".gcl" objects will be used.
  #
  #   minnonmissing - the proportion of loci that a pair must share non missing in order to be reported (passed on to rubias function)
  #
  #   minproportion - the proportion of shared non-missing loci that must be shared between the indivdiuals to be reported as a matching pair, this is passed through to rubias as well 
  #
  #   ncores - the number of cores to use in a foreach %dopar% loop. If the number of core exceeds the number on your device, then ncores defaults to detectCores()
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  #   Returns a tibble of with the following variables:  
  #       Keysillyvial - e.g., KeySilly_ID 
  #       Betweensillyvial - e.g., BetweenSilly_ID,
  #       Keymissing - the number of loci without genotypes for each KeySilly_ID
  #       Betweenmissing - the number of loci without genotypes for each BetweenSilly_ID, 
  #       DuplicateRate - the proportion of duplicate genotypes between each KeySilly_ID and BetweenSilly_ID.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # password <- "************"

    CreateLocusControl.GCL(markersuite = "CookInletChinook2013_43SNPs", username = "awbarclay", password = password)
    LOKI2R.GCL(sillyvec =  c("KKILL05","KKILL06", "KFUNN05", "KFUNN06"), username = "awbarclay", password = password)
    
    #Creating new sillys with duplicates for an example
    PoolCollections.GCL(collections = c("KKILL05", "KFUNN05"), IDs = list(KKILL05 = KKILL05.gcl$FK_FISH_ID, KFUNN05 = c(9, 30)), newname = "KKILL05")
    PoolCollections.GCL(collections = c("KFUNN06", "KKILL06"), IDs = list(KFUNN06 = KFUNN06.gcl$FK_FISH_ID, KKILL06 = c(101, 176)), newname = "KFUNN06")
    KKILL05.gcl <- KKILL05.gcl %>% mutate(SillySource = case_when(FK_FISH_ID %in% c(69, 70) ~ paste(SILLY_CODE, FK_FISH_ID, sep = "_"), TRUE~SillySource))
    KFUNN06.gcl <- KFUNN06.gcl %>% mutate(SillySource = case_when(FK_FISH_ID %in% c(184, 185) ~ paste(SILLY_CODE, FK_FISH_ID, sep = "_"), TRUE~SillySource))
    
    KeySillys <- c("KKILL06", "KFUNN05")
    KeySillyIDs <- list(KKILL06 = c(101, 176), KFUNN05 = c(9, 30))
    BetweenSillys <- c("KKILL05", "KFUNN06")
    loci <- LocusControl$locusnames
    results <- DupCheckBetweenSillys.GCL(KeySillys = c("KKILL06", "KFUNN05"), KeySillyIDs = list(KKILL06 = c(101, 176), KFUNN05 = c(9, 30)), BetweenSillys = c("KKILL05", "KFUNN06"), loci = LocusControl$locusnames, minnonmissing = 0.6, minproportion = 0.9, ncores = 4)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  start.time <- Sys.time() 
  
  if(!all(loci %in% LocusControl$locusnames)){
    
    stop(paste0("The following `loci` were not found in `LocusControl`:\n", paste(setdiff(loci, LocusControl$locusnames), collapse = "\n")))
    
  }
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse, doParallel, parallel, foreach, rubias)  # Install packages, if not in library and then load them.
  
  if(ncores > parallel::detectCores()) {
    
    ncores <- parallel::detectCores()
    
  }
  
  #while(!require(abind)){install.packages("abind")}

  #while(!require(outliers)){install.packages("outliers")}

  ploidy <- LocusControl$ploidy[loci]
  
  nloci <- length(loci)
  
  scores_cols <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
    as.vector()  # This keeps the scores columns in the correct order when there are loci with similar names.
  
  attr <-
    c(
      "FK_FISH_ID",
      "COLLECTION_ID",
      "SILLY_CODE",
      "PLATE_ID",
      "PK_TISSUE_TYPE",
      "CAPTURE_LOCATION",
      "CAPTURE_DATE",
      "END_CAPTURE_DATE",
      "MESH_SIZE",
      "MESH_SIZE_COMMENT",
      "LATITUDE",
      "LONGITUDE",
      "AGENCY",
      "VIAL_BARCODE",
      "DNA_TRAY_CODE",
      "DNA_TRAY_WELL_CODE",
      "DNA_TRAY_WELL_POS",
      "CONTAINER_ARRAY_TYPE_ID",
      "SillySource"
    )

  # my.key
  if(!is.null(KeySillyIDs)){
    
    my.key <- lapply(KeySillys, function(silly){
      
      get(paste0(silly, ".gcl")) %>% 
        filter(FK_FISH_ID %in% KeySillyIDs[[silly]])
      
    }) %>% bind_rows()
      
  } else{
    
    
    my.key <- lapply(KeySillys, function(silly){
      
      get(paste0(silly, ".gcl"))
      
    }) %>% bind_rows()
    
  }
  
  # my.between
  my.between <- lapply(BetweenSillys, function(silly){
    
    get(paste0(silly, ".gcl"))
    
  }) %>% set_names(BetweenSillys)
  
  # Added this if statement code for haploid markers, rubias::close_matching_samples() was counting them as missing loci because of the NAs in the allele2 column. 
  # This can be removed Eric Anderson fixes the function. 
  if(any(ploidy == 1)) {
    
    haploci <- names(ploidy[ploidy == 1])
    
    for(locus in haploci){
      
      my.key[[paste0(locus, ".1")]] <- my.key[[locus]]
      
      my.between <- lapply(BetweenSillys, function(silly){
        
         my.between[[silly]][[paste0(locus, ".1")]] <- my.between[[silly]][[locus]]
        
      })
      
    }
  } 
    
  ## Loop through between sillys
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  dupcheck0 <- foreach::foreach(silly = BetweenSillys, .export = c("loci"), .packages = c("tidyverse","rubias")) %dopar% {
    
    new.gcl <- my.between[[silly]] %>%
      bind_rows(my.key) %>% 
      dplyr::mutate(
        sample_type = "reference",
        repunit = NA_character_,
        collection = SILLY_CODE
      ) %>%
      tidyr::unite(col = "indiv", SILLY_CODE, FK_FISH_ID, sep = "_") %>% 
      dplyr::select(sample_type,
                    repunit,
                    collection,
                    indiv,
                    tidyselect::all_of(scores_cols))
    
    dups <- rubias::close_matching_samples(D = new.gcl, gen_start_col = 5, min_frac_non_miss = minnonmissing, min_frac_matching = minproportion) %>% 
      filter(collection_1 == silly, collection_2 %in% KeySillys)
    
    dups %>% 
      filter(collection_1 == silly, collection_2 %in% KeySillys) %>% 
      rename(Betweensillyvial = indiv_1, Keysillyvial = indiv_2, Betweensilly = collection_1, Keysilly = collection_2) %>% 
      mutate(DuplicateRate = num_match/num_non_miss) %>% 
      select(Keysillyvial, Betweensillyvial, DuplicateRate)
    
  } %>% dplyr::bind_rows()  # End multicore loop
  
  parallel::stopCluster(cl) 
  
  if(dim(dupcheck0)[1]>=1){
    
    lapply(dupcheck0 %>% pull(Keysillyvial), function(ind){
      
      
      keygcls %>% 
        filter(SillySource == KeySillyIDs)
      
    })
    
    
  }
  
  

  newKeySillys=paste0("new",KeySillys)

  newname=paste0(newKeySillys,collapse=".")

  EMPTY=PoolCollections.GCL(newKeySillys,loci,IDs=NULL,newname=newname)

  newkeygcl=get(paste0(newname,".gcl")) 

  newscores=newkeygcl$scores[,loci,,drop=FALSE]

  newsillyvials=newkeygcl$attributes$SillySource

  rownames(newscores)=newsillyvials

  betweensillys=paste0("between",BetweenSillys)
   
  betweenname=paste(betweensillys,collapse=".")

  EMPTY=PoolCollections.GCL(BetweenSillys,loci,IDs=NULL,newname="betweenname")

  betweengcl=get(paste0("betweenname",".gcl"))

  betweenscores=betweengcl$scores[,loci,,drop=FALSE]

  betweensillyvials=betweengcl$attributes$SillySource

  rownames(betweenscores)=betweensillyvials

  combs=expand.grid(Keysillyvial=newsillyvials,Betweensillyvial=betweensillyvials,stringsAsFactors=FALSE)

  scores=abind(newscores,betweenscores,along=1)

  scores[scores%in%c("0","Unk","XXX")]=NA

  allsillyvials=rownames(scores)

  scores.df=data.frame(t(sapply(allsillyvials,function(id){sapply(loci,function(locus){ifelse(is.na(scores[id,locus,1]),NA,paste(scores[id,locus,1:ploidy[locus]],collapse="."))})})),stringsAsFactors=FALSE)

  duplicaterate=rowMeans(scores.df[combs[,1],]==scores.df[combs[,2],],na.rm=TRUE)

  
  
  duplicatesummary <- sapply(unique(combs$Keysillyvial), function(KeySillyID) {
    
    combs.log <- combs$Keysillyvial == KeySillyID
    
    hist(duplicaterate[combs.log], breaks = seq(0, 1, 0.01), col = 8, xlab = "Duplicate Rate", main = KeySillyID)
    
    abline(v = 0.9, col = "red", lwd = 2)
    
    maxduplicaterate.log <- duplicaterate[combs.log] == max(duplicaterate[combs.log])
    
    maxduplicaterate <- setNames(object = duplicaterate[combs.log][maxduplicaterate.log], nm = combs[combs.log, "Betweensillyvial"][maxduplicaterate.log])
    
    text(x = maxduplicaterate[1], y = 5, labels = paste(names(maxduplicaterate), collapse = "_"), srt = 90, adj = 0)
    
    keymissing <- apply(is.na(scores[combs[combs.log, "Keysillyvial"][maxduplicaterate.log], , 1, drop = FALSE]), 1, sum)
    
    betweenmissing <- apply(is.na(scores[combs[combs.log, "Betweensillyvial"][maxduplicaterate.log], , 1, drop = FALSE]), 1, sum)
    
    maxdups <- cbind(combs[combs.log, ][maxduplicaterate.log, ], Keymissing = keymissing, Betweenmissing = betweenmissing, DuplicateRate = maxduplicaterate, stringsAsFactors=FALSE,row.names=seq(sum(maxduplicaterate.log)))
    
    
    
    thresholdduplicaterate.log <- duplicaterate[combs.log] > threshold
    
    thresholdduplicaterate <- setNames(object = duplicaterate[combs.log][thresholdduplicaterate.log], nm = combs[combs.log, "Betweensillyvial"][thresholdduplicaterate.log])
    
    keymissing <- apply(is.na(scores[combs[combs.log, "Keysillyvial"][thresholdduplicaterate.log], , 1, drop = FALSE]), 1, sum)
    
    betweenmissing <- apply(is.na(scores[combs[combs.log, "Betweensillyvial"][thresholdduplicaterate.log], , 1, drop = FALSE]), 1, sum)
    
    if(sum(thresholdduplicaterate.log) > 0) {
    
      thresholddups <- cbind(combs[combs.log, ][thresholdduplicaterate.log, ], Keymissing = keymissing, Betweenmissing = betweenmissing, DuplicateRate = thresholdduplicaterate, stringsAsFactors=FALSE,row.names=seq(sum(thresholdduplicaterate.log)))
      
    } else{
      
      thresholddups <- NULL #  paste("No Duplicates at a threshold of", threshold)
      
    }
    
    
    
    if(grep(pattern = "QC", x = combs[1, ]) == 1) {
      projectduplicaterate.log <- combs$Keysillyvial == KeySillyID & combs$Betweensillyvial == gsub(pattern = "QC", replacement = "", x = KeySillyID)
    } else {
      projectduplicaterate.log <- combs$Keysillyvial == KeySillyID & combs$Betweensillyvial == gsub(pattern = "_", replacement = "QC_", x = KeySillyID)
    }
    
    projectduplicaterate <- setNames(object = duplicaterate[projectduplicaterate.log], nm = combs[projectduplicaterate.log, "Betweensillyvial"])

    keymissing <- apply(is.na(scores[combs[projectduplicaterate.log, "Keysillyvial"], , 1, drop = FALSE]), 1, sum)
    
    betweenmissing <- apply(is.na(scores[combs[projectduplicaterate.log, "Betweensillyvial"], , 1, drop = FALSE]), 1, sum)
    
    projectdups <- cbind(combs[projectduplicaterate.log, ], Keymissing = keymissing, Betweenmissing = betweenmissing, DuplicateRate = projectduplicaterate, stringsAsFactors=FALSE,row.names=seq(sum(projectduplicaterate.log)))
    
    
    
    # list(Threshold = thresholddups, MostSimilar = maxdups, Project = projectdups)
    dplyr::bind_rows(thresholddups, maxdups, projectdups) %>% 
      dplyr::distinct(Keysillyvial, Betweensillyvial, .keep_all = TRUE)
  
  }, simplify = FALSE)
  
  duplicatesummary <- dplyr::bind_rows(duplicatesummary)
  
  return(duplicatesummary)
    
}