HoFisFstTable.GCL <- function(sillyvec, loci, fstatdir = NULL, dir = NULL, ncores = 4){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function creates a table of statistics containing obseved hetrozygosity (Ho), Fis, and Fst for each locus in loci.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   loci - a character vector of locus names
  #   
  #   fstatdir - the path to an existing FSTAT .dat file. If none is supplied the user must supply a directory (dir) where 
  #              an FSTAT .dat file can be written.
  #
  #   dir -  directory where the FSTAT .dat file is dumped.
  #
  #   ncores - a numeric vector of length one indicating the number of cores to use
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   A tibble is returned containing locusname, Ho, Fis, and Fst
  #
  #  if is.null(fstatdir) and exists(dir) "fstatfile.dat" is put into "dir"
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  #
  #   HoFisFstTable.GCL(sillyvec = sillyvec31, loci = loci82, fstatdir =  NULL, dir = "FSTAT", ncores = 8)
  #
  #   HoFisFstTable.GCL(sillyvec = sillyvec31, loci = loci82, fstatdir = "FSTAT/fstat.dat", dir = NULL, ncores = 8)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!exists("LocusControl")) {
    
    stop("'LocusControl' is required and not found, please create.")
    
  }
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse, doParallel, foreach, parallel, hierfstat) #Install packages, if not in library and then load them.
  
  if(ncores > detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  start.time <- Sys.time()
  
  ploidy <- LocusControl$ploidy[loci]
  
  if(is.null(fstatdir)) {
    
    fstatdir <- paste0(dir, "/fstat.dat")
    
    gcl2FSTAT.GCL(sillyvec = sillyvec, loci = loci, path = fstatdir, ncores = ncores)
  
  }

  dat <- read.fstat.data(fstatdir) #Read in FSTAT file

  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores) #Start cluster  
  
  #Get variance components
  #Start parallel loop
  MyVC <- foreach::foreach(locus = loci, .packages = "hierfstat") %dopar% {
    
    diploid <- ploidy[locus]==2
    
    if(diploid){
      
      VC <- hierfstat::varcomp(dat[, c("Pop", locus)], diploid = diploid)$overall
      
      table <- tibble::tibble(locus = !!locus, P = VC[1], I = VC[2], G = VC[3])
      
    }
    
    if(!diploid){
      
      VC <- hierfstat::varcomp(dat[, c("Pop", locus)], diploid = diploid)$overall
      
      table <- tibble::tibble(locus = !!locus, P = VC[1], I = NA, G = VC[2])
      
    } 
    
    table
    
  } %>%  dplyr::bind_rows()
  
  parallel::stopCluster(cl)# Stop cluster
  
  #Summarize variance components
  MyTable0 <- MyVC %>% 
    dplyr::rowwise() %>%  
    dplyr::mutate(total = sum(c(P,I,G))) %>% 
    dplyr::mutate(P = P/total, I = I/total) %>% 
    dplyr::select(locus, P, I, total) %>% 
    dplyr::ungroup()
  
  Overall <- MyVC %>% 
    dplyr::summarize(P = sum(P), I = sum(I), total = sum(MyTable0$total)) %>% 
    dplyr::mutate(locus = "Overall", I = I/total, P = P/total)
  
  MyTable <- dplyr::bind_rows(MyTable0, Overall) 
  
  #Heterozygosities
  Hovec <- apply(hierfstat::basic.stats(dat[, c("Pop", loci[ploidy==2])])$Ho, 1, mean)
  
  Ho <- tibble::tibble(locus = c(loci[ploidy==2], "Overall"), Ho = c(Hovec, mean(Hovec)))
  
  #Join varcomp summary with Ho  
  output <- MyTable %>% 
    dplyr::left_join(Ho, by = "locus") %>% 
    dplyr::mutate(Fis = I, Fst = P) %>% 
    dplyr::select(locus, Ho, Fis, Fst)
  
  print(Sys.time() - start.time)

  return(output)
  
}
