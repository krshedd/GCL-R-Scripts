GenepopTOgcl.GCL <- function(filename, CreateLocusControl = TRUE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function creates ".gcl" and LocusControl objects from a GENEPOP file. 
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   filename - the full file path of the GENEPOP file.
  #
  #   CreateLocusControl - logical; if TRUE (default) a LocusControl object will be created
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #  This function assigns the ".gcl" and LocusControl objects to the current workspace
  # 
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # GenepopTOgcl.GCL(filename = "example/CI67pops413loci.gen", CreateLocusControl = TRUE)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(strsplit(filename, ".gen")==filename){
    
    stop("Genepop file needs a '.gen' extention, hoser!!!")
    
    }

  if(!CreateLocusControl&!exists("LocusControl")){
    
    stop("LocusControl must exist if you don't create one!")
    
  }

  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(adegenet, tidyverse)  # Install packages, if not in library and then load them.
  
  rawdat <- scan(filename, what = "", sep = "\n")

  len <- length(rawdat)

  popind <- sapply(rawdat, function(lin){
    
    strsplit(lin, " ")[[1]][1]=="Pop"
    
    }) %>% as.vector()

  npops <- sum(popind)
 
  ORD <- order(popind, decreasing = TRUE)[1:npops]

  loci <- sapply(rawdat[2:(ORD[1]-1)], function(str){
    
    strsplit(str, " ")[[1]][1]
    
    })

  if(!CreateLocusControl){
    
    if(sum(is.na(match(loci, LocusControl$locusnames)))){
      
      stop("Loci in Genepop file are absent from LocusControl, hoser!!!")
      
    }
    
    loci <- LocusControl$locusnames
  }

  nloc <- length(loci)
  
  # Attribute variables
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
  
  # locus variables
  loc_vars <- c(loci, paste0(loci, ".1")) %>% 
    as.character() %>% 
    sort()
  
  # Put data in .gcl object tibble format
  # create initial data tibble
  dat0 <- rawdat[-(1:(nloc+1))][rawdat[-(1:(nloc+1))]!="Pop"] %>%
    as_tibble() %>% 
    separate(col = value, sep = " ,  ", into = c("SillySource", "geno")) %>% 
    separate(col = geno, sep = " ", into = loci) %>% 
    separate(col = SillySource, sep = "_", into = c("SILLY_CODE", "FK_FISH_ID"), remove = FALSE) %>% 
    mutate(COLLECTION_ID = NA, 
           PLATE_ID = NA, 
           PK_TISSUE_TYPE = NA, 
           CAPTURE_LOCATION = NA, 
           CAPTURE_DATE = NA, 
           END_CAPTURE_DATE = NA, 
           MESH_SIZE = NA, 
           MESH_SIZE_COMMENT = NA,
           LATITUDE = NA,
           LONGITUDE = NA,
           AGENCY = NA,
           VIAL_BARCODE = NA,
           DNA_TRAY_CODE = NA,
           DNA_TRAY_WELL_CODE = NA,
           DNA_TRAY_WELL_POS = NA,
           CONTAINER_ARRAY_TYPE_ID = NA
           )
  
  # Separate genotypes into two variables
  nchar <- dat0 %>% 
    select(all_of(loci)) %>% 
    summarize(across(everything(), nchar)) %>% 
    max()
  
  for(locus in loci){
    
    dat0 <- dat0 %>%
      separate(!!sym(locus), into = c(locus, paste0(locus, ".1")), sep = nchar/2, remove = FALSE) %>% 
      mutate(across(all_of(c(locus, paste0(locus, ".1"))), .fns = as.numeric))
    
  }
  
  # Final data tibble
  dat <- dat0 %>% 
    mutate(across(all_of(loc_vars), ~na_if(x = ., y=0)))
  
  # Assign silly objects to workspace
  sillyvec <- dat$SILLY_CODE %>% unique()
  
  for(silly in sillyvec){
    
    my.dat <-dat %>% filter(SILLY_CODE == silly)
    
    assign(x = paste0(silly, ".gcl"), value = my.dat, pos = 1)
    
  }
  
  message(paste0("A total of ", length(sillyvec), " '.gcl' objects were created from the GENEPOP file and assgined to the current workspace."))
  
  # Make LocusControl if none exists

  if(!exists("LocusControl") & CreateLocusControl){
    
    # alleles
    alleles <- lapply(loci, function(locus){
      
     tibble(allele = dat %>% 
              select(starts_with(locus)) %>% 
              unlist() %>% 
              unique() %>% 
              na.omit() %>% 
              as.numeric() %>% 
              sort() %>% 
              as.integer(), 
            call = NA_character_)
    }) %>% set_names(loci) 
    
    # nalleles
    nalleles <- lapply(loci, function(locus){
      
      dim(alleles[[locus]])[[1]]
      
    }) %>% 
      unlist %>% 
      set_names(loci)
    
    # ploidy
    ploidy <-  if_else (
      condition = dat0[ ,paste0(loci, ".1")] %>%is.na() %>% apply(., 2, sum)  == 0, 
      true = 2, 
      false = 1) %>% set_names(loci)
    
    assign("LocusControl", tibble::tibble(MarkerSuite = "From GENEPOP", 
                                          locusnames = loci,
                                          Publishedlocusnames = NA_character_, 
                                          nalleles = nalleles, 
                                          ploidy = ploidy,
                                          alleles = alleles), pos = 1)
    
    message("")
    }# End LocusControl
  
}