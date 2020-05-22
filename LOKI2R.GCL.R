LOKI2R.GCL <- function(sillyvec, username, password){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function connects to LOKI and creates a "*.gcl" object for each silly in sillyvec.  
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   sillyvec - a character vector of silly codes (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART09"))
  #
  #   username - your state user name
  #
  #   password - your password used to access LOKI - see Eric Lardizabal if you don't have a passord for LOKI
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function assigns a tibble with the following columns for each silly.
  #
  #               Columns 1-19
  #                     FK_FISH_ID <character>: fish ID numbers for each individual  
  #                     COLLECTION_ID <double>: the unique collection ID nubmer for each individual
  #                     SILLY_CODE <character>: the silly code for each individual
  #                     PLATE_ID <character>: the extraction plate ID for each individual
  #                     PK_TISSUE_TYPE <character>: the tissue type extracted for DNA for each individual
  #                     CAPTURE_LOCATION <character>: the location were each individual was captured for sampling
  #                     CAPTURE_DATE <date>: the date each individual was captured (e.g. May 5, 2020 = "2020-05-05") 
  #                     END_CAPTURE_DATE <date>: the last collection date for a silly (e.g. May 5, 2020 = "2020-05-05") 
  #                     MESH_SIZE <character>: the mesh size of the net used to capture (harvest) each individual 
  #                     MESH_SIZE_COMMENT <character>: comments about mesh size
  #                     LATITUDE <double>: the latitude where each individual was captured in decimal degrees
  #                     LONGITUDE <double>: the longitude where each individual was captured in decimal degrees
  #                     AGENCY <character>: the name of the agency or organization that collected each individual
  #                     VIAL_BARCODE <character>: the barcode on the collection vial
  #                     DNA_TRAY_CODE <character>: the barcode on the collection tray/card
  #                     DNA_TRAY_WELL_CODE <double>: the unique number assigned to each postion in the collection tray/card for each individual(e.g postions A1-A10 = codes 1-10, )
  #                     DNA_TRAY_WELL_POS <double>: the postion in the collection tray/card (e.g. A1, A2, B1, B2, etc.)
  #                     CONTAINER_ARRAY_TYPE_ID <doupble>: the number code for the collection container (e.g. tray or card)    
  #                     SillySource <doupble>: the original silly code and fish ID for each individual (e.g. KQUART06_1). When pulled from loki this will be the SILLY_CODE and FK_FISH_ID
  #                               
  #               Columns 20+
  #                     The remaining columns in the object will be the scores for all loci in the LocusControl object. 
  #                     Each locus will have a column for each dose. The columns will be named after the locus with a number added to the locus name after dose 1 (e.g. dose 1 = GTH2B-550; dose 2 = GTH2B-550.1) 
  #
  #
  #   The tibbles will be named after the silly code with a .gcl extention (e.g. KQUART06.gcl)
  #
  #   If one or more *.gcl objects have no data for one or more loci, the function returns a tibble of loci with missing data for each *.gcl object.
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  #   password = "************"
  #
  #   CreateLocusControl.GCL(markersuite = "UCI_Chinook_GTSeq_557SNPs", username = "awbarclay", password = password)
  # 
  #   LOKI2R.GCL(sillyvec = c("KCURRY13", "KDEEP10", "CHDEC912", "KCHICK10", "KALEX16"), username = "awbarclay", password = password)
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function requires a LocusControl object. Run CreateLocusControl.GCL prior to this function.
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  if(!exists("LocusControl")){
    
    stop("'LocusControl' not yet built.")
    
  }
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(RJDBC, tidyverse, lubridate) #Install packages, if not in library and then load them.
  
  # This copies the "odbc8.jar" file to the R folder on your computer if it doesn't exist there. This file contains the java odbc drivers needed for RJDBC
  
  if(!file.exists(path.expand("~/R"))){
    
    dir <- path.expand("~/R")
    
    dir.create(dir)
    
    bool <- file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))
    
  } else {
    
    if(!file.exists(path.expand("~/R/ojdbc8.jar"))){
      
      bool <- file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))
      
    }
    
  }
  
  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {
    
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = "C:/Program Files/R/RequiredLibraries/ojdbc8.jar", " ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
    
  } else {
    
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = path.expand("~/R/ojdbc8.jar"), " ")
    
  }
  
  url <- LOKI_URL.GCL() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password) #The database connection
  
  loci <- LocusControl$locusnames
  
  nloci <- length(loci)
  
  ploidy <- LocusControl$ploidy
  
  alleles <- LocusControl$alleles
  
  nalleles <- LocusControl$nalleles 
  
  gnoqry <- paste("SELECT * FROM AKFINADM.V_GNOQRY WHERE LOCUS IN (", paste0("'", loci, "'", collapse = ","), ") AND SILLY_CODE IN (", paste0("'", sillyvec, "'", collapse = ","), ")", sep = "") #Gentoype query
  
  dataAll0 <- RJDBC::dbGetQuery(con, gnoqry)  #Pulling data from LOKI using the connection and genotype query
  
  dataAllbool <- dataAll0$PK_TISSUE_TYPE == dataAll0$PREFERRED_TISSUE 
  
  bothNAbool <- is.na(dataAll0$PREFERRED_TISSUE) &  is.na(dataAll0$PK_TISSUE_TYPE)
  
  dataAllbool[is.na(dataAllbool)] <- FALSE 
  
  dataAllbool[bothNAbool] <- TRUE  
  
  dataAll <- dataAll0[dataAllbool, ] %>% 
    tibble::as_tibble() %>% 
    dplyr::select(-PREFERRED_TISSUE, -TEST_TYPE)
  
  discon <- RJDBC::dbDisconnect(con)
  
  # what sillys have no data for any of these loci?
  missing_sillys <- setdiff(sillyvec, dataAll$SILLY_CODE %>% unique()) #Find which sillys had no data for any loci in LocusControl
  
  # what indvs are missing loci from LocusControl (i.e. no genotyping attempted)?
  missing_indvs <- dataAll %>% 
    dplyr::count(SILLY_CODE, FISH_ID) %>% 
    dplyr::filter(n < nloci)
  
  # what loci are these indvs missing?
  missing_indvs_loci <- dataAll %>% 
    dplyr::distinct(SILLY_CODE, FISH_ID) %>% 
    tibble::add_column(!!!purrr::set_names(x = rep(NA_real_, nloci), nm = loci)) %>% 
    tidyr::gather(LOCUS, na, -SILLY_CODE, -FISH_ID) %>% 
    dplyr::select(-na) %>% 
    dplyr::anti_join(dplyr::select(.data = dataAll, SILLY_CODE, FISH_ID, LOCUS), by = c("SILLY_CODE", "FISH_ID", "LOCUS")) %>% 
    tidyr::nest(missing_loci = LOCUS) %>% 
    dplyr::right_join(missing_indvs, by = c("SILLY_CODE", "FISH_ID")) %>% 
    dplyr::rename(n_loci_genotyped = n)
  
  # filter out individuals missing loci, replace no calls (0's) with NA
  dataAll <- dataAll %>% 
    dplyr::anti_join(dplyr::select(.data = missing_indvs, SILLY_CODE, FISH_ID), by = c("SILLY_CODE", "FISH_ID")) %>% 
    dplyr::mutate(ALLELE_1 = dplyr::na_if(ALLELE_1, "0"),
                  ALLELE_2 = dplyr::na_if(ALLELE_2, "0"))
  
  # what sillys have complete data?
  dataAll_sillys <- dataAll$SILLY_CODE %>% 
    unique() %>% 
    sort()
  
  # did all indvs from a silly get dropped?
  missing_sillys_indv <- setdiff(setdiff(sillyvec, missing_sillys), dataAll_sillys)
  
  lapply(dataAll_sillys, function(silly){ 
    
    message0 <- paste0(silly, ".gcl created ", match(silly, dataAll_sillys)," of ", length(dataAll_sillys)," completed.") 
    
    sillydata <- dataAll %>% 
      dplyr::filter(SILLY_CODE==silly)
    
    ids <- sillydata$FISH_ID %>% 
      unique() %>% 
      sort() %>% 
      as.character()
    
    sillyvials <- paste(silly, ids, sep = "_")
    
    nind <- length(sillyvials)
    
    silly_df_cols <- rep(NA_real_, nloci*2) %>% 
      purrr::set_names(c(loci, paste0(loci, ".1")) %>% sort()) 
    
    silly_df0 <- sillydata %>%
      dplyr::arrange(LOCUS) %>%
      tidyr::pivot_longer(cols = c("ALLELE_1", "ALLELE_2"), values_to = "Allele") %>% 
      dplyr::mutate(scores_header = case_when(name == "ALLELE_2" ~ paste0(LOCUS, ".1"), 
                                              TRUE ~ LOCUS)) %>% 
      dplyr::select(-LOCUS, -name) %>% 
      tidyr::pivot_wider(names_from = scores_header, values_from = Allele, names_sep="" ) %>% 
      dplyr::mutate(
        CAPTURE_DATE = lubridate::as_date(CAPTURE_DATE),
        END_CAPTURE_DATE = lubridate::as_date(END_CAPTURE_DATE),
        SillySource = paste(SILLY_CODE, FISH_ID, sep = "_")
      )
    
    silly_df <- tibble::add_column(silly_df0, !!!silly_df_cols[setdiff(names(silly_df_cols), names(silly_df0))]) %>%
      dplyr::select(
        FK_FISH_ID = FISH_ID,
        COLLECTION_ID,
        SILLY_CODE,
        PLATE_ID,
        PK_TISSUE_TYPE,
        CAPTURE_LOCATION,
        CAPTURE_DATE,
        END_CAPTURE_DATE,
        MESH_SIZE,
        MESH_SIZE_COMMENT,
        LATITUDE,
        LONGITUDE,
        AGENCY,
        VIAL_BARCODE,
        DNA_TRAY_CODE,
        DNA_TRAY_WELL_CODE,
        DNA_TRAY_WELL_POS,
        CONTAINER_ARRAY_TYPE_ID,
        SillySource,
        tidyselect::all_of(names(silly_df_cols))
      ) %>%
      dplyr::arrange(FK_FISH_ID)
    
    message(message0)
    
    assign(paste0(silly, ".gcl"), silly_df, pos = 1, .GlobalEnv)
    
  })
  
  if(length(missing_sillys) >= 1){
    
    warning(paste0("The following sillys had no data in LOKI for any of the loci in LocusControl:\n", paste0(missing_sillys, collapse = "\n")), call. = FALSE)
    
  }
  
  if(length(missing_sillys_indv) >= 1){
    
    warning(paste0("The following sillys had no individuals with complete data in LOKI for the loci in LocusControl:\n", paste0(missing_sillys_indv, collapse = "\n")), call. = FALSE)
    
  }
  
  if(nrow(missing_indvs_loci) >= 1){
    
    n_missing <- missing_indvs_loci %>% 
      dplyr::mutate(n_loci_missing = nloci - n_loci_genotyped) %>% 
      dplyr::count(SILLY_CODE, n_loci_missing) %>%
      dplyr::rename(n_indv = n) %>% 
      dplyr::mutate(silly_n_miss = paste0(SILLY_CODE, " (", n_indv, " individuals missing ", n_loci_missing, " loci)")) %>% 
      dplyr::pull(silly_n_miss)
    
    warning(paste0("The following sillys had individuals that were missing data for one or more loci:\n", paste(n_missing, collapse = "\n")), call. = FALSE)
    
    warning(paste0("A table of loci missing data for each individual has been assigned to the object 'missing_indvs_loci'"), call. = FALSE)
    
    assign(x = "missing_indvs_loci", value = missing_indvs_loci, pos = 1, .GlobalEnv)
    
  } else { 
    
    print("The *.gcl objects created have data for all loci in LocusControl")
    
  }
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime)
  
}