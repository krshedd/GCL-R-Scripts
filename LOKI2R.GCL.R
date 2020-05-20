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
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  #
  #   password = "************"
  # 
  #   LOKI2R.GCL(sillyvec = sillyvec157, username = "awbarclay", password = password)
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
    
    bool <- file.copy(from="V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar",to=path.expand("~/R/ojdbc8.jar"))
    
  } else {
    
    if(!file.exists(path.expand("~/R/ojdbc8.jar"))){
      
      bool <- file.copy(from="V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar",to=path.expand("~/R/ojdbc8.jar"))
      
    }
    
  }
  
  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {
    
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = "C:/Program Files/R/RequiredLibraries/ojdbc8.jar", " ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
    
  } else {
    
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath=path.expand("~/R/ojdbc8.jar"), " ")
    
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
  
  #Create tibble data fame for each silly in sillyvec
  for(silly in sillyvec){ 
    
    message0 <- paste0(silly, ".gcl created ", match(silly, sillyvec)," of ", length(sillyvec)," completed.") 
    
    sillydata <- dataAll %>% 
      dplyr::filter(SILLY_CODE==silly)
    
    ids <- sillydata$FISH_ID %>% 
      unique() %>% 
      sort() %>% 
      as.character()
    
    sillyvials <- paste(silly, ids, sep = "_")
    
    nind <- length(sillyvials)
    
    if(length(ids) == 0){
      
      message0 <- paste0(silly," is empty.")
      
      message(message0)
      
      next
      
    }
    
    silly_df <- sillydata %>%
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
      ) %>%
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
        everything()
      ) %>%
      dplyr::arrange(FK_FISH_ID)
    
    message(message0)
    
    assign(paste0(silly, ".gcl"), silly_df, pos = 1, .GlobalEnv)
    
  }  # silly 
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime) 
  
  return(fulltime)
  
}