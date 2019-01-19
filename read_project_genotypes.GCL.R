read_project_genotypes.GCL <- function(project_name = NULL, sillyvec = NULL, loci = NULL, username, password) {  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function is intended for use in the QC.R script!!!
  #  This function pulls project genotypes to 
  #    1) create "slim" .gcl objects (only some of attributes table)
  #    2) create `LocusControl`
  #
  #  ** Warning - genotypes (pre October 2016) may not exist in the LOKI lookup table with a project name, 
  #               so genotypes will have to be pulled by sillyvec and loci**  
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  project_name - a character vector of one or more project names as spelled in LOKI example: c("P014","P015","P016")
  #  sillyvec - character vector of SILLYs - example: c("SCIMA18","SCIMA17")
  #  loci - character vector of locus names as they are spelled in LOKI - example: c("One_E2","One_MHC2_251","One_Cytb_17")
  #  username - your user name for accessing LOKI through R example: "awbarclay"
  #  password - your password for accessing LOKI through R example: "password"
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  Returns `ProjectSillys` - a character vector of all sillys in project (sillyvec)
  #  Returns "slim" .gcl objects for each silly (slim = not all attributes table)
  #  Returns LocusControl for all loci used in the project (as well as "loci", "nallales" , "ploidy", "alleles")
  #
  #  Written by Kyle Shedd 1/11/19
  #  Based on "GenotypeReport.GCL" by Andy Barclay 10/5/18
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Recording function start time
  start.time <- Sys.time()
  
  # Checking to make sure the correct combination of arguments is being use. 
  # If wrong, the function will stop and print an error message to the console.
  if(is.null(sillyvec) & is.null(loci) & is.null(project_name) |
     !is.null(sillyvec) & !is.null(loci) & !is.null(project_name) |
     !is.null(sillyvec) & is.null(loci) & !is.null(project_name) |
     is.null(sillyvec) & !is.null(loci) & !is.null(project_name) |
     is.null(sillyvec) & !is.null(loci) & is.null(project_name))
  {
    stop("The user must supply one of the following argument combinations:\n  1) sillyvec (for all loci and individuals for each silly),\n  2) sillyvec and loci (all individuals for supplied locus list), or\n  3) project_name (for all individuals and loci in a given project)")
  }
  
  # Requiring the tidyverse and RJDBC packages. Packages will be installed if the user doesn't have them installed
  while(!require(tidyverse)){install.packages("tidyverse") }
  while(!require(RJDBC)){install.packages("RJDBC")}
  while(!require(tools)){install.packages("tools")}
  
  # Making sure the ojdbc jar file exists on the users C drive. 
  # If no jar file exists it will be copied from the v drive to the appropriate location on the users C drive.
  if(!file.exists(path.expand("~/R"))) {
    dir <- path.expand("~/R")
    dir.create(dir)
    bool <- file.copy(from = "V:/Analysis/R files/Scripts/DEV/jars/ojdbc6.jar", to = path.expand("~/R/ojdbc6.jar"))
  } else {
    if(!file.exists(path.expand("~/R/ojdbc6.jar"))) {
      bool <- file.copy(from = "V:/Analysis/R files/Scripts/DEV/jars/ojdbc6.jar", to = path.expand("~/R/ojdbc6.jar"))
    }
  }
  
  # Setting default java.parameters
  options(java.parameters = "-Xmx10g")
  
  if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc6.jar")) {
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = "C:/Program Files/R/RequiredLibraries/ojdbc6.jar", " ")  # https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
  } else {
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = path.expand("~/R/ojdbc6.jar"), " ")
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Connect to LOKI
  con <- RJDBC::dbConnect(drv, "jdbc:oracle:thin:@(DESCRIPTION=(ADDRESS_LIST=(ADDRESS=(PROTOCOL=TCP)(HOST=db-pcfres.dfg.alaska.local)(PORT=1521)))(CONNECT_DATA=(SID=PCFRES)))", username, password)
  
  #~~~~~~~~~~~~~~~~
  # Get genotypes
  # Creating java query when sillyvec and loci are supplied.  
  if(!is.null(sillyvec) & !is.null(loci)){
    # Pulling haploid and diploid loci separately because haploid only have ALLELE_1.
    gnoqry_d <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1, ALLELE_2 FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE LOCUS IN (", paste0("'", loci, "'", collapse = ","), ") AND SILLY_CODE IN", "(", paste0("'", sillyvec, "'", collapse = ","), ") AND PLOIDY = 'D'")
    gnoqry_h <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1 FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE LOCUS IN (", paste0("'", loci, "'", collapse = ","), ") AND SILLY_CODE IN", "(", paste0("'", sillyvec, "'", collapse = ","), ") AND PLOIDY = 'H'")
  } 
  
  # Creating java query when only sillyvec is supplied
  if(!is.null(sillyvec) & is.null(loci)){
    # Pulling haploid and diploid loci separately because haploid markers only have ALLELE_1.
    gnoqry_d <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1, ALLELE_2 FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE SILLY_CODE IN", "(", paste0("'", sillyvec, "'", collapse = ","), ") AND PLOIDY = 'D'")
    gnoqry_h <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1 FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE SILLY_CODE IN", "(", paste0("'", sillyvec, "'", collapse = ","), ") AND PLOIDY = 'H'")
  }
  
  # Creating java query when only project_name is supplied.
  if(!is.null(project_name)){
    # Pulling haploid and diploid loci separately because haploid only have ALLELE_1.
    gnoqry_d <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1, ALLELE_2 FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE LAB_PROJECT_NAME IN", "(", paste0("'", project_name, "'", collapse = ","), ") AND PLOIDY = 'D'")
    gnoqry_h <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1 FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE LAB_PROJECT_NAME IN", "(", paste0("'", project_name, "'", collapse = ","), ") AND PLOIDY = 'H'")
  }
  
  # Pull diploid data and concatenate alleles into one column with "/" separator
  data_d <- RJDBC::dbGetQuery(con, gnoqry_d) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(ALLELES = paste(pmin(ALLELE_1, ALLELE_2), pmax(ALLELE_1, ALLELE_2), sep = "/")) %>%
    dplyr::select(-ALLELE_1, -ALLELE_2) %>% 
    tidyr::separate(ALLELES, into = c("ALLELE_1", "ALLELE_2"), sep = "/", remove = FALSE)
    
  
  # Pull haploid data
  data_h <- RJDBC::dbGetQuery(con, gnoqry_h) %>% 
    dplyr::as_tibble()%>% 
    dplyr::mutate(ALLELES = ALLELE_1)
  
  # Combine diploid and haploid marker into one data frame and put into crosstab format (i.e., unstacked table)
  dataAll <- dplyr::bind_rows(data_d, data_h)
  
  # Get list of unique sillys and assign `ProjectSillys` this is needed for QC script
  sillyvec <- unique(dataAll$SILLY_CODE)
  assign(x = "ProjectSillys", value = sillyvec, pos = 1)
  
  #~~~~~~~~~~~~~~~~
  # Get PLATE_ID from the extraction table
  extr_qry <- paste0("SELECT FK_PLATE_ID, TISSUETYPE, SILLY_CODE, FISH_NO FROM AKFINADM.GEN_DNA_WELL WHERE SILLY_CODE IN", "(", paste0("'", sillyvec, "'", collapse = ","), ")")
  
  data_ex <- RJDBC::dbGetQuery(con,extr_qry) %>% 
    dplyr::as_tibble() 
  
  # Disconnect from LOKI
  discon <- RJDBC::dbDisconnect(con)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get list of unique loci and assign `LocusControl` this is needed for QC script
  loci <- sort(unique(dataAll$LOCUS))
  nloci <- length(loci)
  
  CreateLocusControl.GCL(locusnames = loci, username = username, password = password)
  
  assign(x = "loci", value = LocusControl$locusnames, pos = 1)
  assign(x = "nalleles", value = LocusControl$nalleles, pos = 1)
  assign(x = "ploidy", value = LocusControl$ploidy, pos = 1)
  assign(x = "alleles", value = LocusControl$alleles, pos = 1)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Join genotypes with extraction information (PLATE_ID) and make .gcl objects
  data_master <- dataAll %>% 
    dplyr::left_join(data_ex, by = c("SILLY_CODE" = "SILLY_CODE", "FK_FISH_ID" = "FISH_NO")) %>%
    tidyr::unite(SillySource, SILLY_CODE, FK_FISH_ID, sep = "_", remove = FALSE) #%>% 
  # dplyr::select(-PLOIDY) %>%
  # tidyr::spread(key = LOCUS, value = ALLELES)
  
  message("Data successfully pulled from LOKI, building SILLY.gcl objects\n")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Make .gcl objects by silly
  
  attnames <- c("FK_FISH_ID", "FK_COLLECTION_ID", "SILLY_CODE", "FK_PLATE_ID", "TISSUETYPE", "SillySource")
  
  for(silly in sillyvec) {  
    
    data_silly <- dplyr::filter(data_master, SILLY_CODE == silly)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Adaptation of LOKI2R to allow pulling all SILLYs just based on ProjectID
    if(nrow(data_silly) == 0) {
      warning(paste(silly, "was not found in LOKI, it may not have been imported!!!", sep = " "))
      message(paste0(silly,".gcl NOT FOUND IN LOKI ", match(silly,sillyvec)," of ",length(sillyvec)," FAILED!!!"))
      next
    }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ids <- as.character(sort(unique(data_silly$FK_FISH_ID)))
    
    nind <- length(ids)
    
    if(!nind){
      
      message0 <- paste0(silly," is empty.")
      
      message(message0)
      
      next
      
    }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Build scores array
    
    scores <- array(data = NA,
                    dim = c(nind, nloci, max(ploidy)), 
                    dimnames = list(ids, loci, paste0("Dose", seq(max(ploidy))))
    )
    
    sillyloci <- sort(unique(data_silly$LOCUS))  # These are the loci available for the current silly, this is needed to subset the tapply
    
    scores_allele_1 <- data_silly %>% 
      dplyr::filter(LOCUS %in% sillyloci) %>% 
      dplyr::select(FK_FISH_ID, LOCUS, ALLELE_1) %>% 
      dplyr::mutate(FK_FISH_ID = as.character(FK_FISH_ID)) %>% 
      tidyr::spread(LOCUS, ALLELE_1) %>% 
      tibble::column_to_rownames(var = "FK_FISH_ID") %>% 
      as.matrix()
    
    scores[ids, sillyloci, "Dose1"] <- scores_allele_1[ids, sillyloci]
    
    scores_allele_2 <- data_silly %>% 
      dplyr::filter(LOCUS %in% sillyloci) %>% 
      dplyr::select(FK_FISH_ID, LOCUS, ALLELE_2) %>% 
      dplyr::mutate(FK_FISH_ID = as.character(FK_FISH_ID)) %>% 
      tidyr::spread(LOCUS, ALLELE_2) %>% 
      tibble::column_to_rownames(var = "FK_FISH_ID") %>% 
      as.matrix()
    
    scores[ids, sillyloci, "Dose2"] <- scores_allele_2[ids, sillyloci]
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Build counts array
    
    counts <- array(data = NA, 
                    dim = c(nind, nloci, max(nalleles)), 
                    dimnames = list(ids, loci, paste0("Allele ", seq(max(nalleles))))
    )
    
    for(ind in ids) {
      
      for(locus in loci) {
        
        for(al in seq(nalleles[locus])) {
          
          counts[ind, locus, al] <- sum(scores[ind, locus, seq(ploidy[locus])] == alleles[[locus]][al])
          
        }  # al          
        
      }  # locus
      
      counts0 = counts[ind, , ]
      
      counts0[is.na(counts0)] <- 0
      
      zeroBOOL <- apply(counts0,1,sum) != ploidy
      
      counts[ind, zeroBOOL, ] <- NA
      
    }  # ind
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Build attributes data.frame
    
    attributes <- data_silly %>% 
      dplyr::select(attnames) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(FISH_ID = as.character(FK_FISH_ID),
                    FK_PLATE_ID = as.character(FK_PLATE_ID)) %>% 
      dplyr::rename(COLLECTION_ID = FK_COLLECTION_ID,
                    PLATE_ID = FK_PLATE_ID,
                    PK_TISSUE_TYPE = TISSUETYPE) %>% 
      dplyr::arrange(FK_FISH_ID) %>% 
      tibble::column_to_rownames(var = "FISH_ID")
    
    message(paste0(silly, ".gcl created ", match(silly, sillyvec), " of ", length(sillyvec), " completed."))
    
    assign(paste0(silly, ".gcl"), list(counts = counts, scores = scores, n = nind, attributes = attributes), pos = 1)
    
  }  # silly
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime) 
  
}