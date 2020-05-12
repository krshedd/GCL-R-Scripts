GTSeq_Metadata_Report.GCL <- function(project_name, dir, file_name, username, password, open.file = FALSE){  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function pulls a genotypes report from LOKI and writes the data to aUTF-8 coded csv file.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - character vector of SILLYs - example: c("SCIMA18","SCIMA17")
  #
  #   loci - character vector of locus names as they are spelled in LOKI - example: c("One_E2","One_MHC2_251","One_Cytb_17")
  #
  #   dir - the directory for writing output file with "\\" separator between folders example: "V:\\Analysis\\Staff\\Andy Barclay\\R"
  #
  #   file_name - the name of the csv file. Reports that exceed the maximum number of columns for Exel will be plit into multilple files with a number added to the file name (e.g. test_metadat_1.csv, test_metadat_2.csv)
  #
  #   username - your user name for accessing LOKI through R example: "awbarclay"
  #
  #   password - your password for accessing LOKI through R example: "password"
  #
  #   project_name - a character vector of one or more project names as spelled in LOKI example: c("P014","P015","P016")
  #
  #   open.file - logical,if set to TRUE the CSV file will open after it has been written. Default is FALSE.
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Returns a tibble of metadata by locus (crosstab/wide format) and writes out to a .csv file.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~  # 
  #  metadata <- GTSeq_Metadata_Report.GCL(project_name = "K158", dir = "C:\\Users\\awbarclay\\Documents", file_name = "test_metadat", username = "awbarclay", password = "************", open.file = FALSE)
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Recording function start time
  start.time <- Sys.time()
  
  # Checking to make sure a file path is supplied
  if(!exists("path")|!length(grep("*.csv",path,value=FALSE))==1){stop("The user must supply a file path with csv extension for writing out the metadata table.")}
  
  # Requiring the tidyverse and RJDBC packages. Packages will be installed if the user doesn't have them installed
  while(!require(tidyverse)){install.packages("tidyverse") }
  while(!require(RJDBC)){install.packages("RJDBC")}
  while(!require(tools)){install.packages("tools")}
  
  # Making sure the ojdbc jar file exists on the users C drive. 
  # If no jar file exists it will be copied from the v drive to the appropriate location on the users C drive.
  if(!file.exists(path.expand("~/R"))){
    dir<-path.expand("~/R")
    dir.create(dir)
    bool <- file.copy(from="V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar",to=path.expand("~/R/ojdbc8.jar"))
  } else {
    if(!file.exists(path.expand("~/R/ojdbc8.jar"))){
      bool <- file.copy(from="V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar",to=path.expand("~/R/ojdbc8.jar"))
    }
  }
  
  # Setting default java.parameters
  options(java.parameters = "-Xmx10g")
  
  if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {
    drv <- JDBC("oracle.jdbc.OracleDriver",classPath="C:/Program Files/R/RequiredLibraries/ojdbc8.jar"," ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
  } else {
    drv <- JDBC("oracle.jdbc.OracleDriver",classPath=path.expand("~/R/ojdbc8.jar")," ")
  }
  
  url <-LOKI_URL.GCL()
  
  # Connect to LOKI
  con <- dbConnect(drv,url=url,user=username,password=password)

  # Creating java query
  
  gnoqry <- paste0("SELECT LAB_PROJECT_NAME, SILLY_CODE, FK_FISH_ID, LOCUS, POSITIONS, HAPLO_ALLELES, HAPLO_COUNTS, GENOTYPE, SNP_ALLELES, PROBES, FWD_PRIMER_SEQ, DNA_PLATE_ID, DNA_PLATE_WELL_POS FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO GENO WHERE EXISTS (SELECT * FROM AKFINADM.V_LAB_PROJECT_WELL LPW WHERE LPW.LAB_PROJECT_NAME IN (", paste0("'", project_name, "'", collapse = ","), ") AND LPW.SILLY_CODE = GENO.SILLY_CODE AND LPW.FISH_NO = GENO.FK_FISH_ID)")

  # Pull genotypes and concatenate alleles into one column with "/" separator
  dataAll0 <- RJDBC::dbGetQuery(con, gnoqry) %>% 
    dplyr::as_tibble() %>% 
      dplyr::filter(LAB_PROJECT_NAME == project_name)

  # Disconnect from LOKI
  discon <- dbDisconnect(con)
  
  # Put in crosstab format
  dataAll <- dataAll0 %>% 
    arrange(LOCUS)

  # Checking to see if no data was pulled for a given project. If true, the function will stop and print an error message in the console.
  if(nrow(dataAll)==0 & exists("project_name")){
    
    stop("No data exist for the supplied lab project name. Check to make sure the name spelled correctly. Also, project names for genotypes imported prior to October 2016 may not be included in LOKI lookup table.")
    
  }
  
  # Write data to an excel csv file with no more than 1,048,575 rows
  loc <- dataAll %>% 
    pull(LOCUS) %>% 
    unique()
  
  n_indiv <- dim(dataAll %>% filter(LOCUS == loc[1]))[1]
  
  xls <- 1048575 #for making sure final file is no longer than the Excel row limit
  
  nloci <- floor(xls/n_indiv)#Maximim number of loci to include in each file
  
  rows <- nloci*n_indiv #Max number of rows per file so loci aren't plit among multiple files
    
  end <- length(dataAll$FK_FISH_ID)
  
  N <- ceiling(end/rows)#Number of files to write
  
  for (i in 1:N){
    
    start <-rows*i-rows+1
    
    stop <- if(i*rows>end){end} else(i*rows)
    
    readr::write_csv(dataAll[start:stop, ], path = paste(dir,"/", file_name, "_", i, ".csv", sep=''))
    
  }

  # Open CSV file in excel if open.file is set to TRUE
  if(open.file == TRUE){shell(path, wait = FALSE)}
  
  # Calculate the time it took the function to pull and write data and print time to console
  stop.time <- Sys.time() 
  
  message(paste0("CSV file writen to:", file_path_as_absolute(path)))
  print(stop.time-start.time)
  
  # Return data to assign to an object
  return(dataAll)
  
}

