GenotypeReport.GCL=function(sillyvec=NULL,loci=NULL,path,username,password,project_name=NULL,open.file=FALSE){  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function pulls a geneotypes report from LOKI and writes the data to aUTF-8 coded csv file. 
  #
  #  **Warning - genotypes (pre October 2016) may not exist in the LOKI lookup table with a project name, so genotypes will have to be pulled by sillyvec and loci**  
  #
  #  **Note: Excel 2010 is limited to 1,048,576 rows by 16,384 columns per worksheet. This limitation may have to be addressed in the future.**  
  #
  #  Argument(s): 
  #
  #  sillyvec <- character vector of SILLYs - example: c("SCIMA18","SCIMA17")
  #
  #  loci <- charcter vector of locus names as they are spelled in LOKI - example: c("One_E2","One_MHC2_251","One_Cytb_17")
  #
  #  path <- file path with csv extension for writting output file with "\\" separator between folders example: "V:\\Analysis\\Staff\\Andy Barclay\\R\Genotypes report"
  #
  #  username <- your user name for accessing LOKI through R example: "awbarclay"
  #
  #  password <- your password for accessing LOKI through R example: "password"
  #
  #  project_name <- a character vector of one or more project names as spelled in LOKI example: c("P014","P015","P016")
  #
  #  open.file <- logical,if set to TRUE the CSV file will open after it has been written. Default is FALSE.
  #
  #  Output:
  #  Crosstab format of silly individuals by loci with both alleles in one column
  #
  #  Written by Andy Barclay 10/5/18
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #Recording function start time
  start.time <- Sys.time()
  
  #Checking to make sure the correct combination of arguments is being use. If wrong, the function will stop and print an error message to the console.
  if(is.null(sillyvec)&is.null(loci)&is.null(project_name)|!is.null(sillyvec)&!is.null(loci)&!is.null(project_name)|!is.null(sillyvec)&is.null(loci)&!is.null(project_name)|is.null(sillyvec)&!is.null(loci)&!is.null(project_name)|is.null(sillyvec)&!is.null(loci)&is.null(project_name)){
    
    stop("The user must supply one of the following argument combinations: sillyvec (for all loci and individuals for each silly), sillyvec and loci (all individuals for supplied locus list), or project_name (for all individuals and loci in a given project)")
    
    }
  
  #Checking to make sure a file path is supplied
  if(!exists("path")|!length(grep("*.csv",path,value=FALSE))==1){stop("The user must supply a file path with csv extension for writing out genotypes table.")}
  
  #Requiring the tidyverse and RJDBC packages. Packages will be installed if the user doesn't have them installed
  while(!require(tidyverse)){install.packages("tidyverse") }
  while(!require(RJDBC)){install.packages("RJDBC")}
  while(!require(tools)){install.packages("tools")}
  
  #Making sure the ojdbc jar file exists on the users C drive. 
  #If no jar file exists it will be copied from the v drive to the appropriate locaiton on the users C drive.
  if(!file.exists(path.expand("~/R"))){
    
    dir<-path.expand("~/R")
    
    dir.create(dir)
    
    bool <- file.copy(from="V:/Analysis/R files/Scripts/DEV/jars/ojdbc6.jar",to=path.expand("~/R/ojdbc6.jar"))
    
  } else {
    
    if(!file.exists(path.expand("~/R/ojdbc6.jar"))){
      
      bool <- file.copy(from="V:/Analysis/R files/Scripts/DEV/jars/ojdbc6.jar",to=path.expand("~/R/ojdbc6.jar"))
      
    }
    
  }
  
  #Setting default java.parameters
  options(java.parameters = "-Xmx10g")
  
  if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc6.jar")) {
    
    drv <- JDBC("oracle.jdbc.OracleDriver",classPath="C:/Program Files/R/RequiredLibraries/ojdbc6.jar"," ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
    
  } else {
    
    drv <- JDBC("oracle.jdbc.OracleDriver",classPath=path.expand("~/R/ojdbc6.jar")," ")
    
  }
  
  #Connect to LOKI
  con <- dbConnect(drv, "jdbc:oracle:thin:@(DESCRIPTION=(ADDRESS_LIST=(ADDRESS=(PROTOCOL=TCP)(HOST=db-pcfres.dfg.alaska.local)(PORT=1521)))(CONNECT_DATA=(SID=PCFRES)))",username,password)
  
  #Creating java query when sillyvec and loci are supplied.  
  if(!is.null(sillyvec) & !is.null(loci)){
    
    #PUlling haploid and diploid loci separetly because haploid only have ALLELE_1.
    
    gnoqry_d <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1, ALLELE_2 FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE LOCUS IN (" ,paste0("'",loci,"'",collapse=","), ") AND SILLY_CODE IN", "(",paste0("'",sillyvec,"'",collapse=","),") AND PLOIDY = 'D'")
    
    gnoqry_h <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1 FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE LOCUS IN (" ,paste0("'",loci,"'",collapse=","), ") AND SILLY_CODE IN", "(",paste0("'",sillyvec,"'",collapse=","),") AND PLOIDY = 'H'")
    
  } 
  
  #Creating java query when only sillyvec is supplied
  if(!is.null(sillyvec) & is.null(loci)){
    
    #PUlling haploid and diploid loci separetly because haploid only have ALLELE_1.
    
    gnoqry_d <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1, ALLELE_2 FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE SILLY_CODE IN", "(",paste0("'",sillyvec,"'",collapse=","),") AND PLOIDY = 'D'")
    
    gnoqry_h <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1 FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE SILLY_CODE IN", "(",paste0("'",sillyvec,"'",collapse=","),") AND PLOIDY = 'H'")
    
  }
  
  #Creating java query when only project_name is supplied.
  if(!is.null(project_name)){
    
    #PUlling haploid and diploid loci separetly because haploid only have ALLELE_1.
    
    gnoqry_d <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1, ALLELE_2 FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE LAB_PROJECT_NAME IN", "(",paste0("'",project_name,"'",collapse=","),") AND PLOIDY = 'D'")
    
    gnoqry_h <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1 FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE LAB_PROJECT_NAME IN", "(",paste0("'",project_name,"'",collapse=","),") AND PLOIDY = 'H'")
    
  }
  
  #Pull diploid data and concatinate alleles into one column with "/" separator
  data_d <- dbGetQuery(con,gnoqry_d) %>% 
    as_tibble() %>% 
    mutate(ALLELES = paste(ALLELE_1,ALLELE_2,sep="/")) %>% 
    select(-ALLELE_1,-ALLELE_2) 

  #Pull haploid data
  data_h <- dbGetQuery(con,gnoqry_h) %>% 
    as_tibble()%>% 
    mutate(ALLELES=ALLELE_1) %>% 
    select(-ALLELE_1)
  
  #Disconnect from LOKI
  discon <- dbDisconnect(con)
  
  #Combine diploid and haploid marker into one data frame and put into crosstab format (i.e., unstack table)
  dataAll <- bind_rows(data_d,data_h) %>% 
    select(-PLOIDY) %>% 
    spread(key=LOCUS,value=ALLELES)
  
  #Checking to see if no data was pulled for a given project. If true, the function will stop and print an erorr message in the console.
  if(nrow(dataAll)==0 & exists("project_name")){
    
    stop("No data exist for the supplied lab project name. Check to make sure the name spelled correctly. Also, project names for genotypes imported prior to October 2016 may not be included in LOKI lookup table.")
   
     }
    
  #Write data to an excel csv file
  dataAll %>% 
    write_excel_csv(path=path)
  
  #Open CSV file in excel if open.file is set to TRUE
  if(open.file==TRUE){shell(path,wait=FALSE)}
  
  #Calculate the time it took the function to pull and write data and print time to console
  stop.time <- Sys.time() 
  
  message(paste0("CSV file writen to:", file_path_as_absolute(path)))
  print(stop.time-start.time)
  
  #Return data to assign to an object
  return(dataAll)
   
}

