CreateLocusControl.GCL <- function(markersuite = NULL, locusnames = NULL, username, password) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # This function connects to LOKI and creates a "LocusControl" object  
  # User needs to supply either "markersuite" OR "locusnames" (character vector of loci)
  # 
  # "markersuite" is the pre-defined name in LOKI for the set of markers for which you want genotypes (e.g. markersuite="KenaiChinook2010_40SNPs"). 
  #   This set must be pre-defined in LOKI (see Eric Lardizabal). 
  #
  # "locusnames" is a character vector of marker names for which you want genotypes.
  #   i.e. c("One_CO1", "One_GPDH", "One_MHC2_190)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Example/Intended ####
  #
  # markersuite <- "Sockeye2011_96SNPs"
  #  
  # username <- "jjasper"; password <- "********"
  #
  # Written by JJ 10/05/2015
  # updated to new ojdbc6.jar path on V:/Analysis by Kyle Shedd on 05/05/2016
  # Updated to add 'locusnames' argument by Kyle Shedd on 01/10/2018
  # Updated by Andy Barclay 4/15/19; updated driver from ojdbc6.jar to ojdbc8.jar and changed the LOKI connection URL
  # to connect to the new Oracle cloud database
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(exists("LocusControl",where=1)){
        stop("LocusControl already exists")
    }

    while(!require(RJDBC)){install.packages("RJDBC")}

    if(!file.exists(path.expand("~/R"))){
  
     dir<-path.expand("~/R")
                                                   
     dir.create(dir)
      
     bool <- file.copy(from="V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar",to=path.expand("~/R/ojdbc8.jar"))
      
     } else {
      
     if(!file.exists(path.expand("~/R/ojdbc8.jar"))){
         
      bool <- file.copy(from="V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar",to=path.expand("~/R/ojdbc8.jar"))
      
      }
      
     }
 
    options(java.parameters = "-Xmx10g")

    if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {

      drv <- JDBC("oracle.jdbc.OracleDriver",classPath="C:/Program Files/R/RequiredLibraries/ojdbc8.jar"," ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
    
    } else {
      
      drv <- JDBC("oracle.jdbc.OracleDriver",classPath=path.expand("~/R/ojdbc8.jar")," ")
      
    }
    
    url <-LOKI_URL.GCL()
    
    con <- dbConnect(drv,url=url,user=username,password=password)
    
    # Query by 'markersuite', else query by 'locusnames'

    if(is.null(markersuite) & is.null(locusnames)) {stop("Need to provide either 'locusnames' or 'markersuite'")}
    
    if(is.null(locusnames) & !is.null(markersuite)) {
      
      lociqry <- paste("SELECT * FROM AKFINADM.V_LOCUSQRY WHERE SUITE_NAME = '",markersuite,"'",sep="")  # Query locus information of markers in markersuite.
      
      locidata <- dbGetQuery(con,lociqry)
      
      locusnames <- sort(locidata$LOCUS_NAME)

    } else {

      markersuite <- "User defined from locusnames"
      
      lociqry <- paste("SELECT * FROM AKFINADM.V_LOCUSQRY WHERE LOCUS_NAME IN (",paste0("'",locusnames,"'",collapse=","),")",sep="")  # Query locus information of markers in locusnames.
      
      locidata <- dbGetQuery(con,lociqry)
      
    }
    
  
    # Warn user if some 'locusnames' not found
    if(!all(locusnames %in% locidata$LOCUS_NAME)) {
    
      miss_loci <- locusnames[!locusnames %in% locidata$LOCUS_NAME]
      
      nmiss_loci <- length(miss_loci)
      
      message(paste(nmiss_loci, "out of", length(locusnames), "locusnames not found in LOKI!!!"))
        
      sapply(miss_loci, function(locus) {
        message(locus)
      } )
      
      locusnames <- locusnames[locusnames %in% locidata$LOCUS_NAME]
      
    }
    
    Publishedlocusnames <- setNames(locidata$PUBLISHED_NAME,locidata$LOCUS_NAME)[locusnames]
    
    nloci <- length(locusnames)
    
    ploidy <- setNames(locidata$PLOIDY,locidata$LOCUS_NAME)[locusnames]
    
    alleleqry <- paste("SELECT * FROM AKFINADM.V_ALLELEQRY WHERE LOCUS IN (",paste0("'",locusnames,"'",collapse=","),")",sep="")#Get possible alleles from allele lookup table.

    alleles0 <- dbGetQuery(con,alleleqry)

    discon <- dbDisconnect(con)

    alleles <- tapply(X=alleles0$VALUE,INDEX=list(alleles0$LOCUS),FUN=sort,simplify=FALSE)[locusnames]
    
    nalleles <- sapply(alleles,length)
      
    assign("LocusControl",list(MarkerSuite=markersuite,locusnames=locusnames,Publishedlocusnames=Publishedlocusnames,alleles=alleles,nalleles=nalleles,ploidy=ploidy),pos=1) #Assign elements to LocusControl list.	

    ans <- paste0("LocusControl created for markersuite='",markersuite,"'")
    
    return(ans)

}  
   

