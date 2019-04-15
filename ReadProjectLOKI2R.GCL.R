ReadProjectLOKI2R.GCL <- function(projectID, username, password){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ####  :WARNING:WARNING:WARNING:WARNING:WARNING:  ####
  ####  THIS FUNCTION REQUIRES THE PACKAGES "RJDBC"  ####
  #
  #  You need to have LocusControl first. Run "CreateLocusControl.GCL" for your markersuite.
  #
  #  You need to provide your username and password to access LOKI.
  #
  #  This function connects to LOKI and creates a "*.gcl" object for each collection in a lab project.
  #
  #    "projectID"  a unique number assigned to a project when it is built in LOKI. If unsure, check iStrategy.
  #
  #  Written by Andy Barclay 9/27/12
  #
  #  Modified by Zac Grauvogel for use with all Oracle Accounts  on 4/17/14
  #
  #  Updated by Kyle Shedd with code from Eric Lardizabal to change sql query for markersuite on 10/15/15
  #  Also changed from using RODBC (ReadLOKI) to JDBC (LOKI2R), code is direct copy/paste from LOKI2R
  #  Updated to new ojdbc6.jar path on V:/Analysis by Kyle Shedd on 05/05/2016
  #  Updated by Andy Barclay 4/15/19; updated driver from ojdbc6.jar to ojdbc8.jar and changed the LOKI connection URL
  #  to connect to the new Oracle cloud database
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!exists("LocusControl")){
    
    stop("'LocusControl' not yet built.")
    
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
  
  message("2 main function tasks: \n1) Pull data from LOKI\n2) Create SILLY.gcl objects\n")
  
  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {
    
    drv <- JDBC("oracle.jdbc.OracleDriver",classPath="C:/Program Files/R/RequiredLibraries/ojdbc8.jar"," ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
    
  } else {
    
    drv <- JDBC("oracle.jdbc.OracleDriver",classPath=path.expand("~/R/ojdbc8.jar")," ")
    
  }
  
  url <-LOKI_URL.GCL()
  
  con <- dbConnect(drv,url=url,user=username,password=password)
  
  loci <- LocusControl$locusnames
  
  nloci <- length(loci)
  
  ploidy <- LocusControl$ploidy
  
  alleles <- LocusControl$alleles
  
  nalleles <- LocusControl$nalleles
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Adaptation of LOKI2R to allow pulling all SILLYs just based on ProjectID
  projlocqry <- paste("SELECT * FROM AKFINADM.LAB_PROJECT_LOCUS WHERE LAB_PROJECT_ID = '", projectID, "'", sep = "")
  
  projlodata <- dbGetQuery(con, projlocqry)
  
  sillyvec <- unique(projlodata$SILLY_CODE)
  
  assign(x = "ProjectSillys", value = sillyvec, pos = 1)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Adaptation of LOKI2R to subset only fish in SILLYs run for ProjectID
  projfshqry <- paste("SELECT * FROM AKFINADM.V_LAB_PROJECT_WELL WHERE LU_LAB_PROJECT_ID = '", projectID, "'", sep = "")
  
  projfshdata <- dbGetQuery(con, projfshqry)
  
  fish.sillyvec <- sapply(unique(projfshdata$SILLY_CODE), function(silly) {
    as.vector(t(subset(projfshdata, SILLY_CODE == silly, select = FISH_NO)))
  }, simplify = FALSE)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gnoqry <- paste("SELECT * FROM AKFINADM.V_GNOQRY WHERE LOCUS IN (",paste0("'",loci,"'",collapse=","),") AND SILLY_CODE IN (",paste0("'",sillyvec,"'",collapse=","),")",sep="")
  
  dataAll0 <- dbGetQuery(con,gnoqry)   
  
  dataAllbool <- dataAll0$PK_TISSUE_TYPE==dataAll0$PREFERRED_TISSUE 
  
  bothNAbool <- is.na(dataAll0$PREFERRED_TISSUE) &  is.na(dataAll0$PK_TISSUE_TYPE)
  
  dataAllbool[is.na(dataAllbool)] <- FALSE 
  
  dataAllbool[bothNAbool] <- TRUE  
  
  dataAll <- dataAll0[dataAllbool,]
  
  attrCol <- c("TEST_TYPE","LOCUS","ALLELE_1","ALLELE_2","SUITE_NAME","PREFERRED_TISSUE") 
  
  attnames <- names(dataAll)[!names(dataAll)%in%attrCol]
  
  discon <- dbDisconnect(con)
  
  message("Data successfully pulled from LOKI, building SILLY.gcl objects\n")
  
  for(silly in sillyvec){  
    
    sillydata <- subset(dataAll,SILLY_CODE==silly)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Adaptation of LOKI2R to allow pulling all SILLYs just based on ProjectID
    if(dim(sillydata)[1] == 0) {
      warning(paste(silly, "was not found in LOKI, it may not have been imported!!!", sep = " "))
      message(paste0(silly,".gcl NOT FOUND IN LOKI ", match(silly,sillyvec)," of ",length(sillyvec)," FAILED!!!"))
      next
    }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ids <- as.character(sort(unique(fish.sillyvec[[silly]])))  # ids <- as.character(sort(unique(sillydata$FISH_ID)))   
    
    sillyvials <- paste(silly,ids,sep="_")
    
    nind <- length(sillyvials)
    
    if(!nind){
      
      message0 <- paste0(silly," is empty.")
      
      message(message0)
      
      next
      
    }
    
    scores <- array(data=NA,dim=c(nind,length(loci),max(ploidy)),dimnames=list(ids,loci,paste0("Dose",seq(max(ploidy)))))
    
    sillyloci <- sort(unique(sillydata$LOCUS)) #These are the loci available for the current silly, this is needed to subset the tapply
    
    scores[ids,sillyloci,1] <- tapply(X=as.vector(sillydata$ALLELE_1),INDEX=list(sillydata$FISH_ID,sillydata$LOCUS),FUN=unique)[ids,sillyloci]
    
    scores[ids,sillyloci,2] <- tapply(X=as.vector(sillydata$ALLELE_2),INDEX=list(sillydata$FISH_ID,sillydata$LOCUS),FUN=unique)[ids,sillyloci]
    
    counts <- array(NA,c(nind,nloci,max(nalleles)),dimnames=list(ids,loci,paste("Allele ",1:max(nalleles),sep="")))
    
    for(ind in ids){
      
      for(locus in loci){
        
        for(al in seq(nalleles[locus])){
          
          counts[ind,locus,al] <- sum(scores[ind,locus,seq(ploidy[locus])]==alleles[[locus]][al])
          
        }#al          
        
      }#locus
      
      counts0=counts[ind,,]
      
      counts0[is.na(counts0)] <- 0
      
      zeroBOOL <- apply(counts0,1,sum)!=ploidy
      
      counts[ind,zeroBOOL,] <- NA
      
    }#ind
    
    attributes <- sillydata[,attnames]
    
    attributes <- aggregate(attributes[,-match("FISH_ID",attnames)],list(FISH_ID=attributes$FISH_ID),unique)
    
    attributes <- subset(attributes, FISH_ID %in% as.numeric(ids))
    
    attributes <-  data.frame(attributes,SillySource=sillyvials,row.names=ids,stringsAsFactors = FALSE)
    
    names(attributes)[names(attributes)=="FISH_ID"] <- "FK_FISH_ID"   
    
    attributes$CAPTURE_DATE <- as.POSIXct(attributes$CAPTURE_DATE)
    
    attributes$END_CAPTURE_DATE <- as.POSIXct(attributes$END_CAPTURE_DATE)
    
    message(paste0(silly,".gcl created ", match(silly,sillyvec)," of ",length(sillyvec)," completed."))
    
    assign(paste(silly,".gcl",sep=""),list(counts=counts,scores=scores,n=nind,attributes=attributes),pos=1)
    
  }#silly
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime) 
  
  return(fulltime)
}