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
  #                     CAPTURE_DATE <dttm>: the date each individual was captured in POSIXct, format (e.g. May 5, 2020 = "2020-05-05") 
  #                     END_CAPTURE_DATE <dttm>: the last collection date for a silly in POSIXct, format (e.g. May 5, 2020 = "2020-05-05") 
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
  #                     SILLY_SOURCE <doupble>: the original silly code and fish ID for each individual (e.g. KQUART06_1). When pulled from loki this will be the SILLY_CODE and FK_FISH_ID
  #                               
  #               Columns 20+
  #                     The remaining columns in the object will be the scores for all loci in the LocusControl object. 
  #                     Each locus will have a column for each dose. The columns will be named after the locus with a number added to the locus name after dose 1 (e.g. dose 1 = GTH2B-550; dose 2 = GTH2B-550.1) 
  #
  #
  #   The tibbles will be named after the silly code with a .gcl extention (e.g. KQUART06.gcl)
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   password = "************"
  # 
  #   LOKI2R.GCL(sillyvec = sillyvec, username = "awbarclay", password = password)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
######################################################################################################################################################################################
#
#  This function connects to LOKI and creates a "*.gcl" object for each silly in sillyvec.  
# 
#  A "*.gcl" object is a Gene Conservation Laboratory genotypes object with associated sample attributes.  
#
#    "sillyvec" is a vector of silly codes you want to pull from LOKI (e.g. sillyvec=c("KQUART06","KQUART08","KQUART09")).
#
######## Example/Intended ############################################################################################################################################################
#
#  sillyvec <- "KSALM95"
#  
#  username <- "jjasper"; password <- "********"
#
#  markersuite="Chinook2011_52SNPs"
#
#  CreateLocusControl.GCL(markersuite,username,password)
#
#  LOKI2R.GCL(sillyvec,username,password)
#
#  Written by AB, EL, & JJ,  05/05/2016 # from.jar="V:/Analysis/R files/Scripts/DEV/jars/ojdbc6.jar" 2016-02-13 10:40:33 AKST for KS
#  Updated by Andy Barclay 4/15/19; updated driver from ojdbc6.jar to ojdbc8.jar and changed the LOKI connection URL
#  to connect to the new Oracle cloud database.
######################################################################################################################################################################################

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
    
  for(silly in sillyvec){ 

    message0 <- paste0(silly,".gcl created ", match(silly,sillyvec)," of ",length(sillyvec)," completed.") 
  
    sillydata <- subset(dataAll,SILLY_CODE==silly)
   
    ids <- as.character(sort(unique(sillydata$FISH_ID)))   
      
    sillyvials <- paste(silly,ids,sep="_")
      
    nind <- length(sillyvials)

    if(length(ids) == 0){

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
      
    attributes <-  data.frame(attributes,SillySource=sillyvials,row.names=ids,stringsAsFactors = FALSE)

    names(attributes)[names(attributes)=="FISH_ID"] <- "FK_FISH_ID"   
    
    attributes$CAPTURE_DATE <- as.POSIXct(attributes$CAPTURE_DATE)
    
    attributes$END_CAPTURE_DATE <- as.POSIXct(attributes$END_CAPTURE_DATE)
      
    message(message0)
      
    assign(paste(silly,".gcl",sep=""),list(counts=counts,scores=scores,n=nind,attributes=attributes),pos=1)

  }#silly 

  stop.time <- Sys.time()

  fulltime <- stop.time - start.time

  print(fulltime) 
      
  return(fulltime)
}  
   

