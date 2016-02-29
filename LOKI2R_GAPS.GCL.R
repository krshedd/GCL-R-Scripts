LOKI2R_GAPS.GCL=function(sillyvec,username,password){

######################################################################################################################################################################################
################  :WARNING:WARNING:WARNING:WARNING:WARNING:   ########################################################################################################################
################ 	THIS FUNCTION REQUIRES THE PACKAGE "RODBC"  ########################################################################################################################
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
#  Written by AB, EL, & JJ,  10/06/2015
######################################################################################################################################################################################

  if(!exists("LocusControl")){

    stop("'LocusControl' not yet built.")

  }

    if(!file.exists(path.expand("~/R"))){
  
     dir<-path.expand("~/R")
                                                   
     dir.create(dir)
      
     bool <- file.copy(from="V:/DATA/R_GEN/JJs GCL/jars/ojdbc6.jar",to=path.expand("~/R/ojdbc6.jar"))
      
     }
      
    if(file.exists(path.expand("~/R"))){
      
     if(!file.exists(path.expand("~/R/ojdbc6.jar"))){
         
      bool <- file.copy(from="V:/DATA/R_GEN/JJs GCL/jars/ojdbc6.jar",to=path.expand("~/R/ojdbc6.jar"))
      
      }
      
     }
                                                                  
  while(!require(RJDBC)){install.packages("JDBC")}

# while(!require(pryr)){install.packages("pryr")}
    
  start.time <- Sys.time() 
     
  options(java.parameters = "-Xmx10g")

  drv <- JDBC("oracle.jdbc.OracleDriver",classPath="C:/Program Files/R/RequiredLibraries/ojdbc6.jar"," ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
 
  con <- dbConnect(drv, "jdbc:oracle:thin:@(DESCRIPTION=(ADDRESS_LIST=(ADDRESS=(PROTOCOL=TCP)(HOST=db-pcfres.dfg.alaska.local)(PORT=1521)))(CONNECT_DATA=(SID=PCFRES)))",username,password)

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
  
    sillydata <- subset(dataAll,SILLY_CODE==silly)
   
    ids <- as.character(sort(unique(sillydata$FISH_ID)))   
      
    sillyvials <- paste(silly,ids,sep="_")
      
    nind <- length(sillyvials)
      
    scores0 <- array(data=NA,dim=c(nind,length(loci),max(ploidy)),dimnames=list(ids,loci,paste0("Dose",seq(max(ploidy)))))
    
    sillyloci <- sort(unique(sillydata$LOCUS)) #These are the loci available for the current silly, this is needed to subset the tapply
    
    scores0[ids,sillyloci,1] <- tapply(X=as.vector(sillydata$ALLELE_1),INDEX=list(sillydata$FISH_ID,sillydata$LOCUS),FUN=unique)[ids,sillyloci]
      
    scores0[ids,sillyloci,2] <- tapply(X=as.vector(sillydata$ALLELE_2),INDEX=list(sillydata$FISH_ID,sillydata$LOCUS),FUN=unique)[ids,sillyloci]

    
      
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
      
    message(paste0(silly,".gcl created ", match(silly,sillyvec)," of ",length(sillyvec)," completed."))
      
    assign(paste(silly,".gcl",sep=""),list(counts=counts,scores=scores,n=nind,attributes=attributes),pos=1)

  }#silly 

  stop.time <- Sys.time()

  fulltime <- stop.time - start.time

  print(fulltime) 
      
  return(fulltime)
}  
   

