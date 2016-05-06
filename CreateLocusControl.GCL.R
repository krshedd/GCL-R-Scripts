CreateLocusControl.GCL=function(markersuite,username,password)
{
######################################################################################################################################################################################
########################################################################################################################
#
#  This function connects to LOKI and creates a "LocusControl" object  
#   
#
#    "markersuite" is the pre-defined name in LOKI for the set of markers for which you want genotypes (e.g. markersuite="KenaiChinook2010_40SNPs"). This set must be pre-defined in LOKI (see Eric Lardizabal). 
#
######## Example/Intended ############################################################################################################################################################
#
#  markersuite <- "Sockeye2011_96SNPs"
#
#  
#  username <- "jjasper"; password <- "********"
#
#  CreateLocusControl.GCL(markersuite,username,password)
#
#  Written by JJ 10/05/2015; updated to new ojdbc6.jar path on V:/Analysis by Kyle Shedd on 05/05/2016
######################################################################################################################################################################################

    if(exists("LocusControl",where=1)){
        stop("LocusControl already exists")
    }

    while(!require(RJDBC)){install.packages("RJDBC")}

    if(!file.exists(path.expand("~/R"))){
  
     dir<-path.expand("~/R")
                                                   
     dir.create(dir)
      
     bool <- file.copy(from="V:/Analysis/R files/Scripts/DEV/jars/ojdbc6.jar",to=path.expand("~/R/ojdbc6.jar"))
      
     } else {
      
     if(!file.exists(path.expand("~/R/ojdbc6.jar"))){
         
      bool <- file.copy(from="V:/Analysis/R files/Scripts/DEV/jars/ojdbc6.jar",to=path.expand("~/R/ojdbc6.jar"))
      
      }
      
     }
 
    options(java.parameters = "-Xmx10g")

    if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc6.jar")) {

      drv <- JDBC("oracle.jdbc.OracleDriver",classPath="C:/Program Files/R/RequiredLibraries/ojdbc6.jar"," ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
    
    } else {
      
      drv <- JDBC("oracle.jdbc.OracleDriver",classPath=path.expand("~/R/ojdbc6.jar")," ")
      
    }
    
    con <- dbConnect(drv, "jdbc:oracle:thin:@(DESCRIPTION=(ADDRESS_LIST=(ADDRESS=(PROTOCOL=TCP)(HOST=db-pcfres.dfg.alaska.local)(PORT=1521)))(CONNECT_DATA=(SID=PCFRES)))",username,password)

    lociqry <- paste("SELECT * FROM AKFINADM.V_LOCUSQRY WHERE SUITE_NAME = '",markersuite,"'",sep="")#Query locus information of markers in markersuite.

    locidata <- dbGetQuery(con,lociqry)

    locusnames <- sort(locidata$LOCUS_NAME)

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
   

