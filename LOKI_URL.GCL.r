LOKI_URL.GCL <- function(){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function builds the URL to connect to LOKI and is called on by other GCL 
  #  functions that connect to LOKI via RJDBC::dbConnect.
  # 
  #  Output:
  #  The URL string used by RJDBC::dbConnect to connect to the Oricle cloud database (LOKI)
  #
  #  Written by Andy Barclay 4/15/19
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
host <- "db-dfgcfresp.dfg.alaska.local"

port <- "1521"

svc <- "dfgcfresp.us1.ocm.s7134325.oraclecloudatcustomer.com"

protocol <-"jdbc:oracle:thin:"

url <- paste0(protocol,"@",host,":",port,"/",svc)

return(url)

}