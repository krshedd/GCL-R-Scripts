AttributesToIDs.GCL=function(silly,attribute,matching){
########################################################################################################################################################
#
#  This function converts attributes into IDs.
#
#  "silly" is the silly without the ".gcl" extension that the user wishes to query.
#
#  "attribute" should match one of the following:
########################################################################################################################################################
##  "FK_COLLECTION_ID","FK_FISH_ID","PK_TISSUE_TYPE","CAPTURE_LOCATION","CAPTURE_DATE","VIAL_BARCODE","EXHAUSTED_HOW","EXHAUSTED_BY","EXHAUSTED_DATE"  
##
##  "STORAGE_ID","UNIT","SHELF_RACK","SLOT","END_CAPTURE_DATE","MESH_SIZE","plateID","SillySource"     
########################################################################################################################################################
#
#  "matching" is a vector of attributes of the same mode as "attribute" that correspond to the desired IDs. 
#
# Written by Jim Jasper 6/15/2011
########################################################################################################################################################
########################################################################################################################################################


  IDs = NULL
  myattribute=as.character(get(paste(silly,".gcl",sep=""),pos=1)$attributes[[attribute]])
   for(match in 1:length(matching)){
     IND=myattribute==matching[match]
     IDs=c(IDs,as.character(dimnames(get(paste(silly,".gcl",sep=""),pos=1)$scores)[[1]][IND]))
     }
  return(IDs)
}