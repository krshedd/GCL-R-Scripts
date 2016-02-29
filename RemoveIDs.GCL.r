RemoveIDs.GCL=function(silly,IDs){
##################################################################################################################################################################
#
#  Removes "IDs" (vial #'s as oppossed to positions) from the "*.gcl" object associated with "silly".
#
##################################################################################################################################################################

  if(!is.list(IDs)){
    IDs=list(IDs)  
  }
  names(IDs)=silly

  my.gcl=get(paste(silly,".gcl",sep=""),pos=1)

  scores=my.gcl$scores

  counts=my.gcl$counts
 
  attributes=my.gcl$attributes

  vials=as.character(dimnames(scores)[[1]])

  IND=!vials%in%IDs[[silly]]

  counts=counts[IND,,,drop=FALSE]

  scores=scores[IND,,,drop=FALSE]

  n=nrow(scores)

  attributes=attributes[IND,,drop=FALSE]

  assign(paste(silly,".gcl",sep=""),list(counts=counts,scores=scores,n=n,attributes=attributes),pos=1)

  return(NULL)
}