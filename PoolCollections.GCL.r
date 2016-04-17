PoolCollections.GCL=function(collections,loci,IDs=NULL,newname=paste(collections,collapse=".")){
######################################################################################################################################################################################
#
#  This function combines "*.gcl" objects into a new one called "newname.gcl".
#
#  "collections" is a vector of sillys (without ".gcl" extension) you wish to combine.
#
#  "loci" is the set of markers used for pooling.
#
#  "IDs" is a named list of fish ID character vectors, each vector is associated with and named after a member of "collections" 
#
#  "newname" is the name of the new "*.gcl" created. Do not provide ".gcl" extention. 
#
#  ReWritten by JJ 12/17/2015 
#  Edited by Kyle Shedd on 4/12/2016
#
######################################################################################################################################################################################

  while(!require(abind)){install.packages("abind")}

  bbind=function(...){abind(...,along=1)}

  if(nchar(newname)>200){newname=substr(newname, start=1, stop=200)}

  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}

  ncollections=length(collections)

  if(is.null(IDs)){IDs=lapply(collections,function(collection){rownames(get(paste0(collection,".gcl"),pos=1)$scores)})}

  if(!is.list(IDs)){stop("'IDs' must be a list")}

  if(ncollections!=length(IDs)){stop("'IDs' must be same length as 'collections'")}

  if(!is.character(unlist(IDs))) {stop("'IDs' must be a character vector, not a numeric vector")}

  names(IDs)=collections

  scores=Reduce(bbind,lapply(collections,function(collection){get(paste0(collection,".gcl"),pos=1)$scores[IDs[[collection]],loci,,drop=FALSE]}))

  n=nrow(scores)

  counts=Reduce(bbind,lapply(collections,function(collection){get(paste0(collection,".gcl"),pos=1)$counts[IDs[[collection]],loci,,drop=FALSE]}))

  attributes=Reduce(rbind,lapply(collections,function(collection){get(paste0(collection,".gcl"),pos=1)$attributes[IDs[[collection]],,drop=FALSE]}))

  rownames(scores)=rownames(counts)=rownames(attributes)=attributes$FK_FISH_ID=as.character(seq(n))

  assign(paste(newname,".gcl",sep=""),list(counts=counts,scores=scores,n=n,attributes=attributes),pos=1)  

  return(NULL)
}
