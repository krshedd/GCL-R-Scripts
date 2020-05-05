PoolCollections.GCL=function(collections,loci,IDs=NULL,newname=paste(collections,collapse=".")){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function combines "*.gcl" objects into a new one called "newname.gcl".
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   collections - a character vector of silly codes without the ".gcl" extention (e.g. collections <- c("KQUART06","KQUART08","KQUART09")). 
  #                 Collections can be a single silly if you want create a new ".gcl" with only fish supplied in IDs.
  #
  #   loci - a character vector of locus names
  # 
  #   IDs - a named list of fish ID character vectors, each vector is associated with and named after a member of "collections".
  #         These will be used to subset each collection before pooling. If no IDs are supplied all individuals from each collection are used.
  #
  #   newname - is the name of the new "*.gcl" created. Do not provide ".gcl" extention. If no name supplied then the newname defaults to
  #             the collection names collapsed with a period between each name (e.g. "KQUART06.KQUART08.KQUART09")
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Assigns a new "pooled collection" to your workspace
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  # 
  #   removedInd <- (collections = c("KQUART06","KQUART08","KQUART09"), loci = loci, IDs=NULL, newname=paste(collections,collapse="."))
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function is also useful for producing "pooled mixture" objects for mixed stock analysis. 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

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

  attributes=Reduce(rbind,lapply(collections,function(collection){get(paste0(collection,".gcl"),pos=1)$attributes[which(get(paste0(collection,".gcl"),pos=1)$attributes$FK_FISH_ID %in% IDs[[collection]]),,drop=FALSE]}))

  rownames(scores)=rownames(counts)=rownames(attributes)=attributes$FK_FISH_ID=as.character(seq(n))

  assign(paste(newname,".gcl",sep=""),list(counts=counts,scores=scores,n=n,attributes=attributes),pos=1)  

  return(NULL)
}
