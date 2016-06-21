DupCheckBetweenSillys.GCL=function(KeySillys,KeySillyIDs=NULL,BetweenSillys,loci,threshold){
##################################################################################################################################################################################################################################################################################################
#
#  "KeySillyIDs" is a named list, same length and names as "KeySillys" 
#
#  e.g. KeySillys="SAUKE09QC"
#
#       KeySillyIDs=list(SAUKE09QC=c("2587","3007"))
#
#  ReWritten by JJ 12/30/2015
#  Revised by KS on 06/20/2016 in order to accomodate multiple KeySillys for the output  
##################################################################################################################################################################################################################################################################################################

  while(!require(abind)){install.packages("abind")}

  while(!require(outliers)){install.packages("outliers")}

  keygcls=sapply(KeySillys,function(silly){assign(paste("new",silly,".gcl",sep=""),get(paste(silly,".gcl",sep=""),pos=1),pos=1)},simplify=FALSE)
 
  ploidy=LocusControl$ploidy[loci]
  
  nloci=length(loci)

  if(!is.null(KeySillyIDs)){

    names(KeySillyIDs)=KeySillys

    KeySillyIDsALL=sapply(KeySillys,function(silly){rownames(keygcls[[silly]]$counts)},simplify=FALSE)
   
    KeySillyIDs2Remove=sapply(KeySillys,function(silly){list(KeySillyIDsALL[[silly]][is.na(match(KeySillyIDsALL[[silly]],KeySillyIDs[[silly]]))])},simplify=FALSE)

    EMPTY=sapply(KeySillys,function(silly){RemoveIDs.GCL(paste0("new",silly),KeySillyIDs2Remove[[silly]])})

  }

  newKeySillys=paste0("new",KeySillys)

  newname=paste0(newKeySillys,collapse=".")

  EMPTY=PoolCollections.GCL(newKeySillys,loci,IDs=NULL,newname=newname)

  newkeygcl=get(paste0(newname,".gcl")) 

  newscores=newkeygcl$scores[,loci,,drop=FALSE]

  newsillyvials=newkeygcl$attributes$SillySource

  rownames(newscores)=newsillyvials

  betweensillys=paste0("between",BetweenSillys)
   
  betweenname=paste(betweensillys,collapse=".")

  EMPTY=PoolCollections.GCL(BetweenSillys,loci,IDs=NULL,newname=betweenname)

  betweengcl=get(paste0(betweenname,".gcl"))

  betweenscores=betweengcl$scores[,loci,,drop=FALSE]

  betweensillyvials=betweengcl$attributes$SillySource

  rownames(betweenscores)=betweensillyvials

  combs=expand.grid(Keysillyvial=newsillyvials,Betweensillyvial=betweensillyvials,stringsAsFactors=FALSE)

  scores=abind(newscores,betweenscores,along=1)

  scores[scores%in%c("0","Unk","XXX")]=NA

  allsillyvials=rownames(scores)

  scores.df=data.frame(t(sapply(allsillyvials,function(id){sapply(loci,function(locus){ifelse(is.na(scores[id,locus,1]),NA,paste(scores[id,locus,1:ploidy[locus]],collapse="."))})})),stringsAsFactors=FALSE)

  duplicaterate=rowMeans(scores.df[combs[,1],]==scores.df[combs[,2],],na.rm=TRUE)

  
  
  duplicatesummary <- sapply(unique(combs$Keysillyvial), function(KeySillyID) {
    
    combs.log <- combs$Keysillyvial == KeySillyID
    
    hist(duplicaterate[combs.log], breaks = seq(0, 1, 0.01), col = 8, xlab = "Duplicate Rate", main = KeySillyID)
    
    abline(v = 0.9, col = "red", lwd = 2)
    
    maxduplicaterate.log <- duplicaterate[combs.log] == max(duplicaterate[combs.log])
    
    maxduplicaterate <- setNames(object = duplicaterate[combs.log][maxduplicaterate.log], nm = combs[combs.log, "Betweensillyvial"][maxduplicaterate.log])
    
    text(x = maxduplicaterate[1], y = 5, labels = paste(names(maxduplicaterate), collapse = "_"), srt = 90, adj = 0)
    
    keymissing <- apply(is.na(scores[combs[combs.log, "Keysillyvial"][maxduplicaterate.log], , 1, drop = FALSE]), 1, sum)
    
    betweenmissing <- apply(is.na(scores[combs[combs.log, "Betweensillyvial"][maxduplicaterate.log], , 1, drop = FALSE]), 1, sum)
    
    maxdups <- cbind(combs[combs.log, ][maxduplicaterate.log, ], Keymissing = keymissing, Betweenmissing = betweenmissing, DuplicateRate = maxduplicaterate, stringsAsFactors=FALSE,row.names=seq(sum(maxduplicaterate.log)))
    
    
    
    thresholdduplicaterate.log <- duplicaterate[combs.log] > threshold
    
    thresholdduplicaterate <- setNames(object = duplicaterate[combs.log][thresholdduplicaterate.log], nm = combs[combs.log, "Betweensillyvial"][thresholdduplicaterate.log])
    
    keymissing <- apply(is.na(scores[combs[combs.log, "Keysillyvial"][thresholdduplicaterate.log], , 1, drop = FALSE]), 1, sum)
    
    betweenmissing <- apply(is.na(scores[combs[combs.log, "Betweensillyvial"][thresholdduplicaterate.log], , 1, drop = FALSE]), 1, sum)
    
    if(sum(thresholdduplicaterate.log) > 0) {
    
      thresholddups <- cbind(combs[combs.log, ][thresholdduplicaterate.log, ], Keymissing = keymissing, Betweenmissing = betweenmissing, DuplicateRate = thresholdduplicaterate, stringsAsFactors=FALSE,row.names=seq(sum(thresholdduplicaterate.log)))
      
    } else{
      
      thresholddups <- paste("No Duplicates at a threshold of", threshold)
      
    }
    
    
    
    if(grep(pattern = "QC", x = combs[1, ]) == 1) {
      projectduplicaterate.log <- combs$Keysillyvial == KeySillyID & combs$Betweensillyvial == gsub(pattern = "QC", replacement = "", x = KeySillyID)
    } else {
      projectduplicaterate.log <- combs$Keysillyvial == KeySillyID & combs$Betweensillyvial == gsub(pattern = "_", replacement = "QC_", x = KeySillyID)
    }
    
    projectduplicaterate <- setNames(object = duplicaterate[projectduplicaterate.log], nm = combs[projectduplicaterate.log, "Betweensillyvial"])

    keymissing <- apply(is.na(scores[combs[projectduplicaterate.log, "Keysillyvial"], , 1, drop = FALSE]), 1, sum)
    
    betweenmissing <- apply(is.na(scores[combs[projectduplicaterate.log, "Betweensillyvial"], , 1, drop = FALSE]), 1, sum)
    
    projectdups <- cbind(combs[projectduplicaterate.log, ], Keymissing = keymissing, Betweenmissing = betweenmissing, DuplicateRate = projectduplicaterate, stringsAsFactors=FALSE,row.names=seq(sum(projectduplicaterate.log)))
    
    
    
    list(Threshold = thresholddups, MostSimilar = maxdups, Project = projectdups)
    
  }, simplify = FALSE)
  
  return(duplicatesummary)
    
}