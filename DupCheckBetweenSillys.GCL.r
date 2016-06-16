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

  
  sapply(unique(combs$Keysillyvial), function(KeysillyID) {
    hist(duplicaterate[combs$Keysillyvial == KeysillyID], breaks = seq(0, 1, 0.01), col = 8, xlab = "Duplicate Rate", main = KeysillyID)
    abline(v = 0.9, col = "red", lwd = 2)
    maxduplicaterate.log <- duplicaterate[combs$Keysillyvial == KeysillyID] == max(duplicaterate[combs$Keysillyvial == KeysillyID])
    maxduplicaterate <- setNames(object = duplicaterate[combs$Keysillyvial == KeysillyID][maxduplicaterate.log], nm = combs[combs$Keysillyvial == KeysillyID, 2][maxduplicaterate.log])
    text(x = maxduplicaterate[1], y = 5, labels = paste(names(maxduplicaterate), collapse = "_"), srt = 90, adj = 0)
    sort(x = duplicaterate[combs$Keysillyvial == KeysillyID], decreasing = TRUE)[1:5]
  }, simplify = FALSE)
  
  
  
  
  
  outliers=outlier(duplicaterate,opposite=FALSE,logical=TRUE)

  keymissing=apply(is.na(scores[combs[outliers,1],,1,drop=FALSE]),1,sum)

  betweenmissing=apply(is.na(scores[combs[outliers,2],,1,drop=FALSE]),1,sum)

  if(sum(outliers)){

    outliers=cbind(combs[outliers,],Keymissing=keymissing,Betweenmissing=betweenmissing,DuplicateRate=duplicaterate[outliers],stringsAsFactors=FALSE,row.names=seq(sum(outliers)))

  }else{

    outliers="No Outliers"

  }
  dupIND=duplicaterate>threshold

  dupcombs=combs[dupIND,]

  duprate=duplicaterate[dupIND]

  keymissing=apply(is.na(scores[dupcombs[,1],loci,1,drop=FALSE]),1,sum)

  betweenmissing=apply(is.na(scores[dupcombs[,2],loci,1,drop=FALSE]),1,sum)

  if(sum(dupIND)){

    Threshold=data.frame(dupcombs,Keymissing=keymissing,Betweenmissing=betweenmissing,DuplicateRate=duprate,stringsAsFactors=FALSE,row.names=seq(sum(dupIND)))

  }else{

    Threshold="No Duplicates"

  }
  return(list(Threshold=Threshold,Outliers=outliers))

}