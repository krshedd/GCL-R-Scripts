gcl2Genepop.GCL=function(sillyvec,loci,path,VialNums=TRUE){

  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}

  alleles=LocusControl$alleles[loci]
 
  ploidy=LocusControl$ploidy[loci]
 
  nalleles=LocusControl$nalleles[loci]
 
  maxchar=nchar(nalleles)+1

  file="GENEPOP input format"
 
  file=rbind(file,cbind(loci))

  for(silly in sillyvec){    
    my.gcl=get(paste(silly,".gcl",sep=""),pos=1)
    vials=as.character(my.gcl$attributes$SillySource)
    if(is.na(vials)[1]==TRUE){vials=paste(silly,"_",dimnames(my.gcl$scores)[[1]],sep='')}
    if(!VialNums){vials=sapply(vials,function(vial){strsplit(vial,"_")[[1]][1]})}
    IDs=dimnames(my.gcl$scores)[[1]]
    names(vials)=IDs    
    scores=my.gcl$scores[,loci,]
    counts=sapply(IDs,function(ID){paste(sapply(loci,function(locus){paste(sapply(1:ploidy[locus],function(ploid){ifelse(is.na(scores[ID,locus,ploid]),paste(rep(0,maxchar[locus]),collapse=""),paste(c(rep(0,maxchar[locus]-nchar(match(scores[ID,locus,ploid],alleles[[locus]]))),match(scores[ID,locus,ploid],alleles[[locus]])),collapse=""))}),collapse="")}),collapse=" ")})
    counts=as.character(sapply(IDs,function(ID){paste(paste(vials[ID]," , ",sep=""),counts[ID],collapse=" ")}))
    file=rbind(file,cbind(c("Pop",counts)))
  }

  write.table(file,path,quote=FALSE,row.names=FALSE,col.names=FALSE)
  
  return(NULL) 
 
}