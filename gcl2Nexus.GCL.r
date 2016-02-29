gcl2Nexus.GCL=function(sillyvec,loci,path,VialNums=TRUE,PopNames=sillyvec){

###################################################################################################################################################################################
#This function writes a Nexus format file for use in GDA.  
#Haplotype loci are moved to the end columns so the "hapset" is a continuous range. Discontinuous hapsets don't work in GDA even though the documentation says it's possible.
########################################################################### Arguments ############################################################################################# 
#sillyvec - character vector of ".gcl" object names without the .gcl extension.  
#loci - character vector of locus names. Loci with a ploidy of "1" will be resorted to the end of the list. 
#path - directory where the file will be written, include the file name with the ".nex" extension.
#VialNumbers - logical statement, if set to FALSE, no vial numbers will be written.  Default is TRUE
#PopNames - a character vector the same length as sillyvec to give populations new names. If no vector is given, PopNames defaults to "sillyvec".  
###################################################################################################################################################################################  
#This function is a modification of the gcl2Genepop.GCL function written by Jim Jasper.  A.B. 11/4/2011
###################################################################################################################################################################################
                         
  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}
  
 
  ploidy=sort(LocusControl$ploidy[loci],decreasing=TRUE)
  
  loci=names(ploidy)
  
  alleles=LocusControl$alleles[loci]
 
  nalleles=LocusControl$nalleles[loci]
 
  maxchar= max(nchar(unlist(sapply(sillyvec,function(silly){as.character(get(paste(silly,".gcl",sep=''))$attributes$SillySource)}))))

  npops=length(sillyvec)
  
  nloci=length(loci)
  
  if(!length(PopNames)==length(sillyvec)){stop("PopNames is not the same length as sillyvec.")}
  
  names(PopNames)=sillyvec
  
  file="#nexus"
  
  file=rbind(file,"")
  
  file=rbind(file,"begin gdadata; [!GDA input format]")
  
  file=rbind(file,paste("dimensions npops=",npops,"  nloci=",nloci,";",sep=''))
  
  file=rbind(file,paste("format missing=? separator=/;",sep=''))
  
  if(min(ploidy)==1){
  hapset=grep("1",ploidy)
  if(length(hapset)>1){hapset=paste(range(hapset),collapse="-")}
  file=rbind(file,paste("hapset ",hapset,";",sep=''))
  }
  
  file=rbind(file,"locusallelelabels")
  
  locuslabels=cbind(c(paste("   ",format(seq(nloci-1),width=max(nchar(seq(length(loci))))),"  ",loci[seq(nloci-1)],",",sep=""),paste("   ",format(nloci,width=max(nchar(seq(length(loci))))),"  ",loci[nloci],sep="")))
  file=rbind(file,locuslabels)
  
  file=rbind(file,cbind(c(";","matrix")))

  for(silly in sillyvec){    
    my.gcl=get(paste(silly,".gcl",sep=""),pos=1)
    vials=format(as.character(my.gcl$attributes$SillySource),width=maxchar+1,justify="right")
    if(!VialNums){vials=sapply(vials,function(vial){strsplit(vial,"_")[[1]][1]})}
    IDs=dimnames(my.gcl$scores)[[1]]
    names(vials)=IDs    
    scores=my.gcl$scores[,loci,]
    counts=sapply(IDs,function(ID){paste(sapply(loci,function(locus){paste(sapply(1:ploidy[locus],function(ploid){ifelse(is.na(scores[ID,locus,ploid]),paste("?",collapse=""),paste(match(scores[ID,locus,ploid],alleles[[locus]]),collapse=""))}),collapse="/  ")}),collapse="   ")})     
    counts=as.character(sapply(IDs,function(ID){paste(paste(vials[ID]," ",sep=""),counts[ID],collapse="")}))
    
    file=rbind(file,cbind(c(paste(PopNames[silly],":",sep=''),counts)))
    file=rbind(file,cbind(c(",")))
    
  }
  file[length(file):(length(file)+1)]=c(";","end;")
  
  write.table(file,path,quote=FALSE,row.names=FALSE,col.names=FALSE)
  
  return(NULL) 
 
}