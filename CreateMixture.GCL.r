CreateMixture.GCL=function(sillys,loci,IDs=NULL,mixname,samplevec=NULL,dir,type="BAYES",PT=TRUE,prefix=""){
#############################################################################################################################################################################################################  
#
#  This function writes a mixture to file.
#
#  "sillys" is a vector of sillys, without ".gcl" extension, that the user wants included in the mixture.
# 
#  "loci" is the loci to include.
#
#  "IDs" is a list, same length as "sillys", each element corresponds to a silly and is a vector of IDs (not positions) to include in the mixture.
#
#  "mixname" is name of the mixture file without the ".mix" extension.
#
#  "samplevec" this is only used whey type="HWLER".  This is a vector of numbers the same length as sillys. With baseline pops getting numbered 1:npops and mixtures getting zeros.
#  
#  "dir" is the directory where the user wants the file written to.
#
#  "type" is either "BAYES", "HWLER", "BUGS" or "SPAM".
#
#  "PT" whether this mixture is for a prooftest or not.  If not set to FALSE for a normal mixture the original ".gcl" object will be removed. #This needs a better fix! Andy Barclay 10/10/2011
#
#  "prefix" only used type = "BUGS"
#
# Written by Jim Jasper 6/9/2011
#############################################################################################################################################################################################################  
  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}

  if(is.na(match(type,c("BAYES","HWLER","BUGS","SPAM")))){stop(paste("What is type ",type,", hoser?",sep=""))}

  if(is.null(IDs)){
    IDs=lapply(sillys,function(silly){
          mygcl=get(paste(silly,".gcl",sep=""),pos=1);
          counts=mygcl$counts;
          attributes=mygcl$attributes
          mydim=dim(counts)[dim(counts)>0];
          if(length(mydim)>2){
            ids=as.character(dimnames(counts)[[1]]) 
          };
          if(length(mydim)==2){
            ids=as.character(attributes$FK_FISH_ID)
          };
           if(length(mydim)<2){
            ids=NULL
          };
          ids
        })
  }

  if(!is.list(IDs)){stop("'IDs' must be a list")}

  IDs=lapply(IDs,function(ID){
        if(is.null(ID)){
          ids=NULL
        };
        if(!is.null(ID)){
          ids=as.character(ID)
        };
        ids
      })
  names(IDs)=sillys

  filename=paste(dir,"\\",mixname,".mix",sep="")

  nloci=length(loci)

  nalleles=LocusControl$nalleles[loci]
  
  if(!PT==TRUE){newname=substr(paste(paste(sillys,collapse="."),"new",sep=""), start=1, stop=200)}
  
  if(PT==TRUE){newname=substr(paste(sillys,collapse="."), start=1, stop=200)}

  PoolCollections.GCL(collections=sillys,loci,IDs=IDs,newname=newname)

  counts=get(paste(newname,".gcl",sep=""),pos=1)$counts

  counts[is.na(counts)]=0

  ID=as.character(dimnames(counts)[[1]]) 

  if(type=="BUGS"){  
    if(!exists("writeDatafileR")){source("V:\\DATA\\R_GEN\\nonGCL Source Scripts\\writeDatafileR.r")}
    BUGSdir=dir
    if(!file.exists(BUGSdir)){dir.create(BUGSdir);print(paste("Created '",BUGSdir,"', hoser!!!",sep=""))}    
    x=counts[,loci,]
    data=list(x=x,M=nrow(x))
    writeDatafileR(DATA=data,towhere=paste(BUGSdir,"\\",prefix,"BUGSdata_x.txt",sep=""))
    fortranformat=NULL
  }

  if(type=="BAYES"){
    mix=sapply(ID,function(id){paste(unlist(lapply(loci,function(locus){paste(counts[id,locus,1:nalleles[locus]],collapse="")})),collapse=" ")}) 
    mix=cbind(rep("",length(mix)),mix)
    write.table(mix,filename,row.names=FALSE,col.names=FALSE,sep=" ",quote=FALSE)
    rle<-rle(nalleles[loci])
    fortranformat=paste("(", paste(sapply(seq(length(rle$values)),function(locus){paste(rle$lengths[locus],"(1X,",rle$values[locus],"I1)",sep='')} ),collapse=","),")",sep='')
  }
  
   if(type=="HWLER"){
    maxchar=max(nchar(samplevec))
    sample=unlist(sapply(1:length(sillys),function(n){rep(samplevec[n],length(IDs[[n]]))}))
    sample=format(sample,length=maxchar)
    SillyID=as.vector(unlist(sapply(1:length(sillys),function(i){paste(sillys[i],IDs[[i]],sep="_")})))
    mix=sapply(ID,function(id){paste(unlist(lapply(loci,function(locus){paste(counts[id,locus,1:nalleles[locus]],collapse="")})),collapse=" ")}) 
    mix=cbind(sample,mix,SillyID)
    write.table(mix,filename,row.names=FALSE,col.names=FALSE,sep=" ",quote=FALSE)
    rle<-rle(nalleles[loci])
    fortranformat=paste("(","I",maxchar,",", paste(sapply(seq(length(rle$values)),function(locus){paste(rle$lengths[locus],"(1X,",rle$values[locus],"I1)",sep='')} ),collapse=","),",1X,A",max(nchar(SillyID)),")",sep='')
  }

  if(type=="SPAM"){
    mix=sapply(1:length(ID),function(id){paste(c("    ",unlist(lapply(loci,function(locus){paste(counts[ID[id],locus,1:nalleles[locus]],collapse="")})),id),collapse="  ")})
    write.table("* characters",filename,row.names=FALSE,col.names=FALSE,sep="\n",quote=FALSE) 
    write.table(cbind(1:nloci,loci),filename,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE) 
    write.table("* end",filename,row.names=FALSE,col.names=FALSE,sep="\n",quote=FALSE,append=TRUE) 
    write.table("\\",filename,row.names=FALSE,col.names=FALSE,sep="\n",quote=FALSE,append=TRUE) 
    write.table(mix,filename,row.names=FALSE,col.names=FALSE,sep="\n",quote=FALSE,append=TRUE)
    fortranformat=NULL
  }

  suppressWarnings(rm(list=paste(newname,".gcl",sep="")))

  return(fortranformat)  
}