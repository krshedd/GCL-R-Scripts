CreateBaseline.GCL=function(sillyvec,loci,dir,basename,type="BAYES",groupvec=NULL){
##########################################################################################################################################################################################################
#  
#      if 'type'="BUGS" the hoser must supply groupvec.
#
#
#
#
#
#
##########################################################################################################################################################################################################



  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}

  if(is.na(match(type,c("BAYES","BUGS","SPAM")))){stop(paste("What is type ",type,", hoser?",sep=""))}

  filename=paste(dir,"\\",basename,".bse",sep="")

  nsilly=length(sillyvec)

  nloci=length(loci)

  nalleles=LocusControl$nalleles[loci]

  ploidy=LocusControl$ploidy[loci]

  y=FreqPop.GCL(sillyvec,loci)

  n=sapply(loci,function(locus){apply(y[sillyvec,locus,1:nalleles[locus]],1,sum)})
  
  dimnames(n)=list(sillyvec,loci)

  maxchar=max(nchar(n))

  if(type=="BUGS"){  
    if(!exists("writeDatafileR")){source("V:\\DATA\\R_GEN\\nonGCL Source Scripts\\writeDatafileR.r")}
    BUGSdir=dir
    if(!file.exists(BUGSdir)){dir.create(BUGSdir);print(paste("Created '",BUGSdir,"', hoser!!!",sep=""))}    
    ORD=order(groupvec)
    Cg=as.vector(table(groupvec[ORD]))
    cmCg=c(0,cumsum(Cg))
    data=list(C=length(sillyvec),cmCg=cmCg,G=max(groupvec),D=length(loci),Jd=as.numeric(nalleles),ploidy=as.numeric(ploidy),y=y[ORD,,],n=n[ORD,],alpha=rep(1/Cg,Cg))
    writeDatafileR(DATA=data,towhere=paste(BUGSdir,"\\",basename,"BUGSdata.txt",sep=""))
    fortranformat=NULL
  }

  if(type=="BAYES"){  

    baseline=array(paste(c(rep(" ",maxchar-1),0),collapse=""),c(nsilly*nloci,max(nalleles)+3),dimnames=list(1:(nsilly*nloci),c("Pop","Locus","N",paste("Allele ",1:max(nalleles),sep=""))))

    SillyLocusGrid=expand.grid(dimnames(y)[[1]],dimnames(y)[[2]])
  
    sillyIND=rep(1:nsilly,each=nloci) 

    locusIND=rep(1:nloci,times=nsilly)
 
    baseline[,1:2]=cbind(sapply(sillyIND,function(ind){paste(c(rep(" ",maxchar-nchar(ind)+1),ind),collapse="")}),sapply(locusIND,function(ind){paste(c(rep(" ",maxchar-nchar(ind)),ind),collapse="")}))
 
    for(i in 1:nsilly){
      silly=sillyvec[i]
      for(d in 1:nloci){
        locus=loci[d]
        baseline[sillyIND==i & locusIND==d,3:(3+nalleles[locus])]=c(paste(c(rep(" ",maxchar-nchar(n[silly,locus])),n[silly,locus]),collapse=""),sapply(1:nalleles[locus],function(allele){paste(c(rep(" ",maxchar-nchar(y[silly,locus,allele])),y[silly,locus,allele]),collapse="")}))  #Changed silly and locus to i and d
      }
    }                                                                                           

    write.table(baseline,filename,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

    fortranformat=paste("(",max(nalleles)+3,"(1X,I",maxchar,"))",sep="")  

  }

  if(type=="SPAM"){  
   
    baseline=y   
  
    baseline[is.na(baseline)]=""   

    emtylist=sapply(1:nsilly,function(i){write.table(paste("#  ",i,"    ",sillyvec[i],sep=""),filename,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE,append=i>1);write.table(baseline[sillyvec[i],,],filename,sep=" ",row.names=TRUE,col.names=FALSE,quote=FALSE,append=TRUE)})

    fortranformat=NULL

  }
  
  return(fortranformat)

}