HoFisFstTable.GCL=function(sillyvec,loci,fstatdir = NULL, dir = NULL){
  
  require("hierfstat")
  
  nsillys=length(sillyvec)
  
  maxsillychar=nchar(nsillys)+1
  
  nloci=length(loci)
  
  ploidy=LocusControl$ploidy[loci]
  
  nalleles=LocusControl$nalleles[loci]
  
  maxchar=nchar(nalleles)+1
  names(maxchar)=loci
  
  alleles=LocusControl$alleles[loci]
  
  if(is.null(fstatdir)) {
    
    fstatdir=paste(dir,"\\fstatfile.dat",sep="")  
    
    my.gcl=lapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)})
    names(my.gcl)=sillyvec
    
    n=sapply(sillyvec,function(silly){my.gcl[[silly]]$n})
    names(n)=sillyvec
    
    scores=lapply(sillyvec,function(silly){scrs=my.gcl[[silly]]$scores[,loci,];scrs[scrs%in%c("0","Unk","XXX")]=NA;scrs})
    names(scores)=sillyvec
    
    counts=lapply(sillyvec,function(silly){
      sapply(1:n[silly],function(i){
        paste(c(match(silly,sillyvec),sapply(loci,function(locus){
          ifelse(is.na(scores[[silly]][i,locus,1]),paste(rep(0,ploidy[locus]*maxchar[locus]),collapse=""),paste(sapply(1:ploidy[locus],function(allele){
            paste(c(rep(0,maxchar[locus]-nchar(match(scores[[silly]][i,locus,allele],alleles[[locus]]))),match(scores[[silly]][i,locus,allele],alleles[[locus]])),collapse="")
          }),collapse=""))
        })),collapse=" ")
      })
    })      
    names(counts)=sillyvec
    
    fstat=paste(nsillys,nloci,max(nalleles),max(maxchar),sep=" ")
    
    fstat=rbind(fstat,cbind(loci))
    
    fstat=rbind(fstat,cbind(as.vector(unlist(counts))))
    
    write.table(fstat,fstatdir,row.names=FALSE,col.names=FALSE,quote=FALSE) 
  
  }

  dat=read.fstat.data(fstatdir)

  MyVC=array(0,c(nloci,3),dimnames=list(loci,c("P","I","G")))
  
  for(locus in loci){
    diploid=ploidy[locus]==2
    if(diploid){
      MyVC[locus,c("P","I","G")]=hierfstat::varcomp(dat[,c("Pop",locus)],diploid=diploid)$overall
    }
    if(!diploid){
      MyVC[locus,c("P","G")]=hierfstat::varcomp(dat[,c("Pop",locus)],diploid=diploid)$overall
    }    
  }

  MyTable=MyVC[loci,c("P","I")]/apply(MyVC,1,sum)

  MyTable=rbind(MyTable,Overall=apply(MyVC[loci,c("P","I")],2,sum)/sum(apply(MyVC,1,sum)))

  Ho=rep("-",nloci+1)

  names(Ho)=c(loci,"Overall")

  Ho[loci[ploidy==2]]=apply(basic.stats(dat[,c("Pop",loci[ploidy==2])])$Ho,1,mean)

  Ho["Overall"]=mean(as.numeric(Ho[loci[ploidy==2]]))

  MyTable=data.frame(Ho=Ho,Fis=MyTable[,"I"],Fst=MyTable[,"P"])

  return(MyTable)
}
