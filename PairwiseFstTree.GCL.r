PairwiseFstTree.GCL=function(sillyvec,loci,dir,nboots=1000){

  fstatdir=paste(dir,"\\fstatfile.dat",sep="")  

  require("ape")

  require("hierfstat")

  nsillys=length(sillyvec)

  maxsillychar=nchar(nsillys)+1

  nloci=length(loci)

  ploidy=LocusControl$ploidy[loci]

  nalleles=LocusControl$nalleles[loci]

  maxchar=nchar(nalleles)+1
  names(maxchar)=loci

  alleles=LocusControl$alleles[loci]

  my.gcl=lapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)})
  names(my.gcl)=sillyvec

  n=sapply(sillyvec,function(silly){my.gcl[[silly]]$n})
  names(n)=sillyvec

  scores=lapply(sillyvec,function(silly){my.gcl[[silly]]$scores[,loci,]})
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

  dat=read.fstat.data(fstatdir)
  
  Fst=array(0,c(nsillys,nsillys),dimnames=list(sillyvec,sillyvec))

  pairs=combn(sillyvec,2)

  pairnames=apply(pairs,2,function(col){paste(col,collapse=".")})
  
  dimnames(pairs)[[2]]=pairnames

  vc=vector("list",choose(nsillys,2))
  names(vc)=pairnames

  for(pair in pairnames){
    vc[[pair]]=array(NA,c(nloci,3))
    dimnames(vc[[pair]])[[1]]=loci
    sillys=pairs[,pair]
    sillyIND=sapply(sillys,function(silly){match(silly,sillyvec)})
    mydat=dat[dat[,1]==sillyIND[1] | dat[,1]==sillyIND[2],loci]
    myn=n[sillys]
    levels=rep(1:2,myn)
    for(locus in loci){
      mydata=data.frame(levels,mydat[,locus])
      if(ploidy[locus]==1){
        vc[[pair]][locus,1:3]=c(varcomp(data.matrix(mydata),diploid=FALSE)$overall,0)
        if(sum(is.na(vc[[pair]][locus,1:3]))){vc[[pair]][locus,1:3]=c(0,0,0)} 
      }
      if(ploidy[locus]==2){
        vc[[pair]][locus,1:3]=varcomp(data.matrix(mydata),diploid=TRUE)$overall
        if(sum(is.na(vc[[pair]][locus,1:3]))){vc[[pair]][locus,1:3]=c(0,0,0)} 
      }         
    }
    Fst[sillys[1],sillys[2]]=Fst[sillys[2],sillys[1]]=sum(vc[[pair]][,1])/sum(vc[[pair]])       
  }

  tree=nj(Fst)

  trees=bootstrapFst=vector("list",nboots)

  for(boot in 1:nboots){
    tempFst=array(0,c(nsillys,nsillys),dimnames=list(sillyvec,sillyvec))
    temploci=sample(loci,nloci,replace=TRUE)
    for(pair in pairnames){
      sillys=pairs[,pair]
      tempFst[sillys[1],sillys[2]]=tempFst[sillys[2],sillys[1]]=sum(vc[[pair]][temploci,1])/sum(vc[[pair]][temploci,1:3]) 
    }
    bootstrapFst[[boot]]=tempFst
    trees[[boot]]=nj(tempFst)   
  }

  bootstrap=prop.clades(tree,trees)

  return(list(tree=tree,bootstrap=bootstrap,PairwiseFst=Fst,vc=vc,BootstrapFst=bootstrapFst))

}
