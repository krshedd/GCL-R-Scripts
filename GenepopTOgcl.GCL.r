GenepopTOgcl.GCL=function(filename,CreateLocusControl=TRUE,sep="\t"){

  if(strsplit(filename,".gen")==filename){stop("Genepop file needs a '.gen' extention, hoser!!!")}

  if(!CreateLocusControl&!exists("LocusControl")){
    stop("LocusControl must exists if you don't create one, hoser!!!")
  }

  require(hierfstat)

#filename="V:\\DATA\\R_GEN\\DELETEME\\Bill\\ONCORinput.gen"

  rawdat=scan(filename,what="",sep="\n")

  len=length(rawdat)

  popind=as.vector(sapply(rawdat,function(lin){strsplit(lin," ")[[1]][1]=="Pop"}))

  npops=sum(popind)
 
  ORD=order(popind,decreasing=TRUE)[1:npops]

  loci=sapply(rawdat[2:(ORD[1]-1)],function(str){strsplit(str," ")[[1]][1]})

  if(!CreateLocusControl){
    if(sum(is.na(match(loci,LocusControl$locusnames)))){stop("Loci in Genepop file are absent from LocusControl, hoser!!!")}
    loci=LocusControl$locusnames
  }

  nloc=length(loci)

  thisdat=as.vector(sapply(rawdat[-(1:(nloc+1))][rawdat[-(1:(nloc+1))]!="Pop"],function(str){vec=strsplit(str,split=sep)[[1]];vec[ vec!="" & vec!="," ]}))

  thatdat=array(matrix(thisdat,nrow=length(thisdat)/(nloc+1),byrow=TRUE),dim=c(length(thisdat)/(nloc+1),nloc+1),dimnames=list(seq(length(thisdat)/(nloc+1)),c("SillyVial",loci)))

  allpossiblegenotypes=as.vector(sapply(rawdat[(nloc+2):len][!popind[(nloc+2):len]],function(str){vec=as.character(unlist(strsplit(strsplit(str,",")[[1]][2],sep)));vec[vec!=""]}))

  allpossiblegenotypes=allpossiblegenotypes[allpossiblegenotypes!="00"]

  #sum(nchar(allpossiblegenotypes)!=6)

  mynchar=nchar(allpossiblegenotypes)/2

  allpossiblealleles=sapply(1:length(allpossiblegenotypes),function(a){c(substring(allpossiblegenotypes[a],first=1,last=mynchar[a]),substring(allpossiblegenotypes[a],first=mynchar[a]+1,last=2*mynchar[a]))})

  allpossiblealleles=allpossiblealleles[allpossiblealleles!=0] 

  maxallele=max(as.numeric(allpossiblealleles))
  
  maxchar=max(mynchar)

  n=c(ORD[2:npops],len+1)-c(ORD[1:npops])-1

  pops=rep(NA,npops)

  fstat=paste(npops,nloc,maxallele,maxchar,sep=" ")

  fstat=rbind(fstat,cbind(loci))

  for(pop in 1:npops){
    temp=unique(sapply(unique(as.vector(sapply(rawdat[(ORD[pop]+1):(ORD[pop]+n[pop])],function(dat){strsplit(dat,",")[[1]][1]}))),function(str){strsplit(str,"_")[[1]][1]}))
    if(length(temp)>1){temp=strsplit(temp[1],"\\.")[[1]][1]} 
    pops[pop]=temp
    fstat=rbind(fstat,cbind(as.vector(sapply(rawdat[(ORD[pop]+1):(ORD[pop]+n[pop])],function(dat){paste(pop,strsplit(dat,",")[[1]][2],sep=" ")}))))
  }
 
  fstatfile=paste(strsplit(filename,".gen"),".dat",sep="")

  write.table(fstat,fstatfile,row.names=FALSE,col.names=FALSE,quote=FALSE)

  dat=read.fstat.data(fstatfile)

  allelesmat=apply(dat[,2:(nloc+1)],2,function(col){als=unique(as.character(genot2al(col)));als=als[als!=0];als=sort(als);als=als[!is.na(suppressWarnings(as.numeric(als)))];ans=rep(NA,maxchar);ans[1:length(als)]=als;ans})

  alleles=sapply(loci,function(locus){as.character(allelesmat[,locus])},simplify=FALSE)
 
  nalleles=sapply(alleles,function(als){length(als)})

  ploidy=nchar(thatdat[1,2:(nloc+1)])/2
  names(ploidy)=loci

  for(pop in 1:npops){
    mydat=dat[dat[,1]==pop,2:(nloc+1)]
    mypop=pops[pop]   
    scores=array(NA,c(n[pop],nloc,max(ploidy)),dimnames=list(1:n[pop],loci,paste("Dose ",1:max(ploidy),sep="")))
    counts=array(NA,c(n[pop],nloc,max(nalleles)),dimnames=list(1:n[pop],loci,paste("Allele ",1:max(nalleles),sep="")))
    attributes=data.frame(array(NA,c(n[pop],17),dimnames=list(1:n[pop],c("FK_COLLECTION_ID","FK_FISH_ID","PK_TISSUE_TYPE","CAPTURE_LOCATION","CAPTURE_DATE","VIAL_BARCODE","EXHAUSTED_HOW","EXHAUSTED_BY","EXHAUSTED_DATE","STORAGE_ID","UNIT","SHELF_RACK","SLOT","END_CAPTURE_DATE","MESH_SIZE","plateID","SillySource"))))
    attributes$FK_FISH_ID=1:n[pop]
    for(i in 1:n[pop]){
      for(locus in loci){
        if(!length(sort(genot2al(mydat[i,locus]))[sort(genot2al(mydat[i,locus]))!=0])){
          scores[i,locus,1:ploidy[locus]]=rep(NA,ploidy[locus])
        }
        if(length(sort(genot2al(mydat[i,locus]))[sort(genot2al(mydat[i,locus]))!=0])){
          scores[i,locus,1:ploidy[locus]]=sort(genot2al(mydat[i,locus]))[sort(genot2al(mydat[i,locus]))!=0]
        }
        if(sum(is.na(sapply(scores[i,locus,1:ploidy[locus]],function(score){match(score,alleles[[match(locus,loci)]])})))){scores[i,locus,1:ploidy[locus]]=rep(NA,ploidy[locus])}
        for(al in 1:nalleles[match(locus,loci)]){
          counts[i,locus,al]=sum(scores[i,locus,1:ploidy[locus]]==alleles[[match(locus,loci)]][al])
        }
        if(!sum(is.na(counts[i,locus,1:nalleles[match(locus,loci)]]))){
          if(sum(counts[i,locus,1:nalleles[match(locus,loci)]])!=ploidy[locus]){counts[i,locus,1:nalleles[match(locus,loci)]]=rep(NA,nalleles[match(locus,loci)]);scores[i,locus,1:ploidy[locus]]=rep(NA,rep(NA,ploidy[locus]))}
        }
      }
    }
    assign(paste(pops[pop],".gcl",sep=""),list(counts=counts,scores=scores,n=n[pop],attributes=attributes),pos=1)
  }   
 
  if(!exists("LocusControl") & CreateLocusControl){assign("LocusControl",list(MarkerSuite="From GenePop",locusnames=loci,Publishedlocusnames=NULL,alleles=alleles,nalleles=nalleles,ploidy=ploidy),pos=1)}

  return(pops)
}