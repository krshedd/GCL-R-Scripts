ReadBiomarkQC.GCL=function(QCcsvFilepaths,skip=15){
###################################################################################################################################################################################################################################################
#
#QCcsvFilepaths=c("V:/DATA/All SNP data/Chinook/Project K42 Chinook Inseason_Postseason 2012/Project K42 M192 BioMark/QC/Data/Single Chip/K42 Combined Output.csv","V:/DATA/All SNP data/Chinook/Project K45- AYK, BB Chinook baseline/Project K45 192chip BioMark/QC/Data/Single Chip/K45_QC2_1381716192/K45_QC2_ModifiedforR.csv")
#
#
#
###################################################################################################################################################################################################################################################

  sillyvec=as.vector(sapply(objects(pattern="*\\.gcl",pos=1),function(gclname){strsplit(gclname,split="\\.gcl")[[1]][1]}))

  GenotypesQC=NULL

  for(pth in QCcsvFilepaths){
    GenotypesQC=rbind(GenotypesQC,read.table(pth,header=TRUE,sep=",",colClasses="character",stringsAsFactors=FALSE,skip=skip)[,c("Name","Converted","Assay","ID","Allele.X","Allele.Y")])
  }
  
  GenotypesQC=GenotypesQC[GenotypesQC$Name!="NTC",] #Added the removal of NTCs 10/2/2012 A.B.
  
  GenotypesQC$Converted[GenotypesQC$Converted=="No Call"]=NA
  
  GenotypesQC$Converted[GenotypesQC$Converted=="Invalid"]=NA

  MyTab=table(GenotypesQC$Name,GenotypesQC$Assay)

  SillyVialQC=dimnames(MyTab)[[1]]

  nVials=length(SillyVialQC)
  
  lociQC=dimnames(MyTab)[[2]]

  nlociQC=length(lociQC)

  loci=LocusControl$locusnames

  nalleles=LocusControl$nalleles[loci]

  ploidy=LocusControl$ploidy[loci]

  nloci=length(loci)

  alleles=LocusControl$alleles[loci]

  names(loci)=loci

  IND=match(lociQC,loci)

  lociQCmissng=lociQC[is.na(IND)]

  if(sum(is.na(IND))){

    print(paste("QC loci: ",paste("'",lociQCmissng,"'",collapse=",")," do not exsist in 'LocusControl'!!! : Carrying On...!!!",sep=""))

    lociQC=lociQC[!is.na(IND)]

    nlociQC=length(lociQC)
    
  }


  Scores=array(NA,c(nVials,nloci),dimnames=list(SillyVialQC,loci))

  for(vial in SillyVialQC){
    
    for(locus in loci){
     
      Scores[vial,locus]=GenotypesQC[GenotypesQC$Name==vial & GenotypesQC$Assay==locus,"Converted"][1]
   
    }
  
  }
  

  GenotypesQCsillys=as.character(sapply(SillyVialQC,function(nm){strsplit(nm,split="_")[[1]][1]}))

  sillyvecQC0=unique(GenotypesQCsillys)

  sillyvecQC=sillyvecQC0[!is.na(match(sillyvecQC0,sillyvec))]

  GenotypesQCIDs=as.character(sapply(SillyVialQC,function(nm){strsplit(nm,split="_")[[1]][2]}))

  AttributsColNames=dimnames(get(objects(pattern="*\\.gcl",pos=1)[1],pos=1)$attributes)[[2]]

  for(silly in sillyvecQC){
    IND=GenotypesQCsillys==silly
    IDs=GenotypesQCIDs[IND]
    n=length(IDs)
    scores=array(NA,c(n,nloci,max(ploidy)),dimnames=list(IDs,loci,paste("Dose",1:max(ploidy),sep="")))
    counts=array(NA,c(n,nloci,max(nalleles)),dimnames=list(IDs,loci,paste("Allele",1:max(nalleles),sep="")))
    for(id in IDs){
      idsillyvial=paste(silly,id,sep="_")
      for(locus in loci){
        if(!is.na(Scores[idsillyvial,locus])){
          scores[id,locus,1:ploidy[locus]]=sort(strsplit(Scores[idsillyvial,locus],split=":")[[1]][1:ploidy[locus]])
          for(alll in 1:nalleles[locus]){
            counts[id,locus,alll]=sum(scores[id,locus,1:ploidy[locus]]==alleles[[locus]][alll])
          }
        }
      }
    }
    attributes=data.frame(array(NA,c(n,length(AttributsColNames)),dimnames=list(IDs,AttributsColNames)))
    attributes$FK_FISH_ID=IDs
    attributes$SillySource=as.character(sapply(SillyVialQC[IND],function(sllyvl){paste(strsplit(sllyvl,split="_")[[1]],collapse="QC_")}))
    assign(paste(silly,"QC.gcl",sep=""),list(counts=counts,scores=scores,n=n,attributes=attributes),pos=1)
  }
assign(x = "QCSillys", value = paste(sillyvecQC,"QC",sep=""), pos = 1)
# return(paste(sillyvecQC,"QC",sep=""))
}






                         