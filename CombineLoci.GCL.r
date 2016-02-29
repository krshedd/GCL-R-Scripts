CombineLoci.GCL=function(sillyvec,markerset,update=FALSE,delim=c(".","_")[1]){

#########################################################################################################################################################################################
#
#  This function combines a set of markers into a single marker.
#
#  "sillyvec" is the set of "*.gcl" objects (without ".gcl" extension) you wish to combine markers for. 
#
#  "markerset" is the set of loci you wish to combine.
#
#  "update" is a logical switch. If TRUE, the "LocusControl" object is updated and all "*.gcl" objects in "sillyvec" will be updated with the new marker.
#                                If FALSE, the "LocusControl" object is not updated and a temporary object called "*.temp.gcl" with the updated data is created.  
#	
#  "delim" specifies the separator between combined loci, either a period (.)which is the default or an underscore (_) so locus names will work in SPAM 4/6/2012 A.B. 
#	
#  Written by JJ 4/24/2011
# 
#########################################################################################################################################################################################

  if(sum(is.na(match(markerset,LocusControl$locusnames)))){stop(paste("'",markerset[is.na(match(markerset,LocusControl$locusnames))],"' from argument 'markerset' not found in 'LocusControl' object!!!",sep=""))}

  nmarkers=length(markerset)  
 
  myploidy=LocusControl$ploidy[markerset]

  if(sort(unique(myploidy))[1]==1){haploid=TRUE}

  if(!sort(unique(myploidy))[1]==1){haploid=FALSE}

  if(sum(myploidy==myploidy[1])!=nmarkers){stop("'markerset' has different ploidies!!!")}  

  MarkerSuite=LocusControl$MarkerSuite
  
  locusnames=LocusControl$locusnames

  newmarkername=paste(markerset,collapse=delim)

  existnewmarker=!is.na(match(newmarkername,locusnames))

  loci=unique(c(locusnames,newmarkername))

  nloci=length(loci)

  Publishedlocusnames=LocusControl$Publishedlocusnames
  Publishedlocusnames=unique(c(Publishedlocusnames,paste(Publishedlocusnames[match(markerset,locusnames)],collapse=".")))

  newalleles=AllPossiblePhenotypes.GCL(markerset)

  alleles=LocusControl$alleles[locusnames]
  nalleles=LocusControl$nalleles[locusnames]
  ploidy=LocusControl$ploidy[locusnames]

  if(!existnewmarker){
    alleles[[length(locusnames)+1]]=newalleles
    nalleles=c(nalleles,length(newalleles))
    ploidy=c(ploidy,1)
  }  
  names(alleles)=names(nalleles)=names(ploidy)=loci

  mymessage=NULL

  for(silly in sillyvec){
    my.gcl=get(paste(silly,".gcl",sep=""),pos=1)

    counts=my.gcl$counts

    scores=my.gcl$scores

    IND=dimnames(scores)

    n=length(IND[[1]])

    if(!is.na(match(newmarkername,IND[[2]]))){
      mymessage=c(mymessage,paste("'",newmarkername,"'"," already created in silly '", silly,"', hoser!!!",sep=""))
      next()
    }

    attributes=my.gcl$attributes
 
    myscores=scores

    ID=as.character(IND[[1]])

    newscores=array(NA,c(n,nloci,max(ploidy)),dimnames=list(ID,loci,1:max(ploidy)))

    newscores[IND[[1]],IND[[2]],as.character(1:max(ploidy))]=scores

    if(!haploid){    
      newscores[IND[[1]],newmarkername,1]=apply(
                                            sapply(1:length(markerset),function(marker){
                                              apply(myscores[,markerset[marker],1:ploidy[markerset[marker]]],1,function(myscore){
                                                if(!sum(is.na(myscore))){
                                                  ans1=paste(myscore,collapse="")
                                                };
                                                if(sum(is.na(myscore))){
                                                  ans1=NA
                                                };
                                                ans1;
                                              });                                                                                 
                                            }),1,function(x){
                                            if(!sum(is.na(x))){
                                              ans2=paste(x,collapse="")
                                            };
                                            if(sum(is.na(x))){
                                              ans2=NA
                                            };
                                            ans2
                                        })
    }
    if(haploid){    
      newscores[IND[[1]],newmarkername,1]=apply( 
                                            sapply(1:length(markerset),function(marker){
                                              sapply(myscores[,markerset[marker],1:ploidy[markerset[marker]]],function(myscore){
                                                if(!sum(is.na(myscore))){
                                                  ans1=paste(myscore,collapse="")
                                                };
                                                if(sum(is.na(myscore))){
                                                  ans1=NA
                                                };
                                                ans1;
                                             })
                                           }),1,function(x){
                                             if(!sum(is.na(x))){
                                               ans2=paste(x,collapse="")
                                             };
                                             if(sum(is.na(x))){
                                               ans2=NA
                                             };
                                             ans2
                                         })
    }


    newcounts=array(NA,c(n,nloci,max(nalleles)),dimnames=list(ID,loci,paste("Allele ",1:max(nalleles),sep="")))
    for(ind in ID){
      for(locus in loci){
        for(allele in 1:nalleles[locus]){
          newcounts[ind,locus,paste("Allele ",allele,sep="")]=sum(newscores[ind,locus,1:ploidy[locus]]==alleles[[locus]][allele]) 
        }
        if(sum(is.na(match(newscores[ind,locus,1:ploidy[locus]],alleles[[locus]])))){
          newscores[ind,locus,1:ploidy[locus]]=rep(NA,ploidy[locus]) 
          newcounts[ind,locus,1:nalleles[locus]]=rep(NA,nalleles[locus]) 
        }
      }
    }    
    if(update){
      assign(paste(silly,".gcl",sep=""),list(counts=newcounts,scores=newscores,n=n,attributes=attributes),pos=1)
    }
    if(!update){
      assign(paste(silly,".temp.gcl",sep=""),list(counts=newcounts,scores=newscores,n=n,attributes=attributes),pos=1)
    } 
  }
  if(update){
    assign("LocusControl",list(MarkerSuite=MarkerSuite,locusnames=loci,Publishedlocusnames=Publishedlocusnames,alleles=alleles,nalleles=nalleles,ploidy=ploidy),pos=1)
  }   
  return(mymessage)
}
