CheckDupWithinSilly.GCL=function(sillyvec,loci,quantile=0.99,minproportion=0.95){

######################################################################################################################################################################################
#
#  This function checks for duplicate individuals within each silly in "sillyvec".
# 
#  A list is reurned, the same length as sillyvec, each element is a report and distribution of duplication within each silly.
#
#  "sillyvec" is a vector of silly names (without the ".gcl") you wish to check for duplicates.
# 
#  "loci" is a vecor of locus names.
#
#  "quantile" and "minproportion" are used together to determine the cutt-off proportion at which a pair of duplicates is defined: i.e. proportion=max(quantile(duplication,quantile),minproportion)  
#
#  Setting "quantile" equal to NULL will skip the calculation of the duplication distribution and will run much faster--REQUIRES the package "reshape" 
#  
#  Written by JJ 4/24/2011
#  Updated to skip the calculation of the duplication distribution.  JJ 12/19/2011
######################################################################################################################################################################################
 
  while(!require(reshape)){install.packages("reshape")}
  
  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}

  nsilly=length(sillyvec)

  nloci=length(loci)

  ploidy=LocusControl$ploidy[loci]

  resultlist=vector(mode="list",length=nsilly)
  names(resultlist)=sillyvec

  if(is.null(quantile)){
    require("reshape")
    library("dplyr")
    for(silly in sillyvec){
      my.gcl=get(paste(silly,".gcl",sep=""),pos=1)
      scores=my.gcl$scores[,loci,,drop=FALSE]
      ID=as.character(dimnames(scores)[[1]])
      n=length(ID)
      if(n<2){
        resultlist[[silly]]=list(report="No Duplicates",DupDist=NULL) 
        next()
      }
      scores.df=data.frame(t(sapply(ID,function(id){sapply(loci,function(locus){ifelse(is.na(scores[id,locus,1]),NA,paste(scores[id,locus,1:ploidy[locus]],collapse="."))})})))  
      sort.scores.df=reshape::sort_df(scores.df)
      sortIDs=dimnames(sort.scores.df)[[1]]
      duplication=t(sapply(1:(n-1),function(id){ 
                    vec=sapply(1:nloci,function(locus){ 
                          sort.scores.df[id,locus]==sort.scores.df[id+1,locus]
                        });
                    ifelse(sum(!is.na(vec)),sum(vec[!is.na(vec)])/sum(!is.na(vec)),0 )
                    }))
      dupIND=as.vector(duplication>minproportion)
      if(sum(dupIND)){
        dups=t(sapply((1:(n-1))[dupIND],function(id){c(ID1=sortIDs[id],ID2=sortIDs[id+1])}))
        missing=t(sapply(1:nrow(dups),function(dup){
                  vec=match(dups[dup,],sortIDs);
                  c("Missing1"=sum(is.na(sort.scores.df[vec[1],])),"Missing2"=sum(is.na(sort.scores.df[vec[2],])))
                }, simplify = TRUE))
        report = dplyr::bind_cols(dplyr::as_tibble(dups), 
                                  dplyr::as_tibble(missing), 
                                  proportion = duplication[dupIND]
                                  )
      }
      if(!sum(dupIND)){
          report="No Duplicates"
      }
        resultlist[[silly]]=list(report=report,DupDist=NULL)  
    }
    return(resultlist)
  }

  if(!is.null(quantile)){
    for(silly in sillyvec){
      my.gcl=get(paste(silly,".gcl",sep=""),pos=1)
      scores=my.gcl$scores[,loci,]
      mydim=dim(scores)
      if(length(mydim)>2){
        ID=as.character(dimnames(scores)[[1]])
        n=length(ID)
      }
      if(length(mydim)==2){
        ID=as.character(my.gcl$attributes$FK_FISH_ID)
        n=length(ID)
      }
      if(length(mydim)<2){
        ID=NULL
        n=0
      }
     
      if(n<2){resultlist[[silly]]=list(report="No Duplicates",DupDist=NULL) ; next()}


      ncombs=choose(n,2)        

      mycombs=combn(ID,2)

      scores.df=data.frame(t(sapply(ID,function(id){sapply(loci,function(locus){ifelse(is.na(scores[id,locus,1]),NA,paste(scores[id,locus,1:ploidy[locus]],collapse="."))})})))

      
      duplication=sapply(1:ncombs,function(comb){compair=scores.df[mycombs[1,comb],]==scores.df[mycombs[2,comb],];sum(compair[!is.na(compair)])/sum(!is.na(compair))})
 
      names(duplication)=sapply(1:ncombs,function(comb){paste(mycombs[,comb],collapse=".")})
 
      proportion=max(quantile(duplication,quantile),minproportion)
 
      dupIND=duplication>proportion
    
      if(sum(dupIND)){

        dups=data.frame(ID1=mycombs[1,dupIND],ID2=mycombs[2,dupIND])
  
        report=data.frame(dups,t(sapply(seq(nrow(dups)),function(row){df=data.frame(Missing1=sum(is.na(scores.df[as.character(dups[row,1]),]))/nloci,Missing2=sum(is.na(scores.df[as.character(dups[row,2]),]))/nloci);vec=scores.df[as.character(dups[row,1]),]==scores.df[as.character(dups[row,2]),];data.frame(df,proportion=sum(vec[!is.na(vec)])/sum(!is.na(vec)))}) ))
 
      }
      if(!sum(dupIND)){
  
        report="No Duplicates"

      }

        resultlist[[silly]]=list(report=report,DupDist=duplication)  
    }

    detach("package:reshape", unload = TRUE)
    
    return(resultlist)
  }#end if
}