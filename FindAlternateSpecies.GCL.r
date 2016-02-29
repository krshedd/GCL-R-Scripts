FindAlternateSpecies.GCL=function(sillyvec,species="chum"){

  require("abind")

  if(species=="chum"){
    AlternateGenotypesPath="V:\\DATA\\R_GEN\\Alternate Homo Markers Chum.txt"
    FailedMarkersPath="V:\\DATA\\R_GEN\\Failed Chum Markers.txt"
  }

  if(species=="sockeye"){
    AlternateGenotypesPath="V:\\DATA\\R_GEN\\Alternate Homo Markers Sockeye.txt"
    FailedMarkersPath="V:\\DATA\\R_GEN\\Failed Sockeye Markers.txt"
  }


  AlternateGenotypes=read.table(AlternateGenotypesPath,sep="\t",header=TRUE)

  AlternateMarkers=as.vector(AlternateGenotypes[,1])

  loci=LocusControl$locusnames

  if(sum(is.na(match(AlternateMarkers,loci)))){
    stop("'AlternateMarkers' not found in 'loci', hoser!!!") 
  }

  FailedMarkers=scan(FailedMarkersPath,what="")[-1]

  if(sum(is.na(match(FailedMarkers,loci)))){
    stop("'FailedMarkers' not found in 'loci', hoser!!!") 
  }


  gclobjects=sapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)$scores[,c(AlternateMarkers,FailedMarkers),]},simplify=FALSE)

  mysillyvial=sapply(sillyvec,function(silly){paste(silly,dimnames(gclobjects[[silly]])[[1]],sep="_")},simplify=FALSE)

  sillyvial=as.vector(unlist(mysillyvial))

  gclobjectAll=NULL
  for( silly in sillyvec){
    gclobjectAll=abind(gclobjectAll,gclobjects[[silly]],along=1)
  }
  dimnames(gclobjectAll)[[1]]=sillyvial


  Failure=apply(is.na(gclobjectAll[sillyvial,FailedMarkers,1]),1,sum)/length(FailedMarkers)

  ploidy=LocusControl$ploidy[c(AlternateMarkers,FailedMarkers)]

  gclobjectAllAlternate=t(sapply(sillyvial,function(vial){
                            sapply(AlternateMarkers,function(loc){
                              ifelse(is.na(gclobjectAll[vial,loc,1]),NA,paste(gclobjectAll[vial,loc,1:ploidy[loc]],collapse=""))
                            })
                          }))

   Alternate=apply(
               sapply(AlternateMarkers,function(loc){
                 gclobjectAllAlternate[sillyvial,loc]==AlternateGenotypes[AlternateGenotypes$AlternateMarker==loc,"AlternateGenotype"]
               }),1,function(rw){
                   ind=!is.na(rw);
                   ifelse(length(rw[ind]),sum(rw[ind])/length(rw[ind]),0)  
               })
   
   Results=data.frame(Failure=Failure,Alternate=Alternate)

   MyResultsTab=table(Results)

   xx=as.numeric(dimnames(MyResultsTab)[[1]])

   yy=as.numeric(dimnames(MyResultsTab)[[2]])

   xy=expand.grid(xx,yy)

   plot(x=xy[,1],y=xy[,2],xlim=c(0,1),ylim=c(0,1),col="white",xlab="Failed",ylab="Alternate",xaxp=c(0,1,10),yaxp=c(0,1,10))

   Axis(x=c(0,1),side=4,yaxp=c(0,1,10))

   Axis(x=c(0,1),side=3,xaxp=c(0,1,10))

   labels=paste(xx,sep="")

  IND=labels>"0"

  for(i in 1:nrow(MyResultsTab)){
    text(x=rep(xx[i],length(yy)),y=yy,labels=MyResultsTab[i,],col="blue",font=2,cex=0.75)
  }

 
  return(Results)

}




