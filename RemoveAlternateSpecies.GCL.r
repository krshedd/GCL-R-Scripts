RemoveAlternateSpecies.GCL=function(AlternateSpeciesReport,AlternateCutOff,FailedCutOff){
  IND=AlternateSpeciesReport$Alternate>AlternateCutOff & AlternateSpeciesReport$Failure>FailedCutOff
  if(sum(IND)){
    SillyVialsToRemove=dimnames(AlternateSpeciesReport[IND,])[[1]]
    SillyVialsToRemoveMat=matrix(unlist(strsplit(SillyVialsToRemove,split="_")),byrow=TRUE,nrow=length(SillyVialsToRemove))
    MySillys=unique(SillyVialsToRemoveMat[,1])
    IDs=vector("list",length(MySillys))
    names(IDs)=MySillys
    for(silly in MySillys){
      IDs[[silly]]=as.character(SillyVialsToRemoveMat[SillyVialsToRemoveMat[,1]==silly,2])
      RemoveIDs.GCL(silly=silly,IDs=IDs[[silly]])
    }
  }
  return(NULL)
}
