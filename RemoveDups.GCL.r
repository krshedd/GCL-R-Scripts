RemoveDups.GCL=function(dupcheck){
##################################################################################################################################################################
#
#  Removes duplicated "IDs" with the highest proportion of missingness.
#
#  "dupcheck" is an object created by the function "CheckDupWithinSilly.GCL".
#
##################################################################################################################################################################

  allsillys=names(dupcheck)

  dupcheck=lapply(allsillys,function(silly){dupcheck[[silly]]$report})
  names(dupcheck)=allsillys

  results=vector("list",length(allsillys))
  names(results)=allsillys

  results[1:length(allsillys)]="Nothing Removed"

  sillys=allsillys[sapply(dupcheck,function(dc){is.data.frame(dc)})]

  if(!length(sillys)){return(results)}

  for(silly in sillys){
    IND=apply(dupcheck[[silly]][,c("Missing1","Missing2")],1,function(x){order(as.numeric(x),decreasing=TRUE)[1]})
    remove=as.character(unique(unlist(sapply(seq(nrow(dupcheck[[silly]][,c("ID1","ID2")])),function(i){dupcheck[[silly]][i,c("ID1","ID2")][IND[i]]}))))
    RemoveIDs.GCL(silly=silly,IDs=remove)
    results[[silly]]=remove
  }
  return(results)
}