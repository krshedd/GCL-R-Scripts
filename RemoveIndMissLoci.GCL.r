RemoveIndMissLoci.GCL=function(sillyvec,proportion=0.8){

######################################################################################################################################################################################
#
#  This function removes individuals from "*.gcl" objects that have fewer non-missing loci than that specified by "proportion".  
#  
#  A list is returned, the same length as "sillyvec", each element is a vector containing the ID (not position) of each fish removed. 
#
#  "sillyvec" is a vector of silly names (without the ".gcl") you wish to remove individuals from.
# 
#  "proportion" is the cut-off proportion of the number of non-missing loci.  
#
#  ReWritten by JJ 12/17/2015
#  Update by Kyle Shedd 1/18/2019
#
######################################################################################################################################################################################


  results=setNames(vector("list",length(sillyvec)),sillyvec)

  for(silly in sillyvec){

    my.gcl=get(paste0(silly,".gcl"),pos=1)

    counts=my.gcl$counts

    ID=as.character(rownames(counts))

    nloci=sum(apply(counts, 2, function(locus) {!all(is.na(locus))} ))  # exclude SNPs not run for a given silly

    if(my.gcl$n > 1) {
      myproportion=apply(!is.na(counts[,,1]),1,sum)/nloci
    } 
    
    if(my.gcl$n == 1){
      myproportion=setNames(object = sum(!is.na(counts[,,1]))/nloci, nm = dimnames(counts)[[1]])
    }
    
    myind=myproportion>=proportion
 
    counts=counts[myind,,,drop=FALSE]      

    scores=my.gcl$scores[myind,,,drop=FALSE]

    n=nrow(scores)

    attributes=my.gcl$attributes[myind,,drop=FALSE]

    assign(paste0(silly,".gcl"),list(counts=counts,scores=scores,n=n,attributes=attributes),pos=1)

    if(sum(!myind)==0){

      results[[silly]]="None"

    } else{

      results[[silly]]=ID[!myind]

    }
                                                                                                                                                            
  }#silly

  return(results)

}
