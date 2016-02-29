AllPossiblePhenotypes.GCL=function(markerset){

    alleles=LocusControl$alleles[markerset]

    ploidy=LocusControl$ploidy[markerset]

    newgntps=lapply(lapply(markerset,function(marker){
                      t(combn(rep(alleles[[marker]],ploidy[marker]),ploidy[marker]))
                    }),function(gmat){
                                      na=NULL
                                      nr=nrow(gmat)
                                      nc=ncol(gmat)
                                      for(i in 1:nr){
                                        if(nc>1){vec=sort(gmat[i,])}
                                        else{vec=gmat[i]} 
                                        na=c(na,paste(vec,collapse=""))
                                      }
                                      sort(unique(na))
             })

     if(length(unique(markerset))==1){

       newalleles=newgntps[[1]]

     }else{
      
       temp=expand.grid(newgntps,stringsAsFactors=FALSE)     

       nr=nrow(temp)

       na=NULL

       for(i in 1:nr){

         na=c(na,paste(as.character(temp[i,]),collapse=""))

       }      

       newalleles=sort(unique(na))

    }#else 

    return(newalleles)
}