FishersTest.GCL <- function(freq, loci, tests){

#  tests=list(sillyvec[1:2],sillyvec[2:3],sillyvec[c(1,3)],sillyvec[1:3])

  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}
  
  testnames=sapply(tests,function(test){paste(test,collapse=".")})

  fullresults=lapply(tests,function(test){FisherCompute.GCL(freq[test,,],loci)})

  names(fullresults)=testnames

  overall=sapply(fullresults,function(x){x[x[,1]=="Overall",2]})

  return(list(OverallResults=cbind(overall),ResultsByLocus=fullresults))
}