RepeatedProofTest.GCL=function(sillyvec,groupnames,loci,groupvec,samplesizemat,prprtnl=TRUE,type="BAYES",testdir,nrepeats,nreps,burn,nchains=1,RunBUGS=FALSE){


    if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}

    G=length(groupnames)

    if(G!=ncol(samplesizemat)){
      stop(paste("Wrong length for 'groupnames', hoser!!!",sep=""))
    }

    npoints=nrow(samplesizemat)

    priorvec=Prior.GCL(groupvec,rep(1/G,G))

    initmat=array(NA,c(length(sillyvec),nchains),dimnames=list(sillyvec,paste("Chain",1:nchains,sep="")))
    for(chain in 1:nchains){
      for(g in 1:G){
        IND=groupvec==g
        gg=rgamma(sum(IND),1/G,1) 
        initmat[sillyvec[IND],chain]=gg 
      }
      initmat[,chain]=initmat[,chain]/sum(initmat[,chain])   
    }

  if(type=="BUGS"){
    require("BRugs")

    model="V:\\DATA\\R_GEN\\BUGS Scripts\\SimpleGSI.txt"

    resultlist=vector("list",nrepeats)
    names(resultlist)=paste("Repeat",1:nrepeats,sep="") 

    for(peat in 1:nrepeats){
      myrepeat=paste("Repeat",peat,sep="")
      resultlist[[myrepeat]]=vector("list",npoints)
      names(resultlist[[myrepeat]])=paste("Point",1:npoints,sep="")
      for(point in 1:npoints){
        mypoint=paste("Point",point,sep="")
        myprefix=paste(myrepeat,mypoint,"BUGS",sep="")
        mydir=paste(testdir,"\\",myprefix,sep="")
        if(!file.exists(mydir)){dir.create(mydir)}
        ProofTest.GCL(sillyvec=sillyvec,loci=loci,groupnames=groupnames,groupvec=groupvec,samplesize=samplesizemat[point,],prefix=myprefix,dir=mydir,prprtnl=TRUE,type="BUGS",priorvec=priorvec)
        if(RunBUGS==TRUE){
          modelCheck(model)
          modelData(paste(mydir,"\\BUGS\\",myprefix,"BUGSdata.txt",sep=""))
          modelCompile(numChains=nchains) 
          modelGenInits()
          modelUpdate(burn)
          samplesSet("R")
          modelUpdate(nreps)
          R=samplesStats("R")
          dimnames(R)[[1]]=groupnames
          if(nchains>1){R=data.frame(R,GR=cbind(sapply(samplesBgr("R",plot=F),function(lst){lst$bgrRatio})[50,]))}
          resultlist[[myrepeat]][[mypoint]]=R
        }
      }
      return(resultlist)
    }
  }
  if(type=="BAYES"){
    for(peat in 1:nrepeats){
      myrepeat=paste("Repeat",peat,sep="")
      for(point in 1:npoints){
        mypoint=paste("Point",point,sep="")
        myprefix=paste(myrepeat,mypoint,"BAYES",sep="")
        mydir=paste(testdir,"\\",myprefix,sep="")
        if(!file.exists(mydir)){dir.create(mydir)}
        ProofTest.GCL(sillyvec=sillyvec,loci=loci,groupnames=groupnames,groupvec=groupvec,samplesize=samplesizemat[point,],prefix=myprefix,dir=mydir,type="BAYES",nreps=nreps,nchains=nchains,priorvec=priorvec,initmat=initmat)
      }
    }
    return(NULL)       
  }
}
