CreateControlFile.GCL<-function(sillyvec,loci,mixname,basename,suffix="",nreps=40000,nchains,groupvec,priorvec,initmat,dir,seeds=matrix(sample(seq(10000),3*nchains),nrow=3),
                               thin=c(1,1,1),mixfortran,basefortran,switches="F T F T F T F"){


  chains=paste("Chain",1:nchains,sep="")

  dirs=paste(dir,"\\",mixname,suffix,"Chain",1:nchains,".ctl",sep="")
  names(dirs)=chains

  priorvec=cbind(substring(format(round(priorvec,6),nsmall=6),first=2))

  initmat=cbind(substring(format(round(initmat,6),nsmall=6),first=2))
  if(nchains>1){
    dimnames(initmat)[[2]]=chains
  }

  nsillys=length(sillyvec)

  ploidy=LocusControl$ploidy[loci]

  nalleles=LocusControl$nalleles[loci]

  nloci=length(loci)
  
  seeds=cbind(seeds)
  dimnames(seeds)=list(1:3,chains)

  files=lapply(chains,function(chain){paste(mixname,suffix,chain,sep="")})
  names(files)=chains

  files=lapply(files,function(file){rbind(file,paste(basename,".bse",sep=""))})
  names(files)=chains

  files=lapply(files,function(file){rbind(file,paste(mixname,".mix",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(mixname,suffix,chain,"SUM.SUM",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(mixname,suffix,chain,"BOT.BOT",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(mixname,suffix,chain,"FRQ.FRQ",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(mixname,suffix,chain,"BO1.BO1",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(mixname,suffix,chain,"CLS.CLS",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(mixname,suffix,chain,"RGN.RGN",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],nreps)})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],length(sillyvec))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],length(loci))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],cbind(seeds[,chain]))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],cbind(thin))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],mixfortran)})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],basefortran)})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],switches)})
  names(files)=chains


  files=lapply(seq(length(chains)),function(chain){rbind(files[[chain]],cbind(sapply(1:nloci,function(d){paste(sprintf("%3s",d),sprintf("%2s",nalleles[loci[d]]),sprintf("%2s",ifelse(ploidy[loci[d]]==2,"T","F")),"",loci[d],collapse="")})))})
  names(files)=chains


  
  files=lapply(chains,function(chain){rbind(files[[chain]],cbind(sapply(1:nsillys,function(i){paste(sprintf("%3s",i),sprintf("%2s",groupvec[i]),sprintf("%7s",priorvec[i]),"", format(ifelse(nchar(sillyvec[i])<18,sillyvec[i],substr(sillyvec[i],start=1,stop=18)),width=18),sprintf("%7s",initmat[i,chain]),collapse="")})))})
  names(files)=chains

  empty=sapply(chains,function(chain){write.table(files[[chain]],dirs[chain],quote=FALSE,row.names=FALSE,col.names=FALSE)})

  return(NULL)


}
