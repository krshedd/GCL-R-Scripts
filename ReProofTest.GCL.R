ReProofTest.GCL=function(sillyvec,loci,groupvec,groupnames,ProofTestIDs,prefix="",dir,suffix="",nreps=40000,nchains=5,priorvec,initmat,type="BAYES",thin=c(1,1,1),switches="F T F T F T F"){
############################################################################################################################################################################################################################################
#
#  sillyvec=BristolBay146Pops
#
#  loci=loci91
#
#  groupvec=BristolBay146Pop12GroupVec
#
#  groupnames=c(BristolBay11Groups[1:4],"NushNoNuya","NuyaNoNush",BristolBay11Groups[6:11])
#
#  ProofTestIDs=NushNoNuyaProof
#
#  prefix="NushNoNuya"
#
#  dir="V:/Work/Sockeye/Bristol Bay/2013 Baseline/For JJ/BAYES"
#
#  priorvec=BristolBay146Pop12FlatPrior
#
#  initmat=BristolBay146PopsInits
#
#  
#
#  
## Changed by Kyle Shedd 1/4/2016: added "switches" and "thin" argument to allow changing CreateControlFile.GCL within this function and allow for Individual Assignment, changes made to line 144 to reflect
############################################################################################################################################################################################################################################

  if(!file.exists(dir)){

    dir.create(dir)

  }  

  IncludedGroups=names(ProofTestIDs)

  baselinesillyvec=sillyvec

  names(baselinesillyvec)=sillyvec

  mixsillyvec=NULL

  baselineSillysIncluded=NULL

  for(group in IncludedGroups){

    g=match(group,groupnames)

    groupsillys=sillyvec[groupvec==g]

    newnames=paste(group,groupsillys,sep=".")

    names(newnames)=groupsillys

    names(ProofTestIDs[[group]])=groupsillys

    IND=sapply(groupsillys,function(silly){!is.null(ProofTestIDs[[group]][[silly]])})

    newnames=newnames[IND]

    for(silly in groupsillys[IND]){

      newname=newnames[silly]

      gclname=paste(silly,"gcl",sep=".")

      silly.gcl=get(gclname,pos=1)

      baselineSillysIncluded=c(newname,baselineSillysIncluded)

      baselinesillyvec[silly]=newname

      newgclname=paste(newname,"gcl",sep=".")

      newmixname=paste("mix",newname,sep=".")

      mixsillyvec=c(newmixname,mixsillyvec)

      newmixgclname=paste(newmixname,"gcl",sep=".")

      assign(newgclname,silly.gcl,pos=1)

      assign(newmixgclname,silly.gcl,pos=1)

      AllIDs=dimnames(silly.gcl$scores)[[1]]

      baselineIDs2remove=as.character(sort(as.numeric(ProofTestIDs[[group]][[silly]])))

      mixIDs2remove=AllIDs[is.na(match(AllIDs,baselineIDs2remove))]

      RemoveIDs.GCL(silly=newname,IDs=list(baselineIDs2remove))

      RemoveIDs.GCL(silly=newmixname,IDs=list(mixIDs2remove))
  
    }

  }

  mixdir <- paste(dir, "/", type, ".mixture", sep = "")  # mixdir=paste(dir,"Mixture",sep="/")

  if(!file.exists(mixdir)){

    dir.create(mixdir)

  }

  basedir <- paste(dir, "/", type, ".baseline", sep = "")  # basedir=paste(dir,"Baseline",sep="/")
  
  if(!file.exists(basedir)){
    
    dir.create(basedir)
    
  }

  controldir <- paste(dir, "/", type, ".control", sep = "")  # controldir=paste(dir,"Control",sep="/")
  
  if(!file.exists(controldir)){
    
    dir.create(controldir)
    
  }  

  if(type%in%c("Genepop","genepop","GENEPOP")){

    gcl2Genepop.GCL(sillyvec=mixsillyvec,loci,path=paste0(mixdir,"/",prefix,"Mixture.gen"),VialNums=FALSE)

    gcl2Genepop.GCL(sillyvec=baselinesillyvec,loci,path=paste0(basedir,"/",prefix,"Baseline.gen"),VialNums=FALSE)

  }  

  if(type%in%"BAYES"){

    MixFortran=CreateMixture.GCL(mixsillyvec,loci,mixname=prefix,dir=mixdir,type="BAYES",PT=TRUE,prefix="")

    BaseFortran=CreateBaseline.GCL(baselinesillyvec,loci,dir=basedir,basename=prefix,type="BAYES",groupvec=NULL)

    CreateControlFile.GCL(baselinesillyvec,loci,mixname=prefix,basename=prefix,suffix="",nreps,nchains,groupvec,priorvec,initmat,controldir,seeds=matrix(sample(seq(10000),3*nchains),nrow=3),thin=thin,mixfortran=MixFortran,basefortran=BaseFortran,switches=switches)    

  }#if

  rm(silly.gcl)

  rm(list=paste(c(mixsillyvec,baselineSillysIncluded),"gcl",sep="."),pos=1)

  return(NULL)

}
