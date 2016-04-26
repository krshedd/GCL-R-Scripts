PlotLikeProfile.GCL=function(likeprof,popvec,loci,groupvec,groupnames,dir,filename,col=NULL){
#
#** Example 1 -- Kenai Chinook Baseline ********************************************************************************************************************* 
#
#  rm(list=ls(all=TRUE))
#
#  load("V:/WORK/Chinook/Cook Inlet/Kenai River/2011 Baseline/KenaiChinook2011.RData")
#
#  popvec=KenaiNames10;loci=loci40;groupvec=groupvec3;groupnames=groups;dir="V:/WORK/Chinook/Cook Inlet/Kenai River/2011 Baseline"
#
#  likeprof=LeaveOneOutLikeProfile.GCL(popvec,loci,groupvec,groupnames,groupcomps=NULL,ncores=15)
#
#  PlotLikeProfile.GCL(likeprof,popvec,loci,groupvec,groupnames,dir)
#
#** Example 2 -- GAPS Chinook Baseline ********************************************************************************************************************* 
#
#  rm(list=ls(all=TRUE))
#
#  load("V:\\DATA\\WORK\\CHINOOK\\SE Alaska\\GAPS3.0\\GAPSbaseline_Feb2013.RData")
#
#  popvec=SEAKPops356;groupvec=GroupVec25RG;groupnames=Groups25
#
#  likeprof=LeaveOneOutLikeProfile.GCL(popvec=SEAKPops356,loci=GAPSloci,groupvec=GroupVec25RG,groupnames=Groups25,groupcomps=NULL,ncores=16)
#    
#  PlotLikeProfile.GCL(likeprof,popvec,loci,groupvec,groupnames,dir)
#
#

  names(groupvec)=popvec

  x=factor(likeprof$Attributes$FromPop,levels=popvec,ordered=FALSE)

  xlab="Assignment Threshold"

  ylab="Threshold Achievement Fraction"

  if(is.null(col)) {col=setNames(sample(rainbow(max(groupvec)),max(groupvec)),groupnames)}

  xx=as.character(seq(0,100,5)/100)

  xlim=c(1/max(groupvec),1)

  pdf(paste0(dir,"/",filename,"PopToGroupAssignmentThresholdCurves.pdf"), width = 9)

  for(group in groupnames){

    main=paste("Assignment Into",group)

    y=round(2*likeprof$IndividualByGroup[,group],1)/2

    summarytab=array(0,c(length(popvec),length(xx)),list(popvec,xx))

    summarytab0=table(data.frame(x,y))[popvec,]

    summarytab[popvec,colnames(summarytab0)]=summarytab0

    summarytab=summarytab/apply(summarytab,1,sum)

    summarytab=apply(summarytab,1,rev)

    summarytab=apply(summarytab,2,cumsum)

    summarytab=apply(summarytab,2,rev)

    yy=summarytab[,1]

    plot(x=xx,y=yy,type="l",lwd=2,col=col[groupvec[1]],main=main,xlab=xlab,ylab=ylab,xaxp=c(0,1,10),yaxp=c(0,1,10),font=2,font.axis=2,font.lab=2, ylim = c(0,1), xlim = c(0, 1.5))#,xlim=xlim

    for(pop in popvec[-1]){

      yy=summarytab[,pop]

      lines(x=xx,y=yy,lwd=2,col=col[groupvec[pop]])


    }#pop

    legend(x="topright",legend=paste("From",groupnames),fill=col,bty="n",text.font=2)

  }#group

  dev.off()

  return(NULL)
  
}
