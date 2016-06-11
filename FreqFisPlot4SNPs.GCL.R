FreqFisPlot4SNPs.GCL=function(sillyvec,loci,groupvec,alpha=0.05,groupcol=NULL,file=NULL,group.pch=19,dot.cex=1,pval.cex=1,pval.digits=2){

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function will create a pdf file with plots of allele frequency and Fis 
  # for each locus on a separate page. HWE p-values are printed if less than
  # alpha.
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # sillyvec = a vector of sillys, without ".gcl" extension that you want the 
  #   freqencies plotted for
  # loci = a vector of loci that you want included in the frequency plots, 
  #   Note: Only use loci with ploidy=2
  # groupvec = a vector of group numbers corresponding to each silly in sillyvec
  # alpha = critical HWE p-value
  # groupcol = optional vector of colors corresponding to each group with 
  #   length = max(groupvec), can be either color numbers (numeric vector) or
  #   color names (character vector). If left NULL (default) a rainbow of colors
  #   will be used.
  # file = the full file path, with .pdf extension where the file will be written
  #   If no file is supplied, the default is to write the file "FreqPlot.pdf" to 
  #   the current working directory.
  # group.pch = ector of pch values corresponding to each group with length = 
  #   max(groupvec), default is 19
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # A pdf file with allele frequency plots for each locus. Also plots Fis for
  # each locus and indicates the HWE pvalue if less than alpha.
  #
  # Also outputs a list object with three matrices with silly as row and loci as columns
  #   1) allele.freqs
  #   2) Fis
  #   3) HWE.pval
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   FreqFisPlot_SWKodiakCollections46GroupVec8 <- FreqFisPlot4SNPs.GCL(
  #     sillyvec = SWKodiakCollections46, loci = loci96, 
  #     groupvec = SWKodiakCollections46GroupVec8, alpha = 0.05, 
  #     groupcol = c("cyan", "red", "blue", "blue4", "green", "purple", "grey", "orange"), 
  #     file = "HWE/FreqFisPlot_SWKodiakCollections46GroupVec8.pdf")
  #    
  # Created by Jim Jasper November 2015
  # Updated by Kyle Shedd Fri Jun 10 15:49:37 2016 in order to generalize and better document
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(is.null(file)){

    file=paste0(getwd(),"/FreqFisPlot.pdf")

  }

  require(HardyWeinberg)

  C=length(sillyvec)

  ploidy=LocusControl$ploidy[loci]

  nalleles=LocusControl$nalleles[loci]

  loci=loci[ploidy==2 & nalleles==2]

  nloci=length(loci)

  fis=array(NA,c(C,nloci),list(sillyvec,loci))

  pval=array(NA,c(C,nloci),list(sillyvec,loci))

  for(silly in sillyvec){

    gclName=paste0(silly,".gcl")

    gcl=get(gclName,pos=1)

    genecounts0=3-gcl$counts[,loci,1]

    genecounts=t(apply(genecounts0,2,tabulate,nbins=3))

    HWEpval=suppressWarnings(HWExactMat(genecounts))$pvalvec

    pval[silly,loci]=HWEpval
    
    FixBOOL=apply(genecounts,1,function(gcounts){sum(gcounts==0)>1})

    HWEfis=suppressWarnings(apply(genecounts[!FixBOOL,],1,HWf))

    fis[silly,!FixBOOL]=HWEfis
 
    fis[silly,FixBOOL]=0
    
  }#silly

  G=max(groupvec)

  if(is.null(groupcol)){
    groupcol=rainbow(G)
  }
    
  if(is.numeric(groupcol)){
     groupcol=colors()[groupcol]
  }  

  PopCol=groupcol[groupvec] 
 
  Freq=FreqPop.GCL(sillyvec=sillyvec,loci=loci)

  n=apply(Freq,1:2,sum)
  
  q=Freq[,,"Allele 1"]/n
  
  if(length(group.pch)==1){

    group.pch=rep(group.pch,G)

  }
  
  PopPch=group.pch[groupvec]
    
  fisylim0=max(abs(fis))

  fisylim=c(-fisylim0,fisylim0)

  pdf(file,width=11,height=8,family="Helvetica",pointsize=20)

  for(locus in loci){
  
    layout(mat=cbind(1:2))

    par(mar=c(0,4,2,0)+0.1)

    plot(q[,locus],main=locus,col=PopCol,pch=PopPch,ylim=c(0,1),ylab="Freqency",xaxt="n",cex=dot.cex)
    
    lines(q[,locus],lty="dotted")

    par(mar=c(4,4,0,0)+0.1)

    plot(fis[,locus],pch=PopPch,col=PopCol,ylim=fisylim,ylab="Fis",xlab="Population",cex=dot.cex)

    lines(fis[,locus],lty="dotted")

    labels=as.character(round(pval[,locus],pval.digits))

    labels[labels=="0"]=paste(c("<0.",rep(0,pval.digits),"5"),collapse="")

    labels[pval[,locus]>alpha]=NA

    pos=rep(3,C)

    pos[fis[,locus]==fisylim0]=1 

    pos[fis[,locus]<0]=1   

    pos[fis[,locus]==-fisylim0]=3 

    text(x=seq(C),y=fis[,locus],labels=labels,pos=pos,cex=pval.cex,font=2,col=PopCol)
    
    abline(h=0)

  }

  dev.off()

  return(list(allele.freqs=q,Fis=fis,HWE.pval=pval))
 
}

  
  