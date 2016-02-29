CreateSPAMcontrol.GCL <- function(sillyvec,groupnames,mixvec,groupvec,baselinefile,loci,dir,outdir=dir,mixdir=dir,basedir=dir,simulate=FALSE,
                                   options=c("off","off","off","off","off","off","off","off","off","off","off","on","on","off"),
                                   ctlPar=c("1000","520","200","2000","19",".1E-9",".1E-15",".1E-45",".1E-9","400",".1E-9","95")){

#This function produces a SPAM control file for each mixture in "mixvec" or a simulation control file for each group in "groupnames".  
############Below are examples of the required arguments ############################################
#mixvec<-c("TestMixture1","TestMixture2","TestMixture3","TestMixture4")
#sillyvec<-c("Quartz","Crescent","UpperKen","Juneau","Russian","Benjamin","Killey","Funny","LowerKen","Slikok")
#groupnames<-c("UpperTrib","LowerTrib","Mainstem")
#groupvec<-c(1,1,3,3,1,2,2,2,3,2) #same length at sillyvec
#baseline<-"Kenai10Pops41loci"
#loci<-loci40
#dir<-"V:/WORK/Andy/R/SPAMctl" # This is the location where the ctl files will be written. 
############Optional Arguments with defaults #########################################################
#outdir -the directory where the SPAM results will be written, if not specified it defaults to "dir"
#mixdir -the directory where the SPAM mixture file is located, if not specified it defaults to "dir"
#basedir -the directory where the SPAM baseline file is located, if not specified it defaults to "dir"
#options -character vector turning "on" or "off" the 14 Optimization options in the following order:
#          1) use IRLS algorithm                      
#          2) fixed baseline frequencies              
#          3) print mixture file                      
#          4) print baseline relative frequencies     
#          5) print conditional genotype prob.        
#          6) print conditional population prob.      
#          7) print iteration history                 
#          8) print bootstrap estimates               
#          9) compute likelihood confidence intervals 
#          10) compute infinitesimal jacknife std. dev.
#          11) compute likelihood ratio                
#          12) resample using mixture frequencies      
#          13) resample using baseline frequencies     
#          14) use pella-masuda model                  
#ctlPar - character vector of control parameters in the following order:
#          1) number of bootstrap resamplings 
#          2) maximum number of genotypes     
#          3) maximum number of classes       
#          4) maximum number of iterations    
#          5) maximum number of missing loci  
#          6) estimate tolerance              
#          7) likelihood tolerance            
#          8) genotype tolerance              
#          9) algorithm tolerance             
#          10) simulation sample size          
#          11) guaranteed percent of maximum                
#########################################################################################################
# written by A.B 9/28/2011
#########################################################################################################

 if(simulate){names=groupnames} 
 if(!simulate){names<-mixvec}

  for(mix in names) {

    if(simulate){ctl<-paste("* simulate: ",mix,sep='')}
    if(!simulate){ctl<-paste("* Estimate: ",mix,sep='')}
    
    #Options Block
    ctl<-rbind(ctl,"* options selected for optimization")
    ctl<-rbind(ctl,paste("  use IRLS algorithm                      :  ",options[1],sep=''))
    ctl<-rbind(ctl,paste("  fixed baseline frequencies              :  ",options[2],sep=''))
    ctl<-rbind(ctl,paste("  print mixture file                      :  ",options[3],sep=''))
    ctl<-rbind(ctl,paste("  print baseline relative frequencies     :  ",options[4],sep=''))
    ctl<-rbind(ctl,paste("  print conditional genotype prob.        :  ",options[5],sep=''))
    ctl<-rbind(ctl,paste("  print conditional population prob.      :  ",options[6],sep=''))
    ctl<-rbind(ctl,paste("  print iteration history                 :  ",options[7],sep=''))
    ctl<-rbind(ctl,paste("  print bootstrap estimates               :  ",options[8],sep=''))
    ctl<-rbind(ctl,paste("  compute likelihood confidence intervals :  ",options[9],sep=''))
    ctl<-rbind(ctl,paste("  compute infinitesimal jacknife std. dev.:  ",options[10],sep=''))
    ctl<-rbind(ctl,paste("  compute likelihood ratio                :  ",options[11],sep=''))
    ctl<-rbind(ctl,paste("  resample using mixture frequencies      :  ",options[12],sep=''))
    ctl<-rbind(ctl,paste("  resample using baseline frequencies     :  ",options[13],sep=''))
    ctl<-rbind(ctl,paste("  use pella-masuda model                  :  ",options[14],sep=''))
  
    #Control Parameters
    ctl<-rbind(ctl,"* control parameters")
    ctl<-rbind(ctl,paste("  number of stocks in analysis     : ",length(sillyvec),sep=''))
    ctl<-rbind(ctl,paste("  number of characters in analysis : ",length(loci),sep=''))
    ctl<-rbind(ctl,paste("  number of bootstrap resamplings  : ",ctlPar[1],sep=''))
    ctl<-rbind(ctl,paste("  maximum number of genotypes      : ",ctlPar[2],sep=''))
    ctl<-rbind(ctl,paste("  maximum number of classes        : ",ctlPar[3],sep=''))
    ctl<-rbind(ctl,paste("  maximum number of iterations     : ",ctlPar[4],sep=''))
    ctl<-rbind(ctl,paste("  maximum number of missing loci   : ",ctlPar[5],sep=''))
    ctl<-rbind(ctl,paste("  estimate tolerance               : ",ctlPar[6],sep=''))
    ctl<-rbind(ctl,paste("  likelihood tolerance             : ",ctlPar[7],sep=''))
    ctl<-rbind(ctl,paste("  genotype tolerance               : ",ctlPar[8],sep=''))
    ctl<-rbind(ctl,paste("  algorithm tolerance              : ",ctlPar[9],sep=''))
    if(simulate){ctl<-rbind(ctl,paste("  simulation sample size           : ",ctlPar[10],sep=''))}
    ctl<-rbind(ctl,paste("  guaranteed percent of maximum    : ",ctlPar[11],sep=''))
    
    randseeds<-sample(21474836,2)
    ctl<-rbind(ctl,paste("  random seed                      : -",randseeds[1],sep=''))
    ctl<-rbind(ctl,paste("  second random seed               : ",randseeds[2],sep=''))

    #Character Block
    ctl<-rbind(ctl,"* characters used")
    ctl<-rbind(ctl,"[id#]  [kind]  [#types/alleles] [character]")
    ploidy<-gsub("1","mtdna",LocusControl$ploidy[loci])
    ploidy<-gsub("2","locus",ploidy)
    ctl<-rbind(ctl,cbind(sapply(1:length(loci),function(d){paste(sprintf("%4d",d),"  ",ploidy[d],format(as.character(LocusControl$nalleles[loci][d]),width=17,justify="centre"),loci[d])})))

    #Population Block
    maxchar<-max(nchar(sillyvec))+2
    ctl<-rbind(ctl,"* populations used in analysis")
    ctl<-rbind(ctl,paste("[id#]",format("[population]",width=maxchar),"      [lev1]"))
    ctl<-rbind(ctl,cbind(sapply(1:length(sillyvec),function(s){paste(sprintf("%3d",s),"  ", format(sillyvec[s],width=maxchar)," : ",sprintf("%6d",groupvec[s]))})))
 
    #Regions block
    ctl<-rbind(ctl,"* regions")
  
    if(simulate){maxchar<-max(nchar(groupnames))+2
                 ctl<-rbind(ctl,paste("[level]","  ","[label]","  ","[region]","               ","[estimate]",sep=''))
                 }
                 
    if(!simulate){ctl<-rbind(ctl,paste("[level]","  ","[label]","  ","[region]",sep=''))}
                  
    if(simulate){ctl<-rbind(ctl,cbind(sapply(1:length(groupnames),function(g){paste(format("1",width=7,justify="centre"),format(as.character(g),width=9,justify="centre"),format(groupnames[g],jjustify="left",width=maxchar),"  :       ",match(groupnames[g],mix,nomatch=0))})))}
                  
    if(!simulate){ctl<-rbind(ctl,cbind(sapply(1:length(groupnames),function(g){paste(format("1",width=7,justify="centre"),format(as.character(g),width=9,justify="centre"),format(groupnames[g],jjustify="left"))})))}
           
    #Files
    ctl<-rbind(ctl,"* files")
    ctl<-rbind(ctl,paste(" baseline: ",basedir,"/",baselinefile,".bse",sep=''))
    if(!simulate){ctl<-rbind(ctl,paste(" mixture : ",mixdir,"/",mix,".mix",sep=''))}
    ctl<-rbind(ctl,paste(" output  : ",outdir,"/",mix,sep=''))   
    ctl<-rbind(ctl,"* run") 

    write.table(ctl,paste(dir,"/",mix,".ctl",sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE)                         
  }
  return(NULL)
}
