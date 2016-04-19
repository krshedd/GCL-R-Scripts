ReadGenepopDis.GCL<-function(file){
#################################################################################################################
#
# This function reads in the population output from a GENEPOP disequilibrium ("*.dis") file and returns a data.frame.
# P-values with "No contingency table" are replaced with "1" and p-values = "0" are replaced with 1/(#Batches*#Iterations per batch).
# 
# Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# file - "V:/WORK/Sockeye/Kodiak/Kodiak Afognak Baseline Look 2014/Kyle/Genepop/Kodiak50Pops92NuclearLoci.txt.DIS"
#   ~ the full file path, including the ".dis" extension.  Make sure the file has not been modified.
# 
# Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Data.frame
# nrow = number of pairwise combinations of loci (i.e. choose(n = nloci, k = 2)).
# ncol = number of populations + 2
#   ~ The first two columns are the locus names for the pair (Locus1 and Locus2)
#   ~ Subsequent columns contain the LD p-values for that locus pair for a given population
#
# Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How many populations have a p-value of < 0.05 for a given locus pair?
# Create an additional column with the number of pops below a given p-value threshold (i.e. 0.05)
# LD <- ReadGenepopDis.GCL(file = "V:/Analysis/3_AYK/Chum/WASC/Genepop/WASCchum.DIS")
# LD$npopsfail <- apply(LD[, 3:dim(LD)[2]] < 0.05, 1, function(locuspair) {sum(locuspair)} )
#
# Written 2/5/13 Andy Barclay
# updated 8/6/14 Kyle Shedd - data.frame with p-values for each pop + overall
#################################################################################################################

  require("reshape2")

  dis<-scan(file,what='',sep = "\n")

  batches<-as.numeric(strsplit(dis[grep("	Batches              :",dis)],split="\tBatches              : ")[[1]][2])

  iterations<-as.numeric(strsplit(dis[grep("	Iterations per batch : ",dis)],split="\tIterations per batch : ")[[1]][2])

  repzero<-format(1/(batches*iterations),scientific=F,digits=6)

  popstart<-grep("Pop             Locus#1  Locus#2    P-Value      S.E.     Switches",dis)+2

  popend<-grep("P-value for each locus pair across all populations",dis)-1

  poptable<-sapply(popstart:popend,function(row){gsub(pattern="[[:blank:]]+",x=dis[row],replacement="/",fixed=F)})

  poptable<-colsplit(poptable,pattern="/",names=c("Pop","Locus#1","Locus#2","PValue","S.E.","Switches"))

  popscol<-gsub(poptable$Pop,pattern="_\\d+",replacement='')

  poptable$Pop=popscol

  npops<-as.numeric(strsplit(dis[grep("Number of populations detected : ",dis)],split="Number of populations detected : ")[[1]][2])

  ncomps<-length(popscol)/npops

  pops=popscol[c(ncomps*1:npops)]

  if(sum(duplicated(x = pops)) > 0) {pops[duplicated(x = pops)] = paste(pops[duplicated(x = pops)], "dup", sep = '')}

  poptable$PValue[poptable$PValue=="No"]=1
  
  poptable$PValue[as.numeric(poptable$PValue)==0]=repzero

  locstart<-popend+6

  locend<-grep("Normal ending.",dis)-1

  loctable<-sapply(locstart:locend,function(row){gsub(pattern="[[:blank:]]+",x=dis[row],replacement="/",fixed=F)})

  loctable<-gsub("/&/",loctable,replacement="/")

  loctable<-colsplit(loctable,pattern="/",names=c("Locus#1","Locus#2","Chi2","df","PValue"))

  loctable$PValue<-gsub(loctable$PValue,pattern="Highly/sign.",replacement='0')

  summary=data.frame(Locus1=loctable[,"Locus#1"], Locus2=loctable[,"Locus#2"], stringsAsFactors = FALSE)

  for (i in 1:npops){
    summary[,paste(pops[i])]=as.numeric(poptable[(((i-1)*ncomps)+1):(i*ncomps),"PValue"])
  }

  summary[,"Overall"]=as.numeric(loctable[,"PValue"])

  return(summary)

}