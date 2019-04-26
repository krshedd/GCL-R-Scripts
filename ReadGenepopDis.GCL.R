ReadGenepopDis.GCL <- function(file, loci = NULL) {
#################################################################################################################
#
# This function reads in the output from a GENEPOP disequillibrium ("*.dis") file and returns a tibble with p-value columns for each population
# and a column of overall population p-values. 
#
# P-values with "No contingency table" are replaced with "1" and p-values = "0" are replaced with 1/(#Batches*#Iterations per batch).
# 
# Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# file <- "V:/Analysis/2_Central/Coho/Cook Inlet/2018/Baseline/Genepop/CI94pops82loci.DIS"
#   ~ the full file path, including the ".dis" extension.  Make sure the file has not been modified.
#
# loci <- NULL
#   ~ defult is NULL or a vector of locus names the same order as the genepop input file used for the tests. 
# If a vector of locus names is supplied the Locus1 and Locus2 columns will contain the correct locus names
# otherwise the locus names will be the same as in the .DIS file (tuncated to 8 characters)
# 
# Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# tibble
# nrow = number of pairwise combinations of loci (i.e. choose(n = nloci, k = 2)).
# ncol = number of populations + 3
#   ~ The first two columns are the locus names for the pair (Locus1 and Locus2)
#   ~ Subsequent columns contain the LD p-values for that locus pair for a given population
#   ~ The last column contains the overall p-values for all populations. 
#     If only one population is tested the overall p-value column will contain NA's
#
# Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How many populations have a p-value of < 0.05 for a given locus pair?
# Create an additional column with the number of pops below a given p-value threshold (i.e. 0.05)
# LD <- ReadGenepopDis.GCL(file = "V:/Analysis/3_AYK/Chum/WASC/Genepop/WASCchum.DIS",loci=NULL)
# LD$npopsfail <- apply(LD[, 3:dim(LD)[2]] < 0.05, 1, function(locuspair) {sum(locuspair)} )
#
# Written 4/25/19 Andy Barclay 
#################################################################################################################

  require(tidyverse)

  dis <- scan(file,what='',sep = "\n")
  
  npops <- as.numeric(strsplit(dis[grep("Number of populations detected :",dis)],split="Number of populations detected :")[[1]][2])
  
  nloci <- as.numeric(strsplit(dis[grep("Number of loci detected        :",dis)],split="Number of loci detected        :")[[1]][2])

  batches <- as.numeric(strsplit(dis[grep("	Batches              :",dis)],split="\tBatches              : ")[[1]][2])

  iterations <- as.numeric(strsplit(dis[grep("	Iterations per batch : ",dis)],split="\tIterations per batch : ")[[1]][2])

  repzero <- format(1/(batches*iterations),scientific=F,digits=6)

  popstart <- grep("Pop             Locus#1  Locus#2    P-Value      S.E.     Switches",dis)+2
  
  ncomps<-choose(nloci,2)
  
  popend <- popstart+ncomps*npops-1
  
  #Loci
  if(!is.null(loci)) {
    
    loc1 <- sapply(seq(length(loci)-1),function(i){
      
      loci[1:i]
    
    }) %>% unlist()
    
    loc2 <- sapply(seq(length(loci)),function(i){
      
      rep(loci[i],i-1)
    
    }) %>% unlist()
    
    pop_df0 <- separate(data = tibble(dat = dis[popstart:popend]), col = dat,sep ="[[:blank:]]+",into = c("PopG","Locus1","Locus2","PValue",NA,NA), remove = TRUE) %>% 
      mutate(Pop = gsub(PopG, pattern ="_\\d+", replacement = ''),`Locus1`= loc1,`Locus2`= loc2, PValue = gsub(PValue, pattern = "No", replacement = "1")) %>% 
      mutate(PValue = gsub(PValue,pattern = "$$0", replacement = repzero) %>% as.numeric())
    
    
  } else {
    
    pop_df0 <- separate(data = tibble(dat = dis[popstart:popend]), col = dat,sep ="[[:blank:]]+",into = c("PopG","Locus1","Locus2","PValue",NA,NA), remove = TRUE) %>% 
      mutate(Pop = gsub(PopG, pattern ="_\\d+", replacement = ''), PValue = gsub(PValue, pattern = "No", replacement = "1")) %>% 
      mutate(PValue = gsub(PValue,pattern = "$$0", replacement = repzero) %>% as.numeric())
    
  } #end Loci
  
  pop_names <- pop_df0 %>% pull(Pop) %>% unique() %>% sort()
  
  pop_df <- pop_df0 %>%
    select(-Pop) %>% 
    group_by_at(vars(-PValue)) %>%
    mutate(row_id = 1:n()) %>% 
    ungroup() %>%
    spread(key = PopG, value = PValue) %>% 
    select(-row_id)
  
  #npops
  if(npops==1) {
    
      summary_df <- pop_df%>% 
        mutate(Overall = NA) %>% 
        set_names("Locus1","Locus2",pop_names,"Overall")
  
  } else {
    
    pop_df <- pop_df%>% 
      set_names("Locus1","Locus2",pop_names)
    
    locstart<-popend+6
    
    locend<-locstart+ncomps-1
    
    loc_df <- tibble(dat = dis[locstart:locend]) %>%
      mutate(dat = substr(dat, start = 46, stop = 54)) %>% 
      mutate(dat = gsub(pattern = "Highly si", x = dat, replacement = "0.000000") %>% as.numeric()) 
    
    summary_df = pop_df %>% 
      mutate(Overall = loc_df %>% pull(dat))
    
  } #end npops
    
     return(summary_df)

}