CheckDupWithinSilly.GCL <- function(sillyvec, loci, quantile = 0.99, minproportion = 0.95, ncores = 4){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function checks for duplicate individuals within each silly in "sillyvec".
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   loci - vector of locus names; if set to NULL all loci in the ".gcl" obejects will be used.
  #
  #   quantile and minproportion - these arguments are used together to determine the cutt-off proportion at which a pair of duplicates 
  #                                is defined: i.e. proportion = max(quantile(duplication, quantile), minproportion. 
  #Setting "quantile" equal to NULL will skip the calculation of the duplication distribution and will run much faster--REQUIRES the package "reshape" ##### Need to use something other than reshape
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Returns a tibble of duplicate pairs of individuals by silly.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  # 
  #   dupcheck <- CheckDupWithinSilly.GCL(sillyvec = sillyvec157, loci = loci557, quantile = 0.99, minproportion = 0.95)
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(sum(is.na(match(loci, LocusControl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci, LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep=""))
    
  }
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse, reshape, doParallel, foreach) #Install packages, if not in library and then load them.

  nsilly <- length(sillyvec)

  nloci <- length(loci)

  scores_cols <- c(loci, paste0(loci, ".1")) %>% 
    sort()
  
  if (.Platform$OS.type == "windows") flush.console()
  
  my.gcl <- lapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)})
  names(my.gcl)=sillyvec
  
  if(is.null(quantile)){
    
    cl <- parallel::makePSOCKcluster(ncores)
    
    registerDoParallel(cl, cores=ncores)  
    
    # multicore loop
    dupcheck <-  foreach(silly = sillyvec) %dopar% {
      
      new.gcl <- my.gcl[[silly]]
        
      IDs <- new.gcl$FK_FISH_ID
      
      n <- silly_n.GCL("new")$n
      
      nloci <- length(loci)
      
      #Combine dose 1 and 2 into single columns for each locus separated by a period.
      scores0 <- new.gcl[ ,c("FK_FISH_ID", scores_cols)] 
    
      if(n < 2){
        
        report <- tibble::tibble(silly = silly, ID1 = NA, ID2 = NA, Missing1 = NA, Missing2 = NA, proportion = NA) #If only one individual the results will be a tibble of NAs. Can't do next() using applys
        
      } else{
        
        scores1 <- sapply(IDs, function(id){
          
          ind <- scores0 %>% 
            filter(FK_FISH_ID==id)
          
          sapply(loci,function(locus){
            
            ifelse(is.na(ind[ ,locus]), NA, paste(ind[ , locus], ind[ , paste0(locus,".1")], sep = "."))
            
            })
          
          }) %>% t() %>% 
          data.frame ()
        
        dimnames(scores1)[[1]] <- IDs
        
        sort.scores.df <- reshape::sort_df(scores1)
        
        sortIDs <- dimnames(sort.scores.df)[[1]]
        
     #Calculate the porpotion of duplicated loci for each silly
      duplication <- t(sapply(1:(n-1), function(id){ 
            
        vec <- sapply(1:nloci, function(locus){ 
          
          sort.scores.df[id, locus]==sort.scores.df[id+1,locus]
          
        });
        
        ifelse(sum(!is.na(vec)),sum(vec[!is.na(vec)])/sum(!is.na(vec)),0 )
        
      }))
      
      dupIND <- as.vector(duplication>minproportion)
      
      if(sum(dupIND)){
        
        dups <- t(sapply((1:(n-1))[dupIND], function(id){c(ID1 = sortIDs[id], ID2 = sortIDs[id+1])}))
        
        missing <- t(sapply(1:nrow(dups),function(dup){
          
          vec <- match(dups[dup,], sortIDs);
          
          c("Missing1" = sum(is.na(sort.scores.df[vec[1],])), "Missing2" = sum(is.na(sort.scores.df[vec[2], ])))
          
        }, simplify = TRUE))
        
        report <- dplyr::bind_cols(dplyr::as_tibble(dups), 
                                  dplyr::as_tibble(missing), 
                                  proportion = duplication[dupIND]
        ) %>% mutate(silly = silly) %>% 
          select(silly, tidyr::everything())
        
      }
      
      if(!sum(dupIND)){
        
       report <- tibble::tibble(silly = silly, ID1 = NA, ID2 = NA, Missing1 = NA, Missing2 = NA, proportion = NA) 
        
      }
      
      }
      
      report
      
    } #End multicore loop
    
    parallel::stopCluster(cl)
    
    output %>% 
      filter(!is.na(silly))
    
    return(output)
    
  }#End NULL quantile
  
  # Start quantile 

  if(!is.null(quantile)){
    for(silly in sillyvec){
      my.gcl=get(paste(silly,".gcl",sep=""),pos=1)
      scores=my.gcl$scores[,loci,]
      mydim=dim(scores)
      if(length(mydim)>2){
        ID=as.character(dimnames(scores)[[1]])
        n=length(ID)
      }
      if(length(mydim)==2){
        ID=as.character(my.gcl$attributes$FK_FISH_ID)
        n=length(ID)
      }
      if(length(mydim)<2){
        ID=NULL
        n=0
      }
     
      if(n<2){resultlist[[silly]]=list(report="No Duplicates", DupDist=NULL) ; next()}
 

      ncombs=choose(n,2)        

      mycombs=combn(ID,2)

      scores.df=data.frame(t(sapply(ID,function(id){sapply(loci,function(locus){ifelse(is.na(scores[id,locus,1]),NA,paste(scores[id,locus,1:ploidy[locus]],collapse="."))})})))

      
      duplication=sapply(1:ncombs,function(comb){compair=scores.df[mycombs[1,comb],]==scores.df[mycombs[2,comb],];sum(compair[!is.na(compair)])/sum(!is.na(compair))})
 
      names(duplication)=sapply(1:ncombs,function(comb){paste(mycombs[,comb],collapse=".")})
 
      proportion=max(quantile(duplication,quantile),minproportion)
 
      dupIND=duplication>proportion
    
      if(sum(dupIND)){

        dups=data.frame(ID1=mycombs[1,dupIND],ID2=mycombs[2,dupIND])
  
        report=data.frame(dups,t(sapply(seq(nrow(dups)),function(row){df=data.frame(Missing1=sum(is.na(scores.df[as.character(dups[row,1]),]))/nloci,Missing2=sum(is.na(scores.df[as.character(dups[row,2]),]))/nloci);vec=scores.df[as.character(dups[row,1]),]==scores.df[as.character(dups[row,2]),];data.frame(df,proportion=sum(vec[!is.na(vec)])/sum(!is.na(vec)))}) ))
 
      }
      if(!sum(dupIND)){
  
        report="No Duplicates"

      }

        resultlist[[silly]]=list(report=report,DupDist=duplication)  
    }

    return(resultlist)
  }#end if
}