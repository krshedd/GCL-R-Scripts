Plot_Allele_Freq.GCL <- function(freq, file, sillyvec = NULL, groupvec = NULL, loci = NULL, group_col = NULL){

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function uses the output from FreqPop.GCL() and produces allele freqency bubble plots.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   freq - a frequency tibble produced by FreqPop.GCL()
  #
  #   file - the file path including the file name with .pdf extention where the plots will be saved
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   groupvec - numeric vector indicating the group affiliation of each pop in sillyvec
  #
  #   loci - vector of locus names; if set to NULL all loci in the ".gcl" obejects will be used.
  #
  #   ncores - the number of cores to use in a foreach %dopar% loop. If the nubmer of core exceeds the number on your device, then ncores defaults to detectCores()
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  A pdf file containing allele fequency plots for each locus in loci
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  #  old2new_LocCtrl.GCL()  
  #  old2new_gcl.GCL(sillyvec67)
  #  Freq <- FreqPop.GCL(sillyvec = sillyvec67, loci = loci413)
  #
  #  Plot_Allele_Freq.GCL(freq, file, sillyvec = NULL, groupvec = NULL, loci = NULL, group_col = NULL)
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   When quantile is set to NULL this function utilizes rubias::close_matching_samples() to perform the duplicate check and it much faster than when you set a quantile.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(is.null(loci)){loci = freq$locus %>% unique()}
  
  if(is.null(groupvec)){groupvec = seq(length(freq$silly %>% unique()))}
  
  if(is.null(group_col)){group_col = rainbow(max(groupvec))}
  
  if(is.null(sillyvec)){sillyvec = freq$silly %>% unique()}
  
  if(!length(sillyvec)==length(groupvec)){
    
    stop("Make sure sillyvec and groupvec are the same length!!")
    
  }
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse)  # Install packages, if not in library and then load them
  
  freq$proportion[freq$proportion==0] <- NA
  
  freq_df <- freq %>% 
    left_join(tibble(silly = sillyvec, groupvec), by = "silly") %>% 
    mutate(groupvec = factor(groupvec))
  
  pdf(file = file)
  
  lapply(loci, function(locus){
    
    freq_df %>%
      filter(locus == !!locus) %>% 
      mutate(silly = factor(silly, levels = sillyvec), allele = factor(allele, unique(allele))) %>% 
      ggplot(aes(x = silly, y = allele, color = groupvec, size = proportion*100)) +
      geom_point() + 
      scale_color_manual(values = group_col, guide = FALSE) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
      ggtitle(label = locus)
    
    })
  
  dev.off()

}