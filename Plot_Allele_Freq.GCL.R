Plot_Allele_Freq.GCL <- function(freq, file, sillyvec = NULL, groupvec = NULL, loci = NULL, group_col = NULL, popnames = sillyvec, xlab.cex = 6, ylab.cex = 12){

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
  #   group_col - vector of colors corresponding to groupvec
  #
  #   loci - vector of locus names; if set to NULL all loci in the ".gcl" obejects will be used.
  #
  #   popnames - a vector of new population names (e.g. popnames <- c("KQUART","KCUP","KPINT"); Default is just sillyvec ; NULL assigns numbers seq_along(sillyvec), 
  #
  #   xlab.cex  - text size on x axis
  # 
  #   ylab.cex - text size y axis
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
  #  Plot_Allele_Freq.GCL(freq = Freq, file = "./test_freq_plots.pdf", sillyvec = sillyvec67, groupvec = groupvec, loci = loci413, group_col = grcol, popnames = NULL)
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   When quantile is set to NULL this function utilizes rubias::close_matching_samples() to perform the duplicate check and it much faster than when you set a quantile.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(is.null(loci)) {
    loci = freq$locus %>% unique()
  }
  
  if (is.null(groupvec)) {
    groupvec = seq(length(freq$silly %>% unique()))
  }
  
  if (is.null(group_col)) {
    group_col = rainbow(max(groupvec))
  }
  
  if(is.null(sillyvec)) {
    sillyvec = freq$silly %>% unique()
  }
  
  if(!length(sillyvec)==length(groupvec)) {
    stop("Make sure sillyvec and groupvec are the same length!!")
    
  }
  
  # set popnames to numeric if NULL
  if (is.null(popnames)) {
    popnames = seq_along(sillyvec)
  }
  
  # Making it esier to read numbers on the x-axis of plots
  if (is.numeric(popnames)) {
    angle = 45
  } else {
    angle = 90
  }
  
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse)  # Install packages, if not in library and then load them
  
  freq$proportion[freq$proportion == 0] <- NA 
  
  freq_df <- freq %>% 
    dplyr::left_join(dplyr::tibble(silly = sillyvec, groupvec), 
                     by = "silly") %>% 
    dplyr::mutate(groupvec = factor(groupvec))
  
  pdf(file = file)
  
  lapply(loci, function(locus){
    
    plot <- freq_df %>%
      filter(locus == locus) %>% 
      dplyr::filter(locus == !!locus) %>% 
      dplyr::mutate(silly = factor(silly, levels = sillyvec), 
                    allele = factor(allele, unique(allele))) %>% 
      ggplot2::ggplot(aes(x = silly, y = allele, color = groupvec, size = proportion*100)) +
      ggplot2::geom_point() + 
      ggplot2::scale_color_manual(values = group_col, guide = FALSE) +
      ggplot2::scale_x_discrete(labels = popnames) +
      ggplot2::theme(legend.position = "none", 
                     axis.text.x = ggplot2::element_text(angle = angle, hjust = 1, vjust = 0.5, size = xlab.cex), 
                     axis.text.y = ggplot2::element_text(size = ylab.cex)) +
      ggplot2::ggtitle(label = locus)
    
    print(plot)
    
    })
  
  dev.off()

}