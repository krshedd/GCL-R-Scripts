plot_SampSizeByLocus <- function(SampSizeByLocus){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function takes the output from SampSizeByLocus() and creates a plotly heatmap of the proportion of fish with scores for each locus and silly
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   SampSizeByLocus - a tibble of sample sizes by locus produced by SampSizeByLocus(). 
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  #  A plotly heatmap of the proportion of fish with scores for each locus and silly
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  CreateLocusControl.GCL(markersuite = "Sockeye2011_96SNPs", username = "awbarclay", password = password)
  #  sillyvec = c("SMCDO03", "SNEVA13")
  #  password = "************"
  #  LOKI2R.GCL(sillyvec = sillyvec, username = "awbarclay", password = password)
  #  RemoveIndMissLoci.GCL(sillyvec = sillyvec)
  #
  #  Output <- SampSizeByLocus.GCL(sillyvec = sillyvec, loci = LocusControl$locusnames)
  #  plot_SampSizeByLocus(SampSizeByLocus = Output)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse, plotly)  # Install packages, if not in library and then load them.
  
  sillyvec <- SampSizeByLocus$silly
  
  silly_n <- silly_n.GCL(sillyvec)
  
  plotly::ggplotly(
    
    SampSizeByLocus %>% 
      gather(-silly, key = "Locus", value = "count") %>% 
      full_join(silly_n, by = "silly") %>% 
      dplyr::mutate(proportion = count/n) %>% 
      ggplot2::ggplot(aes(x = silly, y = Locus, fill = proportion))+
      ggplot2::geom_tile() +
      ggplot2::theme(axis.text.x = element_text(angle = 90))
    
  )
  
}
  
  
  
