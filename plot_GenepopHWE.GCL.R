plot_GenepopHWE.GCL <- function(GenepopHWE_report, sillyvec = NULL) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function plots results from ReadGenepopHWE.GCL(). You provide the report, a list of sillys, and loci and it will provide a visual of
  #   p-values faceted by silly.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   GenepopHWE_report - raw output from ReadGenepopHWE.GCL()
  #
  #   sillyvec <- vector of sillys you're interested in
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Produces a plot of the results from ReadGenepopHWE.GCL. Specifically, displays p-value, overall p-value, by silly.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  #  attach("V:\\Analysis\\5_Coastwide\\Chum\\NPen2WA_Chum_baseline\\NPen2WA_Chum_baseline.Rdata")
  #  sillyvec234 <- sillyvec234
  #  detach()
  #  
  #  HWEreport <- ReadGenepopHWE.GCL("V:/Analysis/5_Coastwide/Chum/NPen2WA_Chum_baseline/GENEPOP/NAKPen2WA_234pops_93loci.P", sillyvec = sillyvec234)
  #
  #  plot_GenepopHWE.GCL(GenepopHWE_report = HWEreport, sillyvec = sillyvec234[1:10])
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse)  # Install packages, if not in library and then load them
  
  summary_sillys <- dimnames(GenepopHWE_report$SummaryPValues)[[2]] 
  
  # Sillyvec is optional in ReadGenepopHWE.GCL. 
  # Check to make sure all sillys have results if sillyvec is supplied.
  if(!is.null(sillyvec)){
    
    n_sillys <- length(sillyvec)
    
    n_matches <-  sum(summary_sillys%in%sillyvec)
    
    if(n_sillys > n_matches){
      
      stop("GeneopHWE_report does not contain results for all sillys in sillyvec. Was sillyvec supplied for ReadGenepopHWE.GCL?")
      
    }
    
  } else{
  
  sillyvec <- summary_sillys
  
  }
  
  # just setting variables for use in plotting - set ncol to something reasonable in facet_wrap 
  if( length(sillyvec) > 4){
    
    ncols <- round(length(sillyvec) / 3, 0)
    
  } else {
    
    ncols <- NULL
    
  }

  # First convert the Genepop HWE report into a usable table and filter for sillyvec. 
  HWEpval <- GenepopHWE_report$SummaryPValues %>% 
    dplyr::as_tibble(rownames = "locus") %>%
    tidyr::pivot_longer(-locus, names_to = "silly", values_to = "pval")  %>% 
    dplyr::filter(silly %in% sillyvec)
  
  # Now plot the pvals
    ggplot2::ggplot(data = HWEpval %>% dplyr::filter(locus != "Overall Loci"), aes(x = pval)) +
    ggplot2::geom_histogram(binwidth = 0.05) +
    ggplot2::geom_hline(yintercept = (HWEpval %>% dplyr::select(locus) %>% dplyr::n_distinct()) / 20, colour = "red") +
    ggplot2::geom_text(
      data = HWEpval %>% dplyr::filter( locus == "Overall Loci"),
      mapping = aes(x = 0.5, y = 15, label = pval),
      colour = "red",
      size = 6
    ) +
    ggplot2::facet_wrap(~ silly, ncol = ncols) +
    ggplot2::theme_bw()
}