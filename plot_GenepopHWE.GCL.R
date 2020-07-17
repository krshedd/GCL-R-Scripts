plot_GenepopHWE.GCL <- function(GenepopHWE_report, sillyvec, loci) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function plots results from ReadGenepopHWE.GCL(). You provide the report, a list of sillys, and loci and it will provide a visual of
  #   p-values faceted by silly.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   GenepopHWE_report - raw output from ReadGenepopHWE.GCL()
  #   sillyvec <- vector of sillys you're interested in
  #   loci <- vector of loci
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Produces a plot of the results from ReadGenepopHWE.GCL. Specifically, displays p-value, overall p-value, by silly.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #
  #
  #
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  # just setting variables for use in plotting - set ncol to something reasonable in facet_wrap 
  if( length(sillyvec) > 4){
    ncols <- round( length(sillyvec) / 3, 0)
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
    ggplot2::geom_hline(yintercept = length(loci) / 20, colour = "red") +
    ggplot2::geom_text(
      data = HWEpval %>% dplyr::filter( locus == "Overall Loci"),
      mapping = aes(x = 0.5, y = 15, label = pval),
      colour = "red",
      size = 6
    ) +
    ggplot2::facet_wrap(~ silly, ncol = ncols) +
    ggplot2::theme_bw()
}