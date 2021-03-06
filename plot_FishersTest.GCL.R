plot_FishersTest.GCL <- function(pooling_test) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function plots results from FishersTest.GCL(). You can feed alist of sillys and it will facet by sillys.
  #   Null hypothesis: there are no significant differences in allele freqencies.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   pooling_test - raw output from FishersTest.GCL()
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Produces a plot of the results from FishersTest.GCL. Specifically, displays p-value, overall p-value, by silly.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  #   old2new_LocCtrl.GCL()
  #   old2new_gcl.GCL(sillyvec = c("KKILL05","KKILL06", "KFUNN05", "KFUNN06"), save_old = FALSE)
  #
  #   freq <- FreqPop.GCL(sillyvec = c("KKILL05","KKILL06", "KFUNN05", "KFUNN06"), loci = loci443)
  # 
  #   temp_pool <- FishersTest.GCL(freq = freq, loci = loci443, tests = list(c("KKILL05","KKILL06"), c("KFUNN05", "KFUNN06")))
  #   
  #   plot_FishersTest.GCL(pooling_test = temp_pool, loci = loci443)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  # just setting variables for use in plotting - set ncol to something reasonable in facet_wrap 
  if( length(pooling_test$test_sillys) > 4){
    ncols <- round( length(pooling_test$test_sillys) / 3, 0)
  } else {
    ncols <- NULL
  }
  
  pooling_test %>% tidyr::unnest(bylocus) %>%
    dplyr::filter(locus != "Overall") %>% # I don't think we actually need this?
    ggplot2::ggplot(aes(x = pval)) +
    ggplot2::geom_histogram(binwidth = 0.05) +
    ggplot2::geom_hline(yintercept = (pooling_test %>% tidyr::unnest(bylocus) %>% dplyr::select(locus) %>% dplyr::n_distinct()) / 20, colour = "red") + #XXXX okay -how do I reference the 'locus' column in the tibble I pipe to ggplot?? seems silly to remake it
    ggplot2::geom_text(
      mapping = aes(x = 0.5, y = 15, label = overall),
      colour = "red",
      size = 6
    ) +
    ggplot2::facet_wrap(~ test_sillys, ncol = ncols) +
    ggplot2::theme_bw()
}