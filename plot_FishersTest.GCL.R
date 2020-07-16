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
  #   plot_FishersTest.GCL(pooling_test = temp_pool)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  # just setting variables for use in plotting - set ncol to something reasonable in facet_wrap 
  if( length(pooling_test$test_sillys) > 4){
    ncols <- round( length(pooling_test$test_sillys) / 3, 0)
  } else {
    ncols <- NULL
  }
  
  pooling_test %>% unnest(bylocus) %>%
    filter(locus != "Overall") %>%
    ggplot(aes(x = pval)) +
    geom_histogram(binwidth = 0.05) +
    geom_hline(yintercept = length(loci96) / 20, colour = "red") +
    geom_text(
      mapping = aes(x = 0.5, y = 15, label = overall),
      colour = "red",
      size = 6
    ) +
    facet_wrap(~ test_sillys, ncol = ncols) +
    theme_bw()
}