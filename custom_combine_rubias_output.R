custom_combine_rubias_output <- function(rubias_output = NULL, mixvec = NULL, group_names = NULL, groupvec = NULL, 
                                         path = "rubias/output", alpha = 0.1, burn_in = 5000, 
                                         bias_corr = FALSE, threshold = 5e-7) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function computes summary statistics from `rubias` output, similar to `CustomCombineBAYESOutput.
  # However, output is a tibble with `mixture_collection` as a column, instead of each mixture as its own list.
  # It can take either the `rubias_output` list object from `run_rubias_mixture` or `infer_mixture`,
  # OR it can read in the .csv files created by `run_rubias_mixture`.
  #
  # NOTE: Currently this function only allows bias correction for the reporting groups run in the mixture
  # It can not do bias correction for different baseline groupvecs, because current `rubias` output only
  # gives the bias corrected means for each `mixture_collection` and `repunit` (i.e. `rho`, not `pi`)
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   rubias_output - output list object from `run_rubias_mixture` or `infer_mixture`
  #   mixvec - character vector of mixture sillys, used to read in output .csv files if `rubias_output = NULL`
  #   group_names - character vector of group_names, used to sort repunit as a factor
  #   groupvec - numeric vector indicating the group affiliation of each pop in sillyvec, used if resumarizing to new groups
  #   path - character vector of where to find output from each mixture as a .csv (created by `run_rubias_mixture`)
  #   alpha - numeric constant specifying credibility intervals, default is 0.1, which gives 90% CIs (i.e. 5% and 95%)
  #   burn_in - numeric constant specifying how many sweeps were used for burn_in in `run_rubias_mixture` or `infer_mixture`
  #   bias_corr - logical vector indicating whether you want bias corrected values from `method = "PB"` or not, 
  #               currently can NOT do bias correction if not using the same repunits that were run in the mixture.
  #   threshold - numeric constant specifying how low stock comp is before assume 0, used for `P=0` calculation, default is from BAYES
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns a tibble with 8 fields for each mixture and repunit (reporting group)
  #   mixture_collection - factor of mixtures (only a factor for ordering, plotting purposes)
  #   repunit - factor of reporting groups (only a factor for ordering, plotting purposes)
  #   mean - mean stock composition
  #   sd - standard deviation
  #   median - median stock composition
  #   loCI - lower bound of credibility interval
  #   hiCI - upper bound of credibility interval
  #   P=0 - the proportion of the stock comp distribution that was below `threshold` (i.e posterior probability that stock comp = 0)
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # load(".RData")
  # lynncanal_2015_sum <- custom_combine_rubias_output(rubias_output = lynncanal_2015_out)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  while(!require(tidyverse)){install.packages("tidyverse")}
  
  ## Error catching
  if(is.null(rubias_output) & is.null(mixvec)) {
    stop("Need to provide either `rubias_output` tibble to summarize or `mixvec` and `path` so that rubias output can be read, hoser!!!")
  }
  if(!is.null(groupvec) & bias_corr) {
    stop("Can not perform bias correction if you are changing the groupvec from what was originally run.\nUnfortunately since `rubias` only outputs bias corrected means for each `repunit`, we can't compute\nbias corrected summary statistics on a new `groupvec`.")
  }
  if(!is.null(groupvec) & is.null(group_names)) {
    stop("Need to provide `group_names` if introducing a new `groupvec`, hoser!!!")
  }
  
  
  ## If no rubias_output, make from .csv files
  if(is.null(rubias_output)) {
    
  }  # build rubias_output from .csv files
  
  
  ## Define other variables if NULL
  if(is.null(group_names)) {
    group_names <- unique(rubias_output$mix_prop_traces$repunit)
  }  # used to order groups as factor
  
  if(is.null(mixvec)) {
    mixvec <- unique(rubias_output$mix_prop_traces$mixture_collection)
  }  # used to order mixtures as factor
  
  ## Summarize output
  if(bias_corr) {
    
    if(nrow(rubias_output$bootstrapped_proportions) == 0) {stop("There is no bias corrected output, hoser!!!")}
    
    d_rho <- rubias_output$mixing_proportions %>% 
      dplyr::group_by(mixture_collection, repunit) %>% 
      dplyr::summarise(rho = sum(pi)) %>% 
      dplyr::ungroup() %>% 
      dplyr::left_join(rubias_output$bootstrapped_proportions, by = c("mixture_collection", "repunit")) %>% 
      dplyr::mutate(d_rho = rho - bs_corrected_repunit_ppn) %>% 
      dplyr::select(mixture_collection, repunit, d_rho)
    
    
  }
  
  
  loCI = alpha / 2
  hiCI = 1 - (alpha / 2)
  
  out_sum <- rubias_output$mix_prop_traces %>% 
    dplyr::filter(sweep >= burn_in) %>%  # remove burn_in
    dplyr::mutate(mixture_collection = factor(x = mixture_collection, levels = mixvec)) %>%  # order by mixvec
    dplyr::mutate(repunit = factor(x = repunit, levels = group_names)) %>%  # order by group_names
    dplyr::group_by(mixture_collection, sweep, repunit) %>%  # group to summarize across collections
    dplyr::summarise(rho = sum(pi)) %>%  # summarize collections to repunits
    dplyr::group_by(mixture_collection, repunit) %>%  # group by mixture and repunit across sweeps
    dplyr::summarise(mean = mean(rho),
                     sd = sd(rho),
                     median = median(rho),
                     loCI = quantile(rho, probs = loCI),
                     hiCI = quantile(rho, probs = hiCI),
                     `P=0` = sum(rho < threshold) / length(rho)) %>%  # summary statistics to return
    dplyr::ungroup() %>% 
    magrittr::set_colnames(c("mixture_collection", "repunit", "mean", "sd", "median", 
                             paste0(loCI * 100, "%"), paste0(hiCI * 100, "%"), "P=0"))
  
  return(out_sum)
}  # end function