BaselineEvalSampleSizes.GCL <- function(sillyvec, groupnames, groupvec, mixsize, scenarios = round(seq(.01, 1, .01), 2), maxprop = 0.5){
  
  #attach("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")

  #Final_Pops <- Final_Pops %>% mutate(group = factor(group, levels = unique(group)))
  
  # test_samps <- BaselineEvalSampleSizes.GCL(sillyvec = Final_Pops$silly, groupnames = Final_Pops$group %>% levels(), groupvec = Final_Pops$group %>% as.numeric(), scenarios = round(seq(.01, 1, .01), 2), mixsize = 200, maxprop = 0.5)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse) #Install packages, if not in library and then load them.
  
  maxp <- silly_n.GCL(sillyvec) %>% 
    dplyr::mutate(groupvec = !!groupvec) %>% 
    dplyr::group_by(groupvec) %>% 
    dplyr::summarize(groupn = sum(as.numeric(n)), .groups = "drop_last") %>%
    dplyr::mutate(maxn = round(groupn*maxprop, 0), mixsize = mixsize) %>% 
    dplyr::mutate(maxprop = ifelse(round(maxn/!!mixsize, 1)>1, 1, round(maxn/mixsize, 1))) %>% 
    dplyr::pull(maxprop) %>% 
    purrr::set_names(groupnames)
  
  ngroups <- length(groupnames) 
  
  lapply(groupnames, function(g){
    
    scenarios <- round(scenarios, digits = 2)
    
    gscn <- scenarios[scenarios <= maxp[g]]
    
    lapply(gscn, function(p){ 
      
      rs <- mixsize*p
      
      nr <- as.numeric(mixsize - rs)
      
      samps <- c(rs, rmultinom(n = 1, size = nr, prob = rep(1/(ngroups - 1), ngroups - 1))[ , 1]) %>% 
        purrr::set_names(c(g, setdiff(groupnames, g)))
      
      tibble::tibble(test_group = g, scenario = p, repunit = names(samps), samps = samps)
      
    }) %>% dplyr::bind_rows()
    
  }) %>% dplyr::bind_rows()
    
}
  
  
 