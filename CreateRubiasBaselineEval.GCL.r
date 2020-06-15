CreateRubiasBaselineEval.GCL <- function(sillyvec, groupnames,  loci, groupvec, sample_sizes, test_groups = groupnames, prprtnl = FALSE, base.path = "rubias/baseline", mix.path = "rubias/mixture", ncores = 4){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function creates rubias mixture and baseline files for different proof test scenarios.  
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  #
  # sillyvec=
  #
  # groupnames=
  #
  # groupvec<=
  # 
  # loci=
  #
  # samplesizemat=
  #
  # prprtnl<-FALSE
  #
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # N
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  attach("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  #  Final_Pops <- Final_Pops %>% mutate(group = factor(group, levels = unique(group)))
  #  sample_sizes <- BaselineEvalSampleSizes.GCL(sillyvec = Final_Pops$silly, groupnames = Final_Pops$group %>% levels(), groupvec = Final_Pops$group %>% as.numeric(), scenarios = round(seq(.01, 1, .01), 2), mixsize = 200, maxprop = 0.5)
  # 
  #  CreateRubiasBaselineEval.GCL(sillyvec = Final_Pops$silly, groupnames = Final_Pops$group %>% levels(), test_groups = Final_Pops$group %>% levels(), loci = loci80, groupvec = Final_Pops$group %>% as.numeric(), sample_sizes = sample_sizes, prprtnl = TRUE, ncores = 8)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse, parallel, doParallel, foreach) #Install packages, if not in library and then load them.
  
  if(sum(is.na(match(loci, LocusControl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci, LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  start_time <- Sys.time()
  
  group_check <- match(test_groups, sample_sizes$test_group %>% unique()) %>% 
    is.na() %>% 
    sum()
  
  if(group_check > 0){
    
    stop(paste("sample_sizes does not contain scenarios for all test_groups", sep = ""))
    
  }
  
  repunit_check <- match(groupnames, sample_sizes$repunit %>% unique()) %>% 
    is.na() %>% 
    sum()
  
  if(repunit_check > 0){
    
    stop(paste("sample_sizes does not contain sample sizes for all groupnames supplied", sep = ""))
    
  }
 
  #Write out baseline and mixture files for each test_group scenario in sample_sizes
  
  full_base <- create_rubias_baseline(sillyvec = sillyvec, loci = loci, group_names = groupnames, groupvec = groupvec, path = base.path, baseline_name = "full_base")
    
  # Multicore loop
  cl <- parallel::makePSOCKcluster(ncores)
    
  doParallel::registerDoParallel(cl, cores = ncores)  
    
  foreach::foreach(g = test_groups, .packages = c("tidyverse", "rubias")) %dopar% {
    
    tgscn <- sample_sizes %>% 
      dplyr::filter(test_group == g) %>% 
      dplyr::pull(scenario) %>% 
      unique()
    
    sapply(tgscn, function(scn){
      
      my.sample_sizes <- sample_sizes %>% 
        dplyr::filter(test_group == g, scenario == scn)
      
      mixture <- lapply(my.sample_sizes$repunit, function(repu){
        
        n <- my.sample_sizes %>% 
          dplyr::filter(repunit== repu) %>% 
          dplyr::pull(samps)
        
        IDs <- full_base %>% 
          dplyr::filter(repunit==repu) %>% 
          dplyr::select(indiv, collection)
        
        if(prprtnl){
        
          # This function rounds the number of samples to select so they add up to the original n.
          smart.round <- function(x) {
            
            y <- floor(x)
            
            indices <- tail(order(x-y), round(sum(x)) - sum(y))
            
            y[indices] <- y[indices] + 1
            
            y
            
          }
          
          props <- IDs %>% 
            dplyr::group_by(collection) %>% 
            dplyr::summarize(total_n = length(collection), .groups = "drop_last") %>% 
            dplyr::mutate(prop = total_n/sum(total_n), sample_n = smart.round(!!n*prop))
          
          mix_samples <-  lapply(props$collection, function(col){
            
            IDs %>% 
              dplyr::filter(collection==col) %>% 
              dplyr::pull(indiv) %>% 
              sample(size = props %>% filter(collection == col) %>% pull(sample_n))
            
          }) %>% unlist()
          
        } else{mix_samples <- sample(IDs, size = n)}
        
        full_base %>% 
          dplyr::filter(repunit == repu) %>% 
          dplyr::filter(indiv%in%mix_samples)
        
      }) %>% 
        dplyr::bind_rows() %>% 
        dplyr::mutate(sample_type = "mixture")
      
      baseline <- full_base %>% 
        dplyr::filter(!indiv%in%mixture$indiv)
      
      readr::write_csv(mixture, path = paste0(mix.path, "/", g, "_", scn, ".mix.csv"))
      
      readr::write_csv(baseline, path = paste0(base.path, "/", g, "_", scn, ".base.csv"))
      
    })
    
  } #End multicore
  
  parallel::stopCluster(cl)
  
  Sys.time()-start_time
  
}