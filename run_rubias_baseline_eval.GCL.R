run_rubias_baseline_eval.GCL <- function(tests, 
                                         group_names, 
                                         gen_start_col = 5, 
                                         base.path = "rubias/baseline", 
                                         mix.path = "rubias/mixture",
                                         out.path = "rubias/output", 
                                         method = "MCMC", 
                                         alle_freq_prior = list(const_scaled = 1), 
                                         pi_prior = NA, 
                                         pi_init = NULL, 
                                         reps = 25000, 
                                         burn_in = 5000, 
                                         pb_iter = 100, 
                                         prelim_reps = NULL, 
                                         prelim_burn_in = NULL, 
                                         sample_int_Pi = 10, 
                                         sample_theta = TRUE, 
                                         pi_prior_sum = 1, 
                                         seed = 56, 
                                         ncores = 4){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function is a wrapper for run_rubias_mixture(). 
  # This function runs baseline evaluation mixtures produced by CreateRubiasBaselineEval.GCL()
  # and multicores by test_group. (i.e., mixtures for each test group are run on a different core simultaineously).
  # 
  # This function is intended for use on a server that has many cores and a lot of RAM. 
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  #   All inputs required for `run_rubias_mixture` are passed on except for the path argument which is called out.path in this function.
  #
  #   tests - a tibble containing two variables: test_group (Chracter) and scenario (double) 
  #           This tibble can be obtained by summarizing the output of BaselineEvalSampleSizes.GCL(), see example below.
  #          
  #   group_names - character vector of group_names, used to sort repunit as a factor
  #
  #
  #   base.path - the file path where the baseline .csv files will be written
  #
  #   mix.path - the file path where the mixture .csv files will be written
  #
  #   out.path - character vector of where to save output from each mixture as a .csv
  #
  #   seed - integer to set the seed for the MCMC
  #
  #   ncores - a numeric vector of length one indicating the number of cores to use
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Breaks the output into each `mixture_collection` and for each saves as .csv files:
  #     1) collection level trace, wide format, collections in order of baseline (akin to .BOT file from BAYES)
  #     2) repunit level trace, wide format, repunit in order of `group_names` (akin to .RGN file from BAYES)
  #     3) straight dump of the `indiv_posteriors` tibble (without column `missing_loci`)
  #     4) straight dump of the `bootstrapped_proportions`
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   attach("C:/Users/awbarclay/Documents/Analysis/Chinook/Susitna_Chinook_baseline_2020/run baseline eval tests/Susitna_Chinook_baseline_2020_noGTseq.Rdata")
  #   Final_Pops <- Final_Pops %>% mutate(group = factor(group, levels = unique(group)))
  #   sample_sizes <- BaselineEvalSampleSizes.GCL(sillyvec = Final_Pops$silly, group_names = Final_Pops$group %>% levels(), groupvec = Final_Pops$group %>% as.numeric(), scenarios = round(seq(.01, 1, .01), 2), mixsize = 200, maxprop = 0.5, seed = 123)
  #   CreateRubiasBaselineEval.GCL(sillyvec = Final_Pops$silly, group_names = Final_Pops$group %>% levels(), test_groups = Final_Pops$group %>% levels(), loci = loci80, groupvec = Final_Pops$group %>% as.numeric(), sample_sizes = sample_sizes, prprtnl = TRUE, seed = 123, ncores = 8)
  # 
  #   tests <- sample_sizes %>% group_by(test_group, scenario) %>% summarize(test_group = test_group %>% unique(), scenario = scenario %>% unique(), .groups = "drop_last")#Total of 510 tests
  #   run_rubias_baseline_eval.GCL(tests = tests, group_names = Final_Pops$group %>% levels(), gen_start_col = 5, base.path = "rubias/baseline",mix.path = "rubias/mixture", out.path = "rubias/output", seed = 56, ncores = 8)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  start_time <- Sys.time()
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(rubias, tidyverse, parallel, doParallel, foreach) #Install packages, if not in library and then load them.
  
  test_groups = tests$test_group %>% 
    unique()
  
  # Multicore loop
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)
  
  foreach::foreach(g = test_groups, .export = "run_rubias_mixture", .packages = c("tidyverse", "rubias", "readr")) %dopar% {
    
    scenario_names = tests %>% 
      filter(test_group==g) %>% 
      unite(col = "name", test_group, scenario) %>% 
      pull(name)
    
    for(scn in scenario_names){
      
      mixture <- readr::read_csv(file = paste0(mix.path, "/", scn, ".mix.csv"), col_types = cols(.default = "c")) # Had to add col_types to make sure all cols are read as character, columns with T's were getting read as logical.
      
      baseline <- readr::read_csv(file = paste0(base.path, "/", scn, ".base.csv"), col_types = cols(.default = "c"))
      
      run_rubias_mixture(reference = baseline, 
                         mixture = mixture, 
                         group_names = group_names, 
                         gen_start_col = gen_start_col ,
                         method = method, 
                         alle_freq_prior = alle_freq_prior, 
                         pi_prior = pi_prior, 
                         pi_init = pi_init, 
                         reps = reps, 
                         burn_in = burn_in, 
                         pb_iter = pb_iter, 
                         prelim_reps = prelim_reps, 
                         prelim_burn_in = prelim_burn_in, 
                         sample_int_Pi = 10, 
                         sample_theta = sample_int_Pi, 
                         pi_prior_sum = pi_prior_sum, 
                         file = out.path, 
                         seed = seed)
      
    }
    
  }# end multicore loop
  
  parallel::stopCluster(cl)
  
  message(paste0("Analysis complete!! The test mixture output filese are located in this folder: ", out.path))
  
  Sys.time()-start_time
  
}
  