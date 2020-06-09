create_rubias_baseline <- function(sillyvec, loci, group_names, groupvec, path = "rubias/baseline", baseline_name) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function creates the baseline dataframe needed for `rubias`.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   sillyvec - character vector of populations in the baseline
  #   loci - character vector of the loci you wish to include
  #   group_names - character vector of group names
  #   groupvec - numeric vector indicating the group affiliation of each pop in sillyvec
  #   path - character vector of where to save each mixture as a .csv
  #   baseline_name - character vector of what to name the baseline.csv
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns a dataframe in `rubias` baseline format
  #   Saves the baseline as a .csv file for posterity
  #     NOTE: to read in these .csv's use `read_csv(file = file, col_types = cols(.default = "c"))` 
  #     to make sure all columns are character vectors (if homozygous for T, it will become a logical vector).
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # CreateLocusControl.GCL(markersuite = "Coho_Baseline2016_95SNPs", username = "awbarclay", password = password)
  # dat <- readr::read_csv("V:/Analysis/2_Central/Coho/Cook Inlet/2019/2019_Cook_Inlet_coho_baseline/output/Final_Pops.csv")
  # load("V:/Analysis/2_Central/Coho/Cook Inlet/2019/2019_Cook_Inlet_coho_baseline/2019_Cook_Inlet_coho_baseline.Rdata")
  # sillyvec <- dat$collection
  # old2new_gcl.GCL(sillyvec)
  # group_names <- dat$groups_10 %>% unique()
  # groupvec <- dat$groups_10 %>% factor(levels = group_names) %>% as.numeric()
  # UCIcoho104pops_96loci.rubias_base <- create_rubias_baseline(sillyvec = sillyvec, loci = loci, group_names = group_names, groupvec = groupvec, path = "rubias/baseline", baseline_name = "UCIcoho104pops_96loci")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse)  # Install packages, if not in library and then load them.
  
  if(!dir.exists(path)) {stop("`path` to save baseline does not exist!!!")}
  
  scores_cols <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
    as.vector() 
  
  baseline <- lapply(sillyvec, function(silly){
    
    s <- match(silly, sillyvec)
    
    get(paste0(silly, ".gcl")) %>% 
      dplyr::mutate(sample_type = "reference", repunit = group_names[groupvec[s]], collection = silly,  indiv = SillySource) %>% 
      dplyr::select(sample_type, repunit, collection, indiv, tidyselect::all_of(scores_cols))
    
    }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::na_if(0)#I think this can be removed now that we convert all zeros to NAs when using LOKI2R.GCL
  
  readr::write_csv(x = baseline, path = paste0(path, "/", baseline_name, "_base.csv"))
  
  return(baseline)
  
}