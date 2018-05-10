create_rubias_baseline <- function(sillyvec, loci, group_names, groupvec, path = "rubias/baseline", baseline_name) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function creates the baseline dataframe needed for `rubias`.
  # It reformats the "scores" from each individual in the baseline into a two column format used by `rubias`.
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
  # load("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Baseline 2012/Chignik2012SockeyeBaseline.RData")
  # chignik_7pops_22loci.rubias_base <- create_rubias_baseline(sillyvec = Chignik7Populations, loci = loci24, group_names = Groups, groupvec = Groupvec7, path = "rubias/baseline", baseline_name = "chignik_7pops_22loci")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  while(!require(tidyverse)){install.packages("tidyverse")}
  
  if(!dir.exists(path)) {stop("`path` to save baseline does not exist, hoser!!!")}
  
  silly_base.lst <- lapply(sillyvec, function(silly) {
    my.gcl <- get(paste0(silly, ".gcl"))
    scores.mat <- t(apply(my.gcl$scores[, loci, ], 1, function(ind) {c(t(ind))} ))
    colnames(scores.mat) <- as.vector(sapply(loci, function(locus) {c(locus, paste(locus, 1, sep = "."))} ))
    scores.df <- data.frame(scores.mat, stringsAsFactors = FALSE)
    scores.df$sample_type <- "reference"
    scores.df$repunit <- group_names[groupvec[sillyvec == silly]]
    scores.df$collection <- silly
    scores.df$indiv <- as.character(my.gcl$attributes$SillySource)
    silly_base.df <- scores.df[, c("sample_type", "repunit", "collection", "indiv", gsub(pattern = "-", replacement = ".", x = colnames(scores.mat)))] } #silly
  )
  baseline <- do.call("rbind", silly_base.lst)
  write_csv(x = baseline, path = paste0(path, "/", baseline_name, "_base.csv"))
  return(baseline)
}