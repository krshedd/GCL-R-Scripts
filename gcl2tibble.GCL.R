gcl2tibble.GCL <- function(sillyvec, loci, geno_sep = "/"){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function creates a tall tibble of genotypes, 1 row = 1 genotype per ind/locus
  #
  #  Argument(s):  
  #  sillyvec - character vector of ".gcl" object names without the .gcl extension
  #  loci - a character vector of locus names
  #  geno_sep - a character string to separate "Dose1" and "Dose2" genotypes (i.e. "G/T")
  #
  #  Output:
  #  Prints a tibble showing the number of individuals per silly
  #  Returns a tall tibble of genotypes, 1 row = 1 genotype per individual and locus
  #    ind - SillySource
  #    silly - sillyvec
  #    locus - a character vector of locus names
  #    genotype - a character vector of genotypes separated by `geno_sep`
  #               NOTE: NA = no call; haploid data (i.e. mitochondrial) has no geno_sep
  #
  #  Written by Kyle Shedd 10/29/18  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # veryify all loci in LocusControl
  if(!all(loci %in% LocusControl$locusnames)) {
    stop("The following loci are not found in `LocusControl`, hoser!!!\n", paste(setdiff(loci, LocusControl$locusnames), collapse ="\n"))
  }
  
  while(!require(tidyverse)){ install.packages("tidyverse") }

  # loop over sillyvec
  x <- lapply(sillyvec, function(silly) {
    
    my.gcl = get(paste0(silly, ".gcl"))
    
    # dose1 genotypes as tall tibble
    dose1 <- tibble::as_tibble(my.gcl$scores[, loci, "Dose1"]) %>% 
      dplyr::mutate(ind = my.gcl$attributes$SillySource) %>% 
      dplyr::mutate(silly = silly) %>% 
      tidyr::gather(locus, dose1_genotype, - ind, - silly)
    
    # dose2 genotyeps as tall tibble
    dose2 <- tibble::as_tibble(my.gcl$scores[, loci, "Dose2"]) %>% 
      dplyr::mutate(ind = my.gcl$attributes$SillySource) %>% 
      dplyr::mutate(silly = silly) %>% 
      tidyr::gather(locus, dose2_genotype, - ind, - silly)
    
    # join dose1 and dose2 genotypes with sep = geno_sep
    geno <- dplyr::left_join(dose1, dose2, by = c("ind", "locus", "silly")) %>% 
      tidyr::unite(genotype, c(dose1_genotype, dose2_genotype), sep = geno_sep)
    
  })
  
  # bind sillys into single tall tibble
  geno_tibble <- dplyr::bind_rows(x)
  
  # make no call genotypes NA
  geno_tibble <- geno_tibble %>% 
    dplyr::mutate(genotype = dplyr::case_when(
      genotype == paste("0", "0", sep = geno_sep) ~ "NA",
      genotype == paste("0", "NA", sep = geno_sep) ~ "NA",
      genotype == paste("A", "NA", sep = geno_sep) ~ "A",
      genotype == paste("C", "NA", sep = geno_sep) ~ "C",
      genotype == paste("G", "NA", sep = geno_sep) ~ "G",
      genotype == paste("T", "NA", sep = geno_sep) ~ "T",
      genotype == paste("NA", "NA", sep = geno_sep) ~ "NA",
      TRUE ~ as.character(genotype))) %>% 
    dplyr::mutate(genotype = dplyr::na_if(genotype, "NA"))
  
  # tibble of number of individuals per silly
  ind_per_silly <- geno_tibble %>% 
    dplyr::group_by(silly) %>% 
    dplyr::summarise(n = n_distinct(ind))
  
  print(ind_per_silly)
  
  return(geno_tibble)
}