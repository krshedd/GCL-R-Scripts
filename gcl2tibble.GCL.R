gcl2tibble.GCL <- function(sillyvec, loci, geno_sep = "/", n_ind_tbl = TRUE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function creates a tall tibble of genotypes, 1 row = 1 genotype per ind/locus
  #
  #  Argument(s):  
  #  sillyvec - character vector of ".gcl" object names without the .gcl extension
  #  loci - a character vector of locus names
  #  geno_sep - a character string to separate "Dose1" and "Dose2" genotypes (i.e. "G/T")
  #  n_ind_tbl - print table for number of individual per silly
  #
  #  Output:
  #  Prints a tibble showing the number of individuals per silly
  #  Returns a tall tibble of genotypes, 1 row = 1 genotype per individual and locus (now (june 2020) with allele1 and allele2 as dose1 and dose2!)
  #    ind - SillySource
  #    silly - sillyvec
  #    locus - a character vector of locus names
  #    genotype - a character vector of genotypes separated by `geno_sep`
  #               NOTE: NA = no call; haploid data (i.e. mitochondrial) has no geno_sep
  #    allele1 - reading from dose1
  #    allele2 - dose2
  #
  #  Written by Kyle Shedd 10/29/18  
  #  modified by bobby june 2020  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # veryify all loci in LocusControl
  if(!all(loci %in% LocusControl$locusnames)) {
    stop("The following loci are not found in `LocusControl`, hozah!!!\n", paste(setdiff(loci, LocusControl$locusnames), collapse ="\n"))
  }
  
  while(!require(tidyverse)){ install.packages("tidyverse") }

  # loop over sillyvec
  x <- lapply(sillyvec, function(silly) {
    
    my.gcl = get(paste0(silly, ".gcl"))
    
    if(my.gcl$n > 1) {
      dose1 <- tibble::as_tibble(my.gcl$scores[, loci, "Dose1"])
      dose2 <- tibble::as_tibble(my.gcl$scores[, loci, "Dose2"])
    } else { # when n = 1
      dose1 <- my.gcl$scores[, loci, "Dose1"] %>% dplyr::bind_rows()
      dose2 <- my.gcl$scores[, loci, "Dose2"] %>% dplyr::bind_rows()
    }
    
    # dose1 genotypes as tall tibble
    dose1 <- dose1 %>%
      dplyr::mutate(ind = my.gcl$attributes$SillySource) %>% 
      dplyr::mutate(silly = silly) %>% 
      tidyr::pivot_longer(cols= -c(ind, silly), 
                          names_to= "locus",
                          values_to= "allele1")

    # dose2 genotyeps as tall tibble
    dose2 <- dose2 %>% 
      dplyr::mutate(ind = my.gcl$attributes$SillySource) %>% 
      dplyr::mutate(silly = silly) %>% 
      tidyr::pivot_longer(cols= -c(ind, silly), 
                          names_to= "locus",
                          values_to= "allele2")

    # join dose1 and dose2 genotypes with sep = geno_sep
    geno <- dplyr::left_join(dose1, dose2,
                             by = c("ind", "silly", "locus")) %>%
      tidyr::unite(genotype, c(allele1, allele2),
                   sep = geno_sep, remove = FALSE)
  })
  
  # bind sillys into single tall tibble
  geno_tibble <- dplyr::bind_rows(x)
  
  # make no call genotypes NA
  geno_tibble <- geno_tibble %>% 
    dplyr::mutate(genotype = dplyr::case_when(
      genotype == paste("0", "0", sep = geno_sep) ~ "NA",
      # stringr::str_detect(genotype, paste0(geno_sep,"0")) &
      #   allele1 != 0 ~ allele1, # added 6/4/20
      # stringr::str_detect(genotype, paste0("0",geno_sep)) &
      #   allele2 != 0 ~ allele2, # added 6/4/20: commented out so they don't get confused with haplods.
      stringr::str_detect(genotype, paste0(geno_sep,"NA")) &
        allele1 != 0 ~ allele1, # added 6/4/20; this should take care the haploids and alleles with multiple SNPs
      stringr::str_detect(genotype, paste0(geno_sep,"NA")) &
        (allele1 == 0|is.na(allele1)) ~ "NA", # added 6/4/20
      # genotype == paste("0", "NA", sep = geno_sep) ~ "NA",
      # genotype == paste("A", "NA", sep = geno_sep) ~ "A",
      # genotype == paste("C", "NA", sep = geno_sep) ~ "C",
      # genotype == paste("G", "NA", sep = geno_sep) ~ "G",
      # genotype == paste("T", "NA", sep = geno_sep) ~ "T",
      # genotype == paste("NA", "NA", sep = geno_sep) ~ "NA",
      genotype == paste("Y", "Y", sep = geno_sep) ~ "Y", # added 6/4/20
      TRUE ~ as.character(genotype))) %>% 
    # dplyr::mutate(genotype = dplyr::na_if(genotype, "NA"))
    dplyr::mutate(genotype = dplyr::na_if(genotype, "NA"),
                  allele1 = dplyr::na_if(allele1, "0"),
                  allele2 = dplyr::na_if(allele2, "0")) # 6/4/20; this turned "0" to NA for allele 1 and 2 columns (easier for summarizing invariant/fixed loci). But genotype still shows 0 if it's like "A/0" or something like that.
  
  # tibble of number of individuals per silly
  ind_per_silly <- geno_tibble %>% 
    dplyr::group_by(silly) %>% 
    dplyr::summarise(n = n_distinct(ind))
  
  if (n_ind_tbl) print(ind_per_silly, n= Inf)
  
  return(geno_tibble)
}