write_genepop.GCL <- function(sillyvec, loci, path, VialNums = TRUE, usat = FALSE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function creates a GENEPOP file.
  # It reformats the "scores" from each individual into a one column numeric format used by GENEPOP.  
  # SNP alleles at each locus are numbered as they occur in LocusControl$alleles and uSAT alleles remain the same.
  # Missing scores are replaced with zeros.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   sillyvec - character vector of populations in the baseline
  #   loci - character vector of the loci you wish to include
  #   path - character vector of where to save the GENEPOP file
  #   VialNums - if set to TRUE vial numbers are included with the silly for each individual, otherwise just the silly is included
  #   usat - set to TRUE if loci are microsatellites
  #
  # Example - SNP~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # attach("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  # gcl2Genepop.GCL(sillyvec = sillyvec70, loci = loci439, path = "V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/GENEPOP/test.gen", VialNums = TRUE, usat = FALSE)
  # Example - usat~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # attach("V:/Analysis/1_SEAK/Chinook/Baseline/GAPS3.0/GAPSbaseline_Feb2013.RData")
  # gcl2Genepop.GCL(sillyvec = SEAKPops356, loci = GAPSloci, path = "V:/Analysis/1_SEAK/Chinook/Baseline/GAPS3.0/Genepop/test.gen", VialNums = TRUE, usat = TRUE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  while(!require(tidyverse)){install.packages("tidyverse")}
  
  #Duplicate sillys
  silly_dup <- sillyvec %>% duplicated()
  
  if(sum(silly_dup) > 0){stop(paste0("The following sillys occur more than once in sillyvec: ", sillyvec[silly_dup]))}
  
  #Similar sillys
  silly_trunc_dup <- str_trunc(sillyvec, width = 14, side = "right", ellipsis = "") %>% 
    duplicated()
  
  if(silly_trunc_dup %>% sum() > 0){
    
    warning(paste0("GENEPOP truncates pop names to 14 characters, pops with similar sillys and will have the same name in the GENEPOP output."))
      
   }
   
  #Locus check
  if(sum(is.na(match(loci, LocusControl$locusnames)))){
    
    stop(paste0("'", loci[is.na(match(loci, LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!"))
    
    }
  
  #Build file
  file <- "GENEPOP input format"
  
  file <- rbind(file, cbind(loci))
  
  #Initial scores
  scores_df0 <- lapply(sillyvec, function(silly){
    
    mysilly <- get(paste(silly, ".gcl", sep = ""), pos = 1);
    
    mysilly$scores %>% 
      as_tibble(rownames = "FishID") %>% 
      left_join(as_tibble(mysilly$attributes[,c("FK_FISH_ID", "SillySource")]) %>%
                  mutate(FK_FISH_ID = as.character(FK_FISH_ID)), by = c("FishID"="FK_FISH_ID"))
    
  }) %>% set_names(sillyvec) %>%
    bind_rows(.id = "Silly") 
  
  dose <- dimnames(get(paste(sillyvec[1], ".gcl", sep = ""), pos = 1)$scores)[[3]][1] %>% 
   str_split(pattern="1") %>% 
   unlist()#usats have a space for "Dose 1" and "Dose 2" whereas SNPs don't have a space.
  
  ploidy <- LocusControl$ploidy[loci]
  
  locus_order0 <- lapply(loci, function(locus){paste0(locus, ".",dose[1], 1:ploidy[locus])})

  locus_order <- locus_order0 %>% 
    unlist()
  
if(usat){
  
  #add fish ids
  scores_df1 <- scores_df0 %>% 
    mutate(SillySource = if(VialNums){paste(Silly, FishID, sep = "_")}else{silly}) %>% 
    mutate(Silly = factor(Silly, levels = sillyvec), SillySource = factor(SillySource, levels = unique(SillySource))) %>% 
    select(Silly, SillySource, locus_order) 
  
  #df to long format and pad alleles with leading zeros if they have less than three characters
  scores_df2 <- scores_df1 %>% 
    gather(-Silly, -SillySource, key = "locus.dose", value = "call") %>% 
    separate(locus.dose, sep = paste0(".",dose[1]), into = c("locus", "dose")) %>% 
    mutate(call = str_pad(call, width = 3, side = "left", pad = "0"))
  
  #Replace NA's with 000 an concatinate doses into single column
  scores_df3 <- scores_df2 %>% mutate(call = replace_na(call,replace = "000"), locus = factor(locus, levels = loci)) %>% 
    spread(key = dose, value = call) %>% 
    unite(`1`,`2`, col = "call", sep = "")
  
  #Spread loci into individual columns
  scores_df <- scores_df3 %>% 
    spread(key = locus, value = call) %>% 
    unite(loci, col = "Geno", sep = " ") %>% 
    unite(SillySource, Geno, col = "dat", sep = " ,  ")
  
}else{ 
  
  #add fish ids 
  scores_df1 <- scores_df0 %>% 
    mutate(SillySource = if(VialNums){paste(Silly, FishID, sep = "_")}else{silly}) %>% 
    mutate(Silly = factor(Silly, levels = sillyvec), SillySource = factor(SillySource, levels = unique(SillySource))) %>% 
    select(Silly, SillySource, locus_order) %>% 
    mutate_at(.vars = locus_order, .funs = ~sub(0, NA, .)) 
  
  #df to long format
  scores_df2 <- scores_df1 %>% 
    gather(-Silly, -SillySource, key = "locus.dose", value = "call") %>% 
    separate(locus.dose, sep = paste0(".",dose[1]), into = c("locus", "dose")) 
  
  #Creating tibble of alleles for each locus
  alleles <- lapply(loci, function(locus){nalleles = length(LocusControl$alleles[[locus]]); LocusControl$alleles[[locus]]  %>% as_tibble() %>% mutate(allele_num = 1:nalleles)}) %>% 
    set_names(loci) %>% 
    bind_rows(.id = "locus") %>% 
    mutate(allele = value) %>% 
    select(locus, allele, allele_num)
  
  #Make alleles numeric, two digit format  
  scores_df3 <- lapply(loci, function(locus0){
      call = scores_df2 %>% filter(locus==locus0) %>% pull(call);
      all = filter(alleles, locus==locus0) %>% pull(allele); 
      call_num = match(call, all);
      scores_df2 %>% filter(locus==locus0) %>% mutate(call = call_num %>% sprintf("%02d", .))
      }) %>% bind_rows()
  
  #Replace NA's with 00 an concatinate doses into single column 
  scores_df4 <- scores_df3 %>% mutate(call = sub(pattern = "NA", x = call, replacement = "00"), locus = factor(locus, levels = loci)) %>% 
    spread(key = dose, value = call) %>% 
    unite(`1`,`2`, col = "call", sep = "")
   
  #Spread loci into individual columns
  scores_df <- scores_df4 %>% 
    spread(key = locus, value = call) %>% 
    unite(loci, col = "Geno", sep = " ") %>% 
    unite(SillySource, Geno, col = "dat", sep = " ,  ")
  
}#SNP
  
  counts <- sapply(sillyvec, function(silly){
    
    c("Pop", scores_df %>% 
    filter(Silly==silly) %>% 
    pull(dat))
  
  }) %>% unlist() %>% 
    as.character()
  
  file <- rbind(file, cbind( counts))
  
  write.table(file, path, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  return(NULL)   
  
}



