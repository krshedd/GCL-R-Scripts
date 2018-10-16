ReadUSatQC.GCL <- function(QCcsvFilepaths) {
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  This function reads in uSat QC genotypes from .csv files as .gcl objects
#
#  Argument(s):  
#  QCcsvFilepaths <- character vector with relative path for QC genotype .csv files
#
#  Output:
#  Creates QC.gcl list objects in global environment
#  QCSillys - a character vector of QC sillys is assigned to the global environment
#
#  Written by Kyle Shedd 10/15/18  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  while(!require(tidyverse)){ install.packages("tidyverse") }
  
  # Read in .csv files
  QC_genotypes <- suppressMessages(
    suppressWarnings(
      dplyr::bind_rows(
        lapply(QCcsvFilepaths, function(fle) {readr::read_csv(file = fle)[, 1:4]} )
      )  # bind_rows
    )  # supressWarnings
  )  # suppressMessages
  
  # Rename columns, split silly_source
  QC_genotypes <- QC_genotypes %>% 
    dplyr::rename(silly_source = "Sample Name", locus = Marker, allele_1 = "Allele 1", allele_2 = "Allele 2") %>% 
    tidyr::separate(col = silly_source, into = c("silly", "fish_id"), sep = "_", remove = FALSE)
  
  # Verify that all QC silly are in the project
  ProjectSillysQC <- unique(QC_genotypes$silly)
  if(!all(ProjectSillysQC %in% ProjectSillys)){ stop(paste0(ProjectSillysQC[! ProjectSillysQC %in% ProjectSillys], " not found in ProjectSillys.")) }
  
  # Verify that all QC loci are in project loci
  lociQC <- sort(unique(QC_genotypes$locus))
  if(!all(lociQC %in% loci)){ stop(paste0(lociQC[! lociQC %in% loci], " not found in LocusControl.")) }
  
  # attributes table names
  attNames <- colnames(get(paste0(ProjectSillys[1], ".gcl"))$attributes)
  
  # Loop over silly to create .gcl objects
  for(x in ProjectSillysQC){
    
    # subset genotypes by silly
    x_genotypes <- QC_genotypes %>% 
      dplyr::filter(silly == x) %>% 
      dplyr::mutate(locus = factor(locus, levels = loci))
    
    # create scores array
    scores <- Reduce(bbind, lapply(c("allele_1", "allele_2"), function(al){ tapply(pull(x_genotypes, al), list(x_genotypes$fish_id, x_genotypes$locus), c)[,, drop = FALSE] }))
    dimnames(scores)[[3]] <- paste0("Dose", 1:2)
    
    # create counts array
    counts <- array(NA, c(nrow(scores), ncol(scores), max(nalleles)), list(rownames(scores), colnames(scores), paste0("Allele", seq(max(nalleles)))))
    for(locus in loci){
      for(al in seq(nalleles[locus])){
        for(id in x_genotypes$fish_id){
          counts[id, locus, al] <- sum(scores[id, locus, seq(ploidy[locus])] == alleles[[locus]][al])
        }  # id
      }  # al           
    }  # locus
    
    # create attributes data.frame
    attributes <- data.frame(matrix(NA, nrow = nrow(counts), ncol = length(attNames), dimnames = list(rownames(counts), attNames)))
    attributes$FK_FISH_ID <- rownames(counts)
    attributes$SillySource <- paste0(x, "QC_", rownames(counts))
    
    # assign to global environment
    assign(x = paste0(x, "QC.gcl"), value = list(counts = counts, scores = scores, n = nrow(scores), attributes = attributes), pos = 1)
    
  }  # x (silly)
  
  assign(x = "QCSillys", value = paste0(ProjectSillysQC, "QC"), pos = 1)
}