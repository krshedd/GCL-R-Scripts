##Don't Run#
############
if(FALSE){##
############
############

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Title: S187 QC
  #Date: Wed Sep 05 10:42:19 2018
  #Name: Heather Hoyt; Kyle Shedd
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # It is important to run this script in order, or some functions will not provide accurate results.
  # This is by design, as this script only hits LOKI once to save time.
  
  # User input is only required above the '#~~~  GO! ~~~~~~~~~~~~~...`
  
  # Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # DupCheck Results (if applicable)
  
  # Figures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Levelplot of genotyping success rate by silly and locus for ProjectSillys
  # Levelplot of genotyping success rate by silly and locus for QCSillys
  # Histogram of overall QC individual conflict rate
  # Individual histograms of duplicate rate for conflict individuals
  
  # Summary excel ~~~~~~~~~~~~~~~~~~~~~~~
  # Summary by Silly
  # Conflicts by Silly
  # DupCheck Results (if applicable)
  # Conflicts by Locus
  # Conflicts by PlateID
  # Failure Rate by Silly
  # Failure Rate by Locus
  # Failure Rate by PlateID
  # Overall Failure Rate
  # Original Project Sample Size by Locus
  # Duplicates within Silly
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Setup ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  date()
  # rm(list=ls(all=TRUE))
  
  # This sources all of the new GCL functions to this workspace
  # source("C:/Users/krshedd/R/Functions.GCL.R")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Arguments ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  dirQC <- "V:/Lab/Genotyping/SNP Projects/Sockeye/Project S187 Copper River In Season 2018/QC/"

  species <- "sockeye"

  markersuite <- "Sockeye2011_96SNPs"

  project <- "S187"

  projectID <- 2421

  username <- "krshedd"
  
  .password <- ""
  
  QCSummaryfile <- paste("Project", project,"QC Summary Simple.xlsx") #  Do name normal summary file!!! If you do, it will overwrite it, not append it
  
  conflict_rate <- 0.10  # conflict rate at which dupcheck between sillys occurs

#~~~  GO! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  while(!require(pacman)){ install.packages("pacman") }

  p_load(tidyverse, lattice, xlsx, abind)  # use pacman to load or install + load necessary packages
  
  bbind <- function(...) { abind(..., along = 3) }

  source(path.expand("~/R/Functions.GCL.R"))

  setwd(dirQC)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Create Locus Control ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  CreateLocusControl.GCL(markersuite = markersuite, username = username, password = .password)
  
  loci <- LocusControl$locusnames

  nalleles <- LocusControl$nalleles

  ploidy <- LocusControl$ploidy

  alleles <- LocusControl$alleles

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Read in Project Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ReadProjectLOKI2R.GCL(projectID = projectID, username = username, password = .password)
  
  rm(.password)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Failure Rate ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  failure_rate <- FailureRate.GCL(sillyvec = ProjectSillys)
  
  failure_rate
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Read in QC Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  QCfiles <- list.files(path = "Genotype Data Files", pattern = ".csv", full.names = TRUE, recursive = FALSE)
  
  if(max(nalleles) <= 2) {
    
    ReadBiomarkQC.GCL(QCcsvFilepaths = QCfiles)
    
  } else {
    
    # Read in .csv files
    QC_genotypes <- suppressMessages(
      suppressWarnings(
        dplyr::bind_rows(
          lapply(QCfiles, function(fle) {readr::read_csv(file = fle)[, 1:4]} )
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
      assign(paste0(x, "QC.gcl"), list(counts = counts, scores = scores, n = nrow(scores), attributes = attributes))
      
    }  # x (silly)
    
    QCSillys <- paste0(ProjectSillysQC, "QC")
    
  }#else for usat
  
  QCColSize <- sapply(paste(QCSillys, ".gcl", sep = ''), function(x) get(x)$n)
  
  QCColSizeAll <- setNames(rep(0, length(ProjectSillys)),paste0(ProjectSillys, "QC.gcl"))
  
  QCColSizeAll[paste0(QCSillys, ".gcl")] <- QCColSize[paste0(QCSillys, ".gcl")]
  
  QCColSizeAll
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Read in Conflict Report ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  QCConcordanceReportfile <- list.files (path = "Conflict Reports", pattern = "ConcordanceReport", full.names = TRUE)
  
  CombineConflictsWithPlateID.GCL(files = QCConcordanceReportfile)
  
  # Old conflict report has "0" for mitochondrial conflicts, new has " " for mitochondrial conflicts, we will refer to them as "Homo-Homo".
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Conflict summaries
  
  conflicts_by_plate <- combined_conflicts %>% 
    dplyr::group_by(plate_id, concordance_type) %>%
    dplyr::summarise(n = dplyr::n()) %>% 
    tidyr::spread(concordance_type, n, fill = 0, drop = FALSE) %>% 
    dplyr::mutate(Conflict = sum(`Het-Het`, `Het-Homo`, `Homo-Het`, `Homo-Homo`)) %>% 
    dplyr::ungroup()
  
  conflicts_by_silly <- combined_conflicts %>% 
    dplyr::group_by(silly, concordance_type) %>% 
    dplyr::summarise(n = dplyr::n()) %>% 
    tidyr::spread(concordance_type, n, fill = 0, drop = FALSE) %>% 
    dplyr::mutate(Conflict = sum(`Het-Het`, `Het-Homo`, `Homo-Het`, `Homo-Homo`)) %>% 
    dplyr::ungroup()
  
  conflicts_by_locus <- combined_conflicts %>% 
    dplyr::group_by(locus, concordance_type) %>% 
    dplyr::summarise(n = dplyr::n()) %>% 
    tidyr::spread(concordance_type, n, fill = 0, drop = FALSE) %>% 
    dplyr::mutate(Conflict = sum(`Het-Het`, `Het-Homo`, `Homo-Het`, `Homo-Homo`)) %>% 
    dplyr::ungroup()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Sample Size by Locus for Project Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  OriginalProjectSampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = ProjectSillys, loci = loci)

  OriginalProjectPercentbyLocus <- apply(OriginalProjectSampleSizebyLocus, 1, function(row) {row / max(row)} )

  reruns <- which(apply(OriginalProjectPercentbyLocus, 2, min) < 0.8)  

  new_colors <- colorRampPalette(c("black", "white"))

  levelplot(t(OriginalProjectPercentbyLocus), col.regions = new_colors, at = seq(0, 1, length.out = 100), main = "% Genotyped", xlab = "SILLY", ylab = "Locus", scales = list(x = list(rot = 90)), aspect = "fill") # aspect = "iso" will make squares

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Sample Size by Locus for QC Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  OriginalQCSampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = QCSillys, loci = loci)

  OriginalQCPercentbyLocus <- apply(OriginalQCSampleSizebyLocus, 1, function(row) {row / max(row)} )

  rerunsQC <- which(apply(OriginalQCPercentbyLocus, 2, min) < 0.8)

  levelplot(t(OriginalQCPercentbyLocus), col.regions = new_colors, at = seq(0, 1, length.out = 100), main = "% Genotyped", xlab = "SILLY", ylab = "Locus", scales = list(x = list(rot = 90)), aspect = "fill") # aspect = "iso" will make squares

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### QA of Project Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ProjectSillys_SampleSizes <- matrix(data = NA, nrow = length(ProjectSillys), ncol = 5, dimnames = list(ProjectSillys, c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))
  
  ProjectSillys_SampleSizes[, "Genotyped"] <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n)
  
  if(species %in% c("chum", "sockeye")) {
    
    Alternate <- FindAlternateSpecies.GCL(sillyvec = ProjectSillys, species = species) %>% 
      dplyr::as_tibble()
    
    nAltBySilly <- sapply(ProjectSillys, function(silly) {
      AlternateSpeciesReport <- Alternate[grep(pattern = silly, x = rownames(Alternate)), ]
      sum(AlternateSpeciesReport$Alternate > 0.5 & AlternateSpeciesReport$Failure > 0.5)
    })
    # RemoveAlternateSpecies.GCL(AlternateSpeciesReport = Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5)  # Do not remove fish, just note how many per silly. Still want to catch them in conflicts later.
    
  } else {
    
    Alternate = tibble::tibble(x = "Not applicable")
    
  }
  
  ColSizePostAlternate <- ProjectSillys_SampleSizes[, "Genotyped"]
  if(exists(x = "nAltBySilly")) {ColSizePostAlternate <- ColSizePostAlternate - nAltBySilly}
  # ColSizePostAlternate <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n)

  ProjectSillys_SampleSizes[, "Alternate"] <- ProjectSillys_SampleSizes[, "Genotyped"] - ColSizePostAlternate 

  MissLoci <- RemoveIndMissLoci.GCL(sillyvec = ProjectSillys, proportion = 0.8)

  ColSizePostMissLoci <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n) - ProjectSillys_SampleSizes[, "Alternate"]

  ProjectSillys_SampleSizes[, "Missing"] <-  ColSizePostAlternate - ColSizePostMissLoci

  DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = ProjectSillys, loci = loci, quantile = NULL, minproportion = 0.95)

  DuplicateCheckReportSummary <- sapply(ProjectSillys, function(x) DuplicateCheck95MinProportion[[x]]$report, simplify = FALSE)
  
  nDupsBySilly <- sapply(DuplicateCheckReportSummary, function(silly) {ifelse(is.character(silly), 0, nrow(as.matrix(silly)))})
  # RemovedDups <- RemoveDups.GCL(DuplicateCheck95MinProportion)  # Do not remove fish, just note how many per silly. Still want to catch them in conflicts later.

  sapply(DuplicateCheckReportSummary[nDupsBySilly >=1], function(silly) {if(1 %in% abs(as.numeric(levels(silly$ID1)) - as.numeric(levels(silly$ID2)))) {"Sequential IDs found as duplicates, check 'DuplicateCheckReportSummary' for duplicated rows"} else {"Duplicates exist, but IDs do not appear sequential"} } )
  
  DuplicateCheckReportSummary[nDupsBySilly >= 1]  # Show within silly duplicates
  
  ColSizePostDuplicate <- ColSizePostMissLoci - nDupsBySilly
  # ColSizePostDuplicate <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n)

  ProjectSillys_SampleSizes[, "Duplicate"] <- ColSizePostMissLoci - ColSizePostDuplicate

  ProjectSillys_SampleSizes[, "Final"] <- ColSizePostDuplicate

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### QA of QC Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  MissLociQC <- RemoveIndMissLoci.GCL(sillyvec = QCSillys, proportion = 0.8)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Perform Duplicate Check on High Conflict Individuals ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Filter for conflicts, determine conflict rate
  conflicts <- combined_conflicts %>% 
    dplyr::filter(concordance_type %in% c("Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo")) %>% 
    dplyr::count(silly_source) %>% 
    dplyr::mutate(p = n / length(loci))
  
  # Histogram of conflict rate
  conflicts %>% 
    ggplot2::ggplot(aes(x = p)) +
    ggplot2::geom_bar() +
    ggplot2::xlim(0, 1) +
    ggplot2::geom_vline(xintercept = conflict_rate) +
    ggplot2::xlab("Conflict rate") +
    ggplot2::ylab("Frequency") +
    ggplot2::ggtitle("QC individual conflict rate")
  
  # Filter for conflicts > conflict_rate
  conflicts_investigate <- conflicts %>% 
    dplyr::filter(p > conflict_rate)
  
  # Duplicate check if necessary
  if(nrow(conflicts_investigate) == 0) {
    
    message(paste0("No individuals have > ", conflict_rate * 100, "% loci with conflicts between project and QC."))
    
  } else {
    
    message(paste0("The following individuals have > ", conflict_rate * 100, "% loci with conflicts between project and QC:\n"), paste(conflicts_investigate$silly_source, conflicts_investigate$n, "conflicts", collapse = "\n"))
    
    # Loop through individuals to see if missing loci
    conflict_indv <- NULL
    
    for (silly_ind in conflicts_investigate$silly_source) {
      
      silly <- stringr::str_split(string = silly_ind, pattern = "_", simplify = TRUE)[, 1]
      
      ind <- stringr::str_split(string = silly_ind, pattern = "_", simplify = TRUE)[, 2]
      
      # QC fish lost in QA?
      if(ind %in% MissLociQC[[paste0(silly, "QC")]]) {
        message(paste0("\n", silly, "QC_", ind, " does not have at least 80% loci genotyped, not running DupCheck for this individual."))
      }  # if QC fish removed due to missing genotypes
      
      # Project fish lost in QA
      if(ind %in% MissLoci[[silly]]) {
        message(paste0("\n", silly, "_", ind, " does not have at least 80% loci genotyped, not running DupCheck for this individual."))
      }  # if project fish removed due to missing genotypes
      
      conflict_indv <- c(conflict_indv, paste(silly, ind, sep = "_")[!(ind %in% MissLociQC[[paste0(silly, "QC")]] | ind %in% MissLoci[[silly]]) ])  # Confirm QC fish and Project fish were not removed
      
    }  # silly_ind
    
    # If no more, stop
    if(is.null(conflict_indv) | length(conflict_indv) == 0) {
      
      message("\nNo remaining high conflict individuals.")
      
      dup_check_results <- tibble::tibble(x = "Not applicable")
      
    } else {
      
      conflicts_investigate <- conflicts_investigate %>% 
        dplyr::filter(silly_source %in% conflict_indv)
      
      message("\nRunning DupCheckBetweenSillys.GCL on these high conflict individuals, as they have at least 80% loci genotyped for Project and QC extractions.")
      message(paste(conflicts_investigate$silly_source, conflicts_investigate$n, "conflicts", collapse = "\n"))
      
      conflict_silly <- unique(stringr::str_split(string = conflicts_investigate$silly_source, pattern = "_", simplify = TRUE)[, 1])
      
      KeySillyIDs <- setNames(
        lapply(conflict_silly, function(silly) {
          sapply(grep(pattern = silly, x = conflict_indv, value = TRUE), function(ind) {
            stringr::str_split(string = ind, pattern = "_", simplify = TRUE)[, 2]
          }, USE.NAMES = FALSE) 
        }),
        paste0(conflict_silly, "QC"))
      
      DupCheckResults <- sapply(conflict_silly, function(silly) {
        DupCheckBetweenSillys.GCL(KeySillys = paste0(silly, "QC"), 
                                  KeySillyIDs = KeySillyIDs[paste0(silly, "QC")], 
                                  BetweenSillys = ProjectSillys, 
                                  loci = loci, 
                                  threshold = 0.9)
      }, simplify = FALSE)  # FALSE
      
      dup_check_results <- dplyr::bind_rows(DupCheckResults, .id = "silly") %>% 
        tibble::as_tibble()
      
    }  # conflict_ind, post missing individuals
    
  }  # else, conflicts_to_investigate
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Create Summary Tables ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  summary_table_1 <- bind_cols(tibble(Silly = ProjectSillys), as.tibble(ProjectSillys_SampleSizes)) %>% 
    dplyr::left_join(FailureRate$silly_failure_rate, by = c("Silly" = "silly")) %>% 
    dplyr::rename("Failure Rate" = fail) %>% 
    dplyr::mutate("Total QC Fish" = QCColSizeAll)
  
  qc_silly_genotypes <- tibble(Silly.Code = factor(ProjectSillys),
                               qc_genotypes = sapply(ProjectSillys, function(silly) {
                                 qc_silly = paste0(silly, "QC.gcl")
                                 ifelse(qc_silly %in% names(QCColSizeAll), QCColSizeAll[qc_silly] * length(loci), 0)
                               } ))
  
  summary_table_2 <- conflicts_by_silly %>% 
    tidyr::gather(type, number, -Silly.Code) %>%  # make tall
    dplyr::left_join(qc_silly_genotypes) %>%  # join number of QC genotypes by silly
    dplyr::mutate(rate = number / qc_genotypes) %>%  # conflict numbers to rates
    tidyr::gather(variable, value, -Silly.Code, -qc_genotypes, -type) %>%  # make tall
    tidyr::unite(temp, type, variable) %>%  # unite conflict type with both number and rate
    tidyr::spread(temp, value) %>%  # make wide
    dplyr::rename(Silly = Silly.Code, 
                  "Total QC Genotypes" = qc_genotypes, 
                  "Total Discrepancies" = Conflict_number, 
                  "Discrepancy Rate" = Conflict_rate,
                  "DB Zeros" = `DB Zero_number`,
                  "DB Zero Rate" = `DB Zero_rate`,
                  "QC Zeros" = `File Zero_number`,
                  "QC Zero Rate" = `File Zero_rate`,
                  "Total Het-Het" = `Het-Het_number`,
                  "Het-Het Rate" = `Het-Het_rate`,
                  "Total Het-Homo" = `Het-Homo_number`,
                  "Het-Homo Rate" = `Het-Homo_rate`,
                  "Total Homo-Het" = `Homo-Het_number`,
                  "Homo-Het Rate" = `Homo-Het_rate`,
                  "Total Homo-Homo" = `Homo-Homo_number`,
                  "Homo-Homo Rate" = `Homo-Homo_rate`)
    
    summary_table_3 <- conflicts_by_locus %>% 
    tidyr::gather(type, number, -Locus) %>%  # make tall
    dplyr::mutate(qc_genotypes = sum(QCColSizeAll)) %>%  # join number of QC genotypes by locus
    dplyr::mutate(rate = number / qc_genotypes) %>%  # conflict numbers to rates
    tidyr::gather(variable, value, -Locus, -qc_genotypes, -type) %>%  # make tall
    tidyr::unite(temp, type, variable) %>%  # unite conflict type with both number and rate
    tidyr::spread(temp, value) %>%  # make wide
    dplyr::rename("Total QC Genotypes" = qc_genotypes, 
                  "Total Discrepancies" = Conflict_number, 
                  "Discrepancy Rate" = Conflict_rate,
                  "DB Zeros" = `DB Zero_number`,
                  "DB Zero Rate" = `DB Zero_rate`,
                  "QC Zeros" = `File Zero_number`,
                  "QC Zero Rate" = `File Zero_rate`,
                  "Total Het-Het" = `Het-Het_number`,
                  "Het-Het Rate" = `Het-Het_rate`,
                  "Total Het-Homo" = `Het-Homo_number`,
                  "Het-Homo Rate" = `Het-Homo_rate`,
                  "Total Homo-Het" = `Homo-Het_number`,
                  "Homo-Het Rate" = `Homo-Het_rate`,
                  "Total Homo-Homo" = `Homo-Homo_number`,
                  "Homo-Homo Rate" = `Homo-Homo_rate`)
  

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Append Summary Tables to QCSummaryfile.xlsx ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  write.xlsx(x = as.data.frame(summary_table_1), file = QCSummaryfile, sheetName = "Summary by Silly", row.names = FALSE, col.names = TRUE, append = TRUE)

  write.xlsx(x = as.data.frame(summary_table_2), file = QCSummaryfile, sheetName = "Conflicts by Silly", row.names = FALSE, col.names = TRUE, append = TRUE)

  if(exists("DupCheckResults")) {

    print(DupCheckResults)
    
    invisible(lapply(conflict_silly, function(silly) {write.xlsx(x = DupCheckResults[[silly]], file = QCSummaryfile, sheetName = paste("DupCheckBetween", silly, sep = " "), row.names = TRUE, col.names = TRUE, append = TRUE)} ))

  }

  write.xlsx(x = as.data.frame(summary_table_3), file = QCSummaryfile, sheetName = "Conflicts by Locus", row.names = FALSE, col.names = TRUE, append = TRUE)

  write.xlsx(x = as.data.frame(conflicts_by_plate), file = QCSummaryfile, sheetName = "Conflicts by PlateID", row.names = FALSE, col.names = TRUE, append = TRUE)

  write.xlsx(x = as.data.frame(FailureRate$silly_failure_rate), file = QCSummaryfile, sheetName = "Failure Rate by Silly", row.names = FALSE, col.names = TRUE, append = TRUE)

  write.xlsx(x = as.data.frame(FailureRate$locus_failure_rate), file = QCSummaryfile, sheetName = "Failure Rate by Locus", row.names = FALSE, col.names = TRUE, append = TRUE)

  write.xlsx(x = as.data.frame(FailureRate$plate_failure_rate), file = QCSummaryfile, sheetName = "Failure Rate by Plate", row.names = FALSE, col.names = TRUE, append = TRUE)

  write.xlsx(x = as.data.frame(FailureRate$overall_failure_rate), file = QCSummaryfile, sheetName = "Overall Failure Rate", row.names = FALSE, col.names = TRUE, append = TRUE)

  write.xlsx(x = OriginalProjectSampleSizebyLocus, file = QCSummaryfile, sheetName = "OriginalProjectSampleSizebyLocus", row.names = TRUE, col.names = TRUE, append = TRUE)

  if(exists("Alternate")){

    write.xlsx(x = Alternate, file = QCSummaryfile, sheetName = "Alternate Species", row.names = TRUE, col.names = TRUE, append = TRUE)

  }
  
  write.xlsx(x = DuplicateCheckReportSummary, file = QCSummaryfile, sheetName = "Duplicates Within Silly", row.names = TRUE, col.names = TRUE, append = TRUE)
  

#~~~  STOP!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  

##Don't Run#
############
}###########
############
############

