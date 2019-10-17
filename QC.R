##Don't Run#
############
if(FALSE){##
  ############
  ############
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Title: P014 QC
  # Date: Mon Oct 15 17:06:46 2018
  # Name: Heather Hoyt; Kyle Shedd; Emily Lescak
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
  
  dirQC <- "V:/Lab/Genotyping/SNP Projects/Pink/Project P014 AHRP Parentage GTSeq Part 2/QC/"
  
  species <- "pink"
  
  project <- "P014"

  project_type <- "Biomark"  # this can be either "Biomark" for regular SNP projects with Fluidigm, "uSat", or "GT-seq"

  username <- "krshedd"
  
  .password <- ""
  
  QCSummaryfile <- paste("Project", project,"QC Summary R Script.xlsx") #  Do name normal summary file!!! If you do, it will overwrite it, not append it
  
  conflict_rate <- 0.10  # conflict rate at which dupcheck between sillys occurs
  
  #~~~  GO! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  while(!require(pacman)){ install.packages("pacman") }
  
  p_load(tidyverse, lattice, writexl, abind)  # use pacman to load or install + load necessary packages
  
  bbind <- function(...) { abind(..., along = 3) }
  
  source(path.expand("~/../R/Functions.GCL.R"))  # user may need to change depending on where you put this directory
  
  setwd(dirQC)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Read in Project Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  read_project_genotypes.GCL(project_name = project, username = username, password = .password)
  
  rm(.password)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Failure Rate ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  (failure_rate <- FailureRate.GCL(sillyvec = ProjectSillys))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Read in QC Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  QCfiles <- list.files(path = "Genotype Data Files", pattern = ".csv", full.names = TRUE, recursive = FALSE)
  
  if (project_type == "Biomark") {
    ReadBiomarkQC.GCL(QCcsvFilepaths = QCfiles)
  } else {
    if (project_type == "uSat") {
      ReadUSatQC.GCL(QCcsvFilepaths = QCfiles)
    } else {
      if (project_type == "GT-seq") {
        ReadGTseqQC.GCL(QCcsvFilepaths = QCfiles)
      } else {
        message("project_type must be either 'Biomark', 'uSat', or' GT-seq' to read in QC genotypes")
      }  # GT-seq
    }  # uSat
  }  # Biomark

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
    dplyr::summarise(n = n()) %>% 
    tidyr::spread(concordance_type, n, fill = 0, drop = FALSE) %>% 
    dplyr::mutate(Conflict = sum(`Het-Het`, `Het-Homo`, `Homo-Het`, `Homo-Homo`)) %>% 
    dplyr::ungroup()
  
  conflicts_by_silly <- combined_conflicts %>% 
    dplyr::group_by(silly, concordance_type) %>% 
    dplyr::summarise(n = n()) %>% 
    tidyr::spread(concordance_type, n, fill = 0, drop = FALSE) %>% 
    dplyr::mutate(Conflict = sum(`Het-Het`, `Het-Homo`, `Homo-Het`, `Homo-Homo`)) %>% 
    dplyr::ungroup()
  
  conflicts_by_locus <- combined_conflicts %>% 
    dplyr::group_by(locus, concordance_type) %>% 
    dplyr::summarise(n = n()) %>% 
    tidyr::spread(concordance_type, n, fill = 0, drop = FALSE) %>% 
    dplyr::mutate(Conflict = sum(`Het-Het`, `Het-Homo`, `Homo-Het`, `Homo-Homo`)) %>% 
    dplyr::ungroup()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Sample Size by Locus for Project Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  OriginalProjectSampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = ProjectSillys, loci = loci)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Sample Size by Locus for QC Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  OriginalQCSampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = QCSillys, loci = loci)
  
  OriginalQCPercentbyLocus <- apply(OriginalQCSampleSizebyLocus, 1, function(row) {row / max(row)} )
  
  rerunsQC <- which(apply(OriginalQCPercentbyLocus, 2, min) < 0.8)
  
  new_colors <- colorRampPalette(c("black", "white"))
  
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
  
  ProjectSillys_SampleSizes
  
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
    ggplot2::geom_vline(xintercept = conflict_rate, colour = "red", lwd = 1.5) +
    ggplot2::xlab("Conflict rate") +
    ggplot2::ylab("Frequency") +
    ggplot2::ggtitle("QC individual conflict rate")
  
  # Filter for conflicts > conflict_rate
  conflicts_investigate <- conflicts %>% 
    dplyr::filter(p > conflict_rate)
  
  # Duplicate check if necessary
  if(nrow(conflicts_investigate) == 0) {
    
    message(paste0("No individuals have > ", conflict_rate * 100, "% loci with conflicts between project and QC."))
    
    dup_check_results <- tibble::tibble(x = "Not applicable")
    
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
  
  summary_table_1 <- dplyr::bind_cols(tibble(Silly = ProjectSillys), as_tibble(ProjectSillys_SampleSizes)) %>% 
    dplyr::left_join(failure_rate$silly_failure_rate, by = c("Silly" = "silly")) %>% 
    dplyr::rename("Failure Rate" = fail) %>% 
    dplyr::mutate("Total QC Fish" = QCColSizeAll)
  
  qc_silly_genotypes <- tibble::tibble(silly = factor(ProjectSillys),
                               qc_genotypes = sapply(ProjectSillys, function(silly) {
                                 qc_silly = paste0(silly, "QC.gcl")
                                 ifelse(qc_silly %in% names(QCColSizeAll), QCColSizeAll[qc_silly] * length(loci), 0)
                               } ))
  
  summary_table_2 <- conflicts_by_silly %>% 
    tidyr::gather(type, number, -silly) %>%  # make tall
    dplyr::left_join(qc_silly_genotypes, by = "silly") %>%  # join number of QC genotypes by silly
    dplyr::mutate(rate = number / qc_genotypes) %>%  # conflict numbers to rates
    tidyr::gather(variable, value, -silly, -qc_genotypes, -type) %>%  # make tall
    tidyr::unite(temp, type, variable) %>%  # unite conflict type with both number and rate
    tidyr::spread(temp, value) %>%  # make wide
    dplyr::rename(Silly = silly, 
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
    tidyr::gather(type, number, -locus) %>%  # make tall
    dplyr::mutate(qc_genotypes = sum(QCColSizeAll)) %>%  # join number of QC genotypes by locus
    dplyr::mutate(rate = number / qc_genotypes) %>%  # conflict numbers to rates
    tidyr::gather(variable, value, -locus, -qc_genotypes, -type) %>%  # make tall
    tidyr::unite(temp, type, variable) %>%  # unite conflict type with both number and rate
    tidyr::spread(temp, value) %>%  # make wide
    dplyr::rename(Locus = locus,
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
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Append Summary Tables to QCSummaryfile.xlsx ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create list object for all output to simple summary.xlsx file
  summary_lst <- suppressWarnings(list("Summary by Silly" = summary_table_1,
                      "Conflicts by Silly" = summary_table_2,
                      "Conflicts by Locus" = summary_table_3,
                      "Conflicts by PlateID" = conflicts_by_plate,
                      "QC Duplicate Check" = dup_check_results,
                      "Failure Rate by Silly" = failure_rate$silly_failure_rate,
                      "Failure Rate by Locus" = failure_rate$locus_failure_rate,
                      "Failure Rate by Plate" = failure_rate$plate_failure_rate,
                      "Overall Failure Rate" = failure_rate$overall_failure_rate,
                      "Project Sample Size by Locus" = OriginalProjectSampleSizebyLocus %>% tibble::rownames_to_column("silly") %>% tibble::as_tibble(),
                      "Duplicate Check in Project" = dplyr::bind_rows(DuplicateCheckReportSummary[!DuplicateCheckReportSummary == "No Duplicates"], .id = "silly") %>% tibble::as_tibble(),
                      "Alternate Species" = Alternate))
  
  # Write out a "Simple" file, can't update normal Summary File by inserting new tabs
  if(file.exists(QCSummaryfile)) {
    stop(paste0("QC Summary file: '", QCSummaryfile ,"' already exists, change the file name so you don't overwrite it, hoser!!!"))
  } else {
    writexl::write_xlsx(x = summary_lst, path = QCSummaryfile, col_names = TRUE)
  }

  #~~~  STOP!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  ##Don't Run#
  ############
}###########
############
############