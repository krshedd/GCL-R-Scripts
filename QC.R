##Don't Run#
############
if(FALSE){##
############
############

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Title: K102 QC
  #Date: Thu Sep 01 12:01:53 2016
  #Name: Heather Hoyt; Jim Jasper; Kyle Shedd
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
  # Overall Failure Rate
  # Original Project Sample Size by Locus
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Setup ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  date()
  # rm(list=ls(all=TRUE))
  
  # This sources all of the new GCL functions to this workspace
  # source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Arguments ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  dirQC <- "V:/Lab/Genotyping/SNP Projects/Sockeye/Project S161 Chignik Inseason 2016/QC"

  species <- "sockeye"

  markersuite <- "Sockeye2013Chignik_24SNPs"

  project <- "S161"

  projectID <- 2284

  username <- "krshedd"

  password <- ""

  QCSummaryfile <- "Project K102 QC Summary TEST.xlsx"
  
  conflict_rate <- 0.10  # conflict rate at which dupcheck between sillys occurs

#~~~  GO! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  while(!require(abind)){ install.packages("abind") }

  bbind <- function(...) { abind(..., along = 3) }

  while(!require(lattice)){ install.packages("lattice") }

  while(!require(xlsx)){ install.packages("xlsx") }

  source(path.expand("~/R/Functions.GCL.R"))

  setwd(dirQC)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Create Locus Control ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  CreateLocusControl.GCL(markersuite = markersuite, username = username, password = password)

  loci <- LocusControl$locusnames

  nalleles <- LocusControl$nalleles

  ploidy <- LocusControl$ploidy

  alleles <- LocusControl$alleles

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Read in Project Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ReadProjectLOKI2R.GCL(projectID = projectID, username = username, password = password)

  rm(password)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Failure Rate ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  FailureRate <- FailureRate.GCL(sillyvec = ProjectSillys)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Read in QC Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  QCfiles <- list.files(path = "Genotype Data Files", pattern = ".csv", full.names = TRUE, recursive = FALSE)

  if(max(nalleles) <= 2) {
    
    ReadBiomarkQC.GCL(QCcsvFilepaths = QCfiles)
    
  } else {
    
    QCdat <- Reduce(rbind, lapply(QCfiles, function(fle) {
      QCdat.temp <- read.csv(file = fle)[, 1:4]
      colnames(QCdat.temp) <- c("Sample.Name", "Marker", "Allele.1", "Allele.2")
      QCdat.temp} )
    )
    
    levels(QCdat$Marker) <- gsub(pattern = " ", replacement = "", x = levels(QCdat$Marker))
    
    vials <- levels(QCdat$Sample.Name)
    
    lociQC <- levels(QCdat$Marker)  
    
    QCdat <- data.frame(Silly = unlist(lapply(strsplit(as.character(QCdat$Sample.Name), split = "_"), "[", 1)), Vial = unlist(lapply(strsplit(as.character(QCdat$Sample.Name), split = "_"), "[", 2)), QCdat, stringsAsFactors = FALSE)
    
    QCdat$Vial <- as.character(as.numeric(QCdat$Vial))
    
    ProjectSillysQC <- unique(unlist(lapply(strsplit(as.character(QCdat$Sample.Name), split = "_"), "[", 1)))
    
    if(sum(! ProjectSillysQC %in% ProjectSillys)){ stop(paste0(ProjectSillysQC[! ProjectSillysQC %in% ProjectSillys], " not found in ProjectSillys.")) }
    
    if(sum(! lociQC %in% loci)){ stop(paste0(lociQC[! lociQC %in% loci], " not found in LocusControl.")) }
    
    attNames <- colnames(get(paste0(ProjectSillys[1], ".gcl"))$attributes)
    
    for(silly in ProjectSillysQC){
      
      sillydat <- subset(QCdat, Silly == silly)
      
      scores <- Reduce(bbind , lapply(c("Allele.1", "Allele.2"), function(al){ tapply(sillydat[, al], list(sillydat$Vial, sillydat$Marker), c)[,, drop = FALSE] }))
      
      counts <- array(NA, c(nrow(scores), ncol(scores), max(nalleles)), list(rownames(scores), colnames(scores), paste0("Allele", seq(max(nalleles)))))
      
      for(locus in loci){
        
        for(al in seq(nalleles[locus])){
          
          for(id in sillydat$Vial){
            
            counts[id, locus, al] <- sum(scores[id, locus, seq(ploidy[locus])] == alleles[[locus]][al])
            
          }#id
          
        }#al           
        
      }#locus
      
      attributes <- data.frame(matrix(NA, nrow = nrow(counts), ncol = length(attNames), dimnames = list(rownames(counts), attNames)))
      
      attributes$SillySource <- paste0(silly, "QC_", rownames(counts))
      
      assign(paste0(silly, "QC.gcl"), list(counts = counts, scores = scores, n = nrow(scores), attributes = attributes))
      
    }#silly
    
    QCSillys <- paste0(ProjectSillysQC, "QC")
    
  }#else for usat

  QCColSize <- sapply(paste(QCSillys, ".gcl", sep = ''), function(x) get(x)$n)
  
  QCColSizeAll <- setNames(rep(0, length(ProjectSillys)),paste0(ProjectSillys, "QC.gcl"))
  
  QCColSizeAll[paste0(QCSillys, ".gcl")] <- QCColSize[paste0(QCSillys, ".gcl")]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Read in Conflict Report ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  QCConcordanceReportfile <- list.files (path = "Conflict Reports", pattern = "ConcordanceReport", full.names = TRUE)

  CombineConflictsWithPlateID.GCL(files = QCConcordanceReportfile)

  if("Agreement" %in% levels(CombinedConflicts$Conflict)) {CombinedConflicts <- subset(CombinedConflicts, Conflict == "Conflict")}
  
  if("0" %in% levels(CombinedConflicts$Type)) { levels(CombinedConflicts$Type) <- gsub(pattern = "0", replacement = "Homo-Homo", levels(CombinedConflicts$Type)) }  # Old conflict report has "0" for mitochondrial conflicts
  
  if(" " %in% levels(CombinedConflicts$Type)) { levels(CombinedConflicts$Type)[levels(CombinedConflicts$Type) == " "] <- "Homo-Homo" }  # New conclict report has " " for mitochondrial conflicts
  
  QCtypes <- levels(CombinedConflicts$Type)

  types <- c("DB Zero", "File Zero", "Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo", "Conflict")

  ConflictsByPlateID <- matrix(data = 0, nrow = length(unique(CombinedConflicts$PlateID)), ncol = length(types), dimnames = list(sort(unique(CombinedConflicts$PlateID)), types))

  ConflictsByPlateID[, QCtypes] <- table(CombinedConflicts$PlateID, CombinedConflicts$Type)

  ConflictsByPlateID[, "Conflict"] = rowSums(ConflictsByPlateID[, c("Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo"), drop = FALSE])

  ConflictsBySilly <- matrix(data = 0, nrow = length(ProjectSillys), ncol = length(types), dimnames = list(ProjectSillys, types))

  ConflictsBySilly[, QCtypes] <- table(CombinedConflicts$Silly.Code, CombinedConflicts$Type)

  ConflictsBySilly[, "Conflict"] <- rowSums(ConflictsBySilly[, c("Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo"), drop = FALSE])

  ConflictsByLocus <- matrix(data = 0, nrow = length(unique(CombinedConflicts$Locus)), ncol = length(types), dimnames = list(sort(unique(CombinedConflicts$Locus)), types))

  ConflictsByLocus[, QCtypes] <- table(CombinedConflicts$Locus, CombinedConflicts$Type)

  ConflictsByLocus[, "Conflict"] <- rowSums(ConflictsByLocus[, c("Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo"), drop = FALSE])

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
  #### QC of Project Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ProjectSillys_SampleSizes <- matrix(data = NA, nrow = length(ProjectSillys), ncol = 5, dimnames = list(ProjectSillys, c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))

  ProjectSillys_SampleSizes[, "Genotyped"] <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n)

  if(species %in% c("chum", "sockeye")) {
  
    Alternate <- FindAlternateSpecies.GCL(sillyvec = ProjectSillys, species = species)
  
    nAltBySilly <- sapply(ProjectSillys, function(silly) {
      AlternateSpeciesReport <- Alternate[grep(pattern = silly, x = rownames(Alternate)), ]
      sum(AlternateSpeciesReport$Alternate > 0.5 & AlternateSpeciesReport$Failure > 0.5)
    })
    # RemoveAlternateSpecies.GCL(AlternateSpeciesReport = Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5)  # Do not remove fish, just note how many per silly. Still want to catch them in conflicts later.

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
  #### QC of QC Genotypes ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  MissLociQC <- RemoveIndMissLoci.GCL(sillyvec = QCSillys, proportion = 0.8)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Perform Duplicate Check on High Conflict Individuals ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  conflicts <- which(CombinedConflicts$Type %in% c("Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo"))

  conflictstab <- table(CombinedConflicts[conflicts, "Silly.Source"])

  hist(conflictstab / length(loci), main = "QC individual conflict rate", xlab = "Conflict rate", col = 8, breaks = seq(from = 0, to = 1, by = 0.02))

  conflict_bool <- conflictstab / length(loci) > conflict_rate

  if(sum(conflict_bool)){ 

    conflict_indv_bool <- table(CombinedConflicts[conflicts, "Silly.Source"]) / length(loci) > conflict_rate

    conflict_indv <- names(conflict_indv_bool)[conflict_indv_bool]
 
    conflict_indv_numconflicts <- table(CombinedConflicts[conflicts, "Silly.Source"])[conflict_indv_bool]   

    message(paste0("The following individuals have > ", conflict_rate * 100, "% loci with conflicts between project and QC:\n"), paste(conflict_indv, conflict_indv_numconflicts[conflict_indv], "conflicts", collapse = "\n"))
    
    new_conflict_indv <- NULL

    for(silly_indv in conflict_indv) {

      conflict_indv_split <- unlist(strsplit(x = silly_indv, split = "_"))

      id <- conflict_indv_split[2]

      silly <- conflict_indv_split[1]

      if(id %in% MissLociQC[[paste0(silly, "QC")]]) {message(paste0("\n", silly, "QC_", id, " does not have at least 80% loci genotyped, not running DupCheck for this individual."))}
      
      if(! id %in% get(paste0(silly, ".gcl"))$attributes$FK_FISH_ID) {message(paste0("\n", silly, "_", id, " does not have at least 80% loci genotyped, not running DupCheck for this individual."))}
      
      new_conflict_indv <- c(new_conflict_indv, paste(silly, id, sep = "_")[ ! id %in% MissLociQC[[paste0(silly, "QC")]] & id %in% get(paste0(silly, ".gcl"))$attributes$FK_FISH_ID ])  # Confirm QC fish and Project fish were not removed

    }#silly_indv

    original_conflict_indv <- conflict_indv

    conflict_indv <- new_conflict_indv

    conflict_silly <- unique(unlist(lapply(conflict_indv, function(ind) {strsplit(x = ind, split = "_")[[1]][1]} )))

    if(is.null(conflict_silly)) {
      
      message("\nNo remaining high conflict individuals.")
      
    } else {#conflict_silly
      
      message("\nRunning DupCheckBetweenSillys.GCL on these high conflict individuals, as they have at least 80% loci genotyped for Project and QC extractions.")
      
      message(paste(conflict_indv, conflict_indv_numconflicts[conflict_indv], "conflicts", collapse = "\n"))
      
      KeySillyIDs <- setNames(lapply(conflict_silly, function(silly) {sapply(grep(pattern = silly, x = conflict_indv, value = TRUE), function(ind) {unlist(strsplit(x = ind, split = paste(silly, "_", sep = '')))[2]}, USE.NAMES = FALSE) }), paste0(conflict_silly, "QC"))
      
      DupCheckResults <- setNames(lapply(conflict_silly, function(silly) {DupCheckBetweenSillys.GCL(KeySillys = paste0(silly, "QC"), KeySillyIDs = KeySillyIDs[paste0(silly, "QC")], BetweenSillys = ProjectSillys, loci = loci, threshold = 0.9)} ), nm = conflict_silly)
      
    }
    
  } else {#conflict_bool

    message(paste("No individuals have > ", conflict_rate * 100, "% loci with conflicts between project and QC.", sep = ''))
    
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Create Summary Tables ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  SummaryTable1 <- cbind(ProjectSillys_SampleSizes, "Failure Rate" = FailureRate$Silly_Failure_Rate, "Total QC Fish" = QCColSize)

  tab_names <- c("Total QC Genotypes", "Discrepancy Rate", "Total Het-Het", "Het-Het Rate", "Total Het-Homo", "Het-Homo Rate", "Total Homo-Het", "Homo-Het Rate", "Total Homo-Homo", "Homo-Homo Rate", "DB Zeros", "DB Zero Rate", "QC Zeros", "QC Zero Rate")

  SummaryTable2 <- matrix(NA, nrow = length(ProjectSillys), ncol = length(tab_names), dimnames = list(ProjectSillys, tab_names))

  SummaryTable2[, "Total QC Genotypes"] <- QCColSizeAll * length(loci)

  for(silly in ProjectSillys) {

    if(silly %in% rownames(ConflictsBySilly)){

      SummaryTable2[silly, c("Total Het-Het", "Total Het-Homo", "Total Homo-Het", "Total Homo-Homo", "DB Zeros", "QC Zeros")] <- ConflictsBySilly[silly, c("Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo", "DB Zero", "File Zero")]

      SummaryTable2[silly, "Discrepancy Rate"] <- ConflictsBySilly[silly, "Conflict"] / SummaryTable2[silly, "Total QC Genotypes"]

      SummaryTable2[silly, c("Het-Het Rate", "Het-Homo Rate", "Homo-Het Rate", "Homo-Homo Rate", "DB Zero Rate", "QC Zero Rate")] <- SummaryTable2[silly, c("Total Het-Het", "Total Het-Homo", "Total Homo-Het", "Total Homo-Homo", "DB Zeros", "QC Zeros")] / SummaryTable2[silly, "Total QC Genotypes"]

    } else {

      SummaryTable2[silly, -1] <- 0

    }

  }#silly
  
  SummaryTable2[is.na(SummaryTable2)] <- 0

  SummaryTable3 <- matrix(data = NA, nrow = length(loci), ncol = length(tab_names), dimnames = list(loci, tab_names))

  SummaryTable3[, "Total QC Genotypes"] <- sum(QCColSizeAll)

  for(locus in loci) {

    if(locus %in% rownames(ConflictsByLocus)){

      SummaryTable3[locus, c("Total Het-Het", "Total Het-Homo", "Total Homo-Het", "Total Homo-Homo", "DB Zeros", "QC Zeros")] <- ConflictsByLocus[locus, c("Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo", "DB Zero", "File Zero")]

      SummaryTable3[locus, "Discrepancy Rate"] <- ConflictsByLocus[locus, "Conflict"] / SummaryTable3[locus, "Total QC Genotypes"]

      SummaryTable3[locus, c("Het-Het Rate", "Het-Homo Rate", "Homo-Het Rate", "Homo-Homo Rate", "DB Zero Rate", "QC Zero Rate")] <- SummaryTable3[locus, c("Total Het-Het", "Total Het-Homo", "Total Homo-Het", "Total Homo-Homo", "DB Zeros", "QC Zeros")] / SummaryTable3[locus, "Total QC Genotypes"]

    } else {

      SummaryTable3[locus, -1] <- 0

    }

  }#locus

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Append Summary Tables to QCSummaryfile.xlsx ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  write.xlsx(x = SummaryTable1, file = QCSummaryfile, sheetName = "Summary by Silly", row.names = TRUE, col.names = TRUE, append = TRUE)

  write.xlsx(x = SummaryTable2, file = QCSummaryfile, sheetName = "Conflicts by Silly", row.names = TRUE, col.names = TRUE, append = TRUE)

  if(exists("DupCheckResults")) {

    print(DupCheckResults)
    
    invisible(lapply(conflict_silly, function(silly) {write.xlsx(x = DupCheckResults[[silly]], file = QCSummaryfile, sheetName = paste("DupCheckBetween", silly, sep = " "), row.names = TRUE, col.names = TRUE, append = TRUE)} ))

  }

  write.xlsx(x = SummaryTable3, file = QCSummaryfile, sheetName = "Conflicts by Locus", row.names = TRUE, col.names = TRUE, append = TRUE)

  write.xlsx(x = ConflictsByPlateID, file = QCSummaryfile, sheetName = "Conflicts by PlateID", row.names = TRUE, col.names = TRUE, append = TRUE)

  write.xlsx(x = sort(FailureRate$Silly_Failure_Rate, decreasing = TRUE), file = QCSummaryfile, sheetName = "Failure Rate by Silly", row.names = TRUE, col.names = TRUE, append = TRUE)

  write.xlsx(x = sort(FailureRate$Locus_Failure_Rate, decreasing = TRUE), file = QCSummaryfile, sheetName = "Failure Rate by Locus", row.names = TRUE, col.names = TRUE, append = TRUE)

  write.xlsx(x = FailureRate$Overall_Failure_Rate, file = QCSummaryfile, sheetName = "Overall Failure Rate", row.names = TRUE, col.names = TRUE, append = TRUE)

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










