##Don't Run#
############
if(FALSE){##
############
############


#~~~  Arguments  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  dirQC <- "V:/Lab/Genotyping/Microsatellite Projects/Chinook/Project K102 SEAK Origins Winter Troll 2016 Part 1/QC"

  species <- "chinook"

  markersuite <- "GAPS_Chinook_uSATs"

  project <- "K102"

  projectID <- 2281

  username <- "jjasper"

  password <- ""

  QCSummaryfile <- "Project K095 QC Summary.xlsx"

#~~~  GO! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  while(!require(abind)){ install.packages("abind") }

  bbind <- function(...) { abind(..., along = 3) }

  while(!require(lattice)){ install.packages("lattice") }

  while(!require(xlsx)){ install.packages("xlsx") }

  source(path.expand("~/R/Functions.GCL.R"))

  setwd(dirQC)

  CreateLocusControl.GCL(markersuite = markersuite, username = username, password = password)

  loci <- LocusControl$locusnames

  nalleles <- LocusControl$nalleles

  ploidy <- LocusControl$ploidy

  alleles <- LocusControl$alleles

  ReadProjectLOKI2R.GCL(projectID = projectID, username = username, password = password)

  rm(password)

  FailureRate <- FailureRate.GCL(sillyvec = ProjectSillys)

  QCfiles <- list.files(path = "Genotype Data Files", pattern = ".csv", full.names = TRUE, recursive = FALSE)

  QCdat <- Reduce(rbind, lapply(QCfiles, read.csv))

  vials <- levels(QCdat$Name)

  lociQC <- levels(QCdat$Assay)  

  QCdat <- data.frame(Silly = unlist(lapply(strsplit(as.character(QCdat$Name), split = "_"), "[", 1)), Vial = unlist(lapply(strsplit(as.character(QCdat$Name), split = "_"), "[", 2)), QCdat, stringsAsFactors = FALSE)

  ProjectSillysQC <- unique(unlist(lapply(strsplit(as.character(QCdat$Name), split = "_"), "[", 1)))

  if(sum(! ProjectSillysQC %in% ProjectSillys)){ stop(paste0(ProjectSillysQC[! ProjectSillysQC %in% ProjectSillys], " not found in ProjectSillys.")) }

  if(sum(! lociQC %in% loci)){ stop(paste0(lociQC[! lociQC %in% loci], " not found in LocusControl.")) }

  attNames <- colnames(get(paste0(ProjectSillys[1], ".gcl"))$attributes)

  for(silly in ProjectSillysQC){

    sillydat <- subset(QCdat, Silly == silly)

    scores <- Reduce(bbind , lapply(c("Allele.1", "Allele.2"), function(al){ tapply(sillydat[, al], list(sillydat$Vial, sillydat$Assay), c)[,, drop = FALSE] }))

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

  QCColSize <- sapply(paste(QCSillys, ".gcl", sep = ''), function(x) get(x)$n)

  QCColSizeAll <- setNames(rep(0, length(ProjectSillys)),paste0(ProjectSillys, "QC.gcl"))

  QCColSizeAll[paste0(QCSillys, ".gcl")] <- QCColSize[paste0(QCSillys, ".gcl")]

  QCConcordanceReportfile <- list.files ("Conflict Reports", pattern = "ConcordanceReport", full.names = TRUE)

  CombineConflictsWithPlateID.GCL(files = QCConcordanceReportfile)

  if("0" %in% CombinedConflicts$Type) { CombinedConflicts$Type <- gsub(pattern = "0", replacement = "Homo-Homo", CombinedConflicts$Type) }

  QCtypes <- levels(CombinedConflicts$Type)

  types <- union(QCtypes, c("DB Zero", "File Zero", "Het-Homo", "Homo-Het", "Homo-Homo", "Conflict"))

  ConflictsByPlateID <- matrix(data = 0, nrow = length(unique(CombinedConflicts$PlateID)), ncol = length(types), dimnames = list(unique(CombinedConflicts$PlateID), types))

  ConflictsByPlateID[, QCtypes] <- table(CombinedConflicts$PlateID, CombinedConflicts$Type)

  ConflictsByPlateID[, "Conflict"] = rowSums(ConflictsByPlateID[, c("Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo"), drop = FALSE])

  ConflictsBySilly <- matrix(data = 0, nrow = length(unique(CombinedConflicts$Silly.Code)), ncol = length(types), dimnames = list(unique(CombinedConflicts$Silly.Code), types))

  ConflictsBySilly[, QCtypes] <- table(CombinedConflicts$Silly.Code, CombinedConflicts$Type)

  ConflictsBySilly[, "Conflict"] <- rowSums(ConflictsBySilly[, c("Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo"), drop = FALSE])

  ConflictsByLocus <- matrix(data = 0, nrow = length(unique(CombinedConflicts$Locus)), ncol = length(types), dimnames = list(unique(CombinedConflicts$Locus), types))

  ConflictsByLocus[, QCtypes] <- table(CombinedConflicts$Locus, CombinedConflicts$Type)

  ConflictsByLocus[, "Conflict"] <- rowSums(ConflictsByLocus[, c("Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo"), drop = FALSE])

  OriginalProjectSampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = ProjectSillys, loci = loci)

  OriginalProjectPercentbyLocus <- apply(OriginalProjectSampleSizebyLocus, 1, function(row) {row / max(row)} )

  reruns <- which(apply(OriginalProjectPercentbyLocus, 2, min) < 0.8)  

  new_colors <- colorRampPalette(c("black", "white"))

  levelplot(t(OriginalProjectPercentbyLocus), col.regions = new_colors, at = seq(0, 1, length.out = 100), main = "% Genotyped", xlab = "SILLY", ylab = "Locus", scales = list(x = list(rot = 90)), aspect = "fill") # aspect = "iso" will make squares

  OriginalQCSampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = QCSillys, loci = loci)

  OriginalQCPercentbyLocus <- apply(OriginalQCSampleSizebyLocus, 1, function(row) {row / max(row)} )

  rerunsQC <- which(apply(OriginalQCPercentbyLocus, 2, min) < 0.8)

  levelplot(t(OriginalQCPercentbyLocus), col.regions = new_colors, at = seq(0, 1, length.out = 100), main = "% Genotyped", xlab = "SILLY", ylab = "Locus", scales = list(x = list(rot = 90)), aspect = "fill") # aspect = "iso" will make squares

  ProjectSillys_SampleSizes <- matrix(data = NA, nrow = length(ProjectSillys), ncol = 5, dimnames = list(ProjectSillys, c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))

  ProjectSillys_SampleSizes[, "Genotyped"] <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n)

  if(species == "chum" | species == "sockeye") {
  
    Alternate <- FindAlternateSpecies.GCL(sillyvec = ProjectSillys, species = species)
  
    RemoveAlternateSpecies.GCL(AlternateSpeciesReport = Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5)

  }

  ColSizePostAlternate <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n)

  ProjectSillys_SampleSizes[, "Alternate"] <- ProjectSillys_SampleSizes[, "Genotyped"] - ColSizePostAlternate 

  MissLoci <- RemoveIndMissLoci.GCL(sillyvec = ProjectSillys, proportion = 0.8)

  ColSizePostMissLoci <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n)

  ProjectSillys_SampleSizes[, "Missing"] <-  ColSizePostAlternate - ColSizePostMissLoci

  DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = ProjectSillys, loci = loci, quantile = NULL, minproportion = 0.95)

  DuplicateCheckReportSummary <- sapply(ProjectSillys, function(x) DuplicateCheck95MinProportion[[x]]$report)

  RemovedDups <- RemoveDups.GCL(DuplicateCheck95MinProportion)

  ColSizePostDuplicate <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n)

  ProjectSillys_SampleSizes[, "Duplicate"] <- ColSizePostMissLoci - ColSizePostDuplicate

  ProjectSillys_SampleSizes[, "Final"] <- ColSizePostDuplicate

  MissLociQC <- RemoveIndMissLoci.GCL(sillyvec = QCSillys, proportion = 0.8)

  conflicts <- which(CombinedConflicts$Type %in% c("Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo"))

  conflictstab <- table(CombinedConflicts[conflicts, "Silly.Source"])

  hist(conflictstab / length(loci), main = "QC individual conflict rate", xlab = "Conflict rate", col = 8, breaks = seq(from = 0, to = 0.5, by = 0.02))

  conflict_bool <- conflictstab / length(loci) > 0.10

  if(sum(conflict_bool)){ 

    conflict_indv_bool <- table(CombinedConflicts[conflicts, "Silly.Source"]) / length(loci) > 0.10

    conflict_indv <- names(conflict_indv_bool)[conflict_indv_bool]
 
    conflict_indv_numconflicts <- table(CombinedConflicts[conflicts, "Silly.Source"])[conflict_indv_bool]   

    new_conflict_indv <- NULL

    for(silly_indv in conflict_indv) {

      conflict_indv_split <- unlist(strsplit(x = silly_indv, split = "_"))

      id <- conflict_indv_split[2]

      silly <- conflict_indv_split[1]

      new_conflict_indv <- c(new_conflict_indv, paste(silly, id, sep = "_")[ ! id %in% MissLociQC[[paste0(silly, "QC")]] ])     

    }#silly_indv

    original_conflict_indv <- conflict_indv

    conflict_indv <- new_conflict_indv

    conflict_silly <- unique(unlist(lapply(conflict_indv, function(ind) {strsplit(x = ind, split = "_")[[1]][1]} )))

    message(paste("The following individuals have > 10% loci with conflicts between project and QC:\n"), paste(conflict_indv, conflict_indv_numconflicts[conflict_indv], "conflicts", collapse = "\n"))

    message(paste("Running DupCheckBetweenSillys.GCL on these SILLYs"))

    KeySillyIDs <- setNames(lapply(conflict_silly, function(silly) {sapply(grep(pattern = silly, x = conflict_indv, value = TRUE), function(ind) {unlist(strsplit(x = ind, split = paste(silly, "_", sep = '')))[2]}, USE.NAMES = FALSE) }), paste0(conflict_silly, "QC"))

    DupCheckResults <- setNames(lapply(conflict_silly, function(silly) {DupCheckBetweenSillys.GCL(KeySillys = paste0(silly, "QC"), KeySillyIDs = KeySillyIDs[paste0(silly, "QC")], BetweenSillys = ProjectSillys, loci = loci, threshold = 0.9)} ), nm = conflict_silly)

  }#conflict_bool

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

  write.xlsx(x = SummaryTable1, file = QCSummaryfile, sheetName = "Summary by Silly", row.names = TRUE, col.names = TRUE, append = TRUE)

  write.xlsx(x = SummaryTable2, file = QCSummaryfile, sheetName = "Conflicts by Silly", row.names = TRUE, col.names = TRUE, append = TRUE)

  if(exists("DupCheckResults")) {

    lapply(conflict.silly, function(silly) {write.xlsx(x = DupCheckResults[[silly]], file = QCSummaryfile, sheetName = paste("DupCheckBetween", silly, sep = " "), row.names = TRUE, col.names = TRUE, append = TRUE)} )

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

#~~~  STOP!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  

##Don't Run#
############
}###########
############
############










