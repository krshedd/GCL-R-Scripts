write_snps_gsi_sim.GCL <- function(collections, loci, path.name) {
# 
## write mixture file
#
# load("Z:\\Analysis\\2_Central\\Chinook\\Cook Inlet\\2014\\Mixture\\2015 ESSN Chinook postseason analysis\\ESSN2015PostseasonAnalysis.RData")
#
# path.name <- "C:/Users/jjasper/Documents/R/gsi_sim/test_data/KasilofEarly2015.txt"
#
# collections <- Mixtures[1]
#
# loci <- loci
#
## write baseline file
#
# attach("V:/Analysis/2_Central/Chinook/Cook Inlet/2014/Baseline/CIChinook2014Baseline.RData")
#
# path.name <- "C:/Users/jjasper/Documents/R/gsi_sim/test_data/CIbaselineChinook2014.txt"
#
# collections <- sillyvec

  while(!require(pryr)){install.packages("pryr")}

  filelines <- Reduce(rbind,lapply(collections, function(collection){cbind(c(paste("POP", collection), paste(get(paste0(collection, ".gcl"))$attributes$SillySource, apply(sapply(loci, function(locus){gsub(apply(apply(get(paste0(collection, ".gcl"))$scores[, locus, seq(LocusControl$ploidy[locus]), drop = FALSE], 1, as.numeric %.% factor, levels = c("A", "C", "G", "T", "-")), 2, paste, collapse = " "), pattern = "NA", replacement = "0")}), 1, paste, collapse = "  "), sep="   ")))}))

  write.table(paste(length(filelines) - length(collections), length(loci)), path.name, row.names = FALSE, col.names = FALSE, quote = FALSE, append=FALSE)

  write.table(cbind(loci), path.name, row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)

  write.table(filelines, path.name, row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)

  return(NULL)

}
