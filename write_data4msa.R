write_data4msa <- function(mix_silly, base_sillys, loci){
##########################################################################################################################################################################################################################################
#
# wd <- "C:/Users/jjasper/Documents"
#
# setwd(wd)
#
# attach("V:\\Analysis\\2_Central\\Chinook\\Cook Inlet\\2014\\Mixture\\2016 ESSN Chinook postseason analysis\\ESSN2016PostseasonAnalysis.RData")
#
# attach("V:/Analysis/2_Central/Chinook/Cook Inlet/2014/Baseline/CIChinook2014Baseline.RData")
#
# mix_silly <- "KEASTSC16"
#
# base_sillys <- sillyvec
#
# loci <- loci39
#
# write_data4msa(mix_silly, base_sillys, loci)
#
# lapply(search()[2:3], detach, character.only = TRUE)
#
##########################################################################################################################################################################################################################################

  nalleles <- LocusControl$nalleles[loci]

  if(! missing(mix_silly)){

    mix <- get(paste0(mix_silly, ".gcl"))$counts[, loci, ]

    mix <- Reduce(cbind,lapply(loci, function(locus){ mix[,locus, seq(nalleles[locus])] }))

    mix[is.na(mix)] <- 0
 
    write.table(mix, "mix.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

  }

  if(! missing(base_sillys)){

    base <- FreqPop.GCL(base_sillys, loci) 

    base <- Reduce(cbind, lapply(loci, function(locus){base[,locus, seq(LocusControl$nalleles[locus])] }))

    write.table(base, "base.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

  }
  
  return(NULL)

}