create_rubias_mixture <- function(sillyvec, loci) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function creates the mixture dataframe needed for `rubias`.
  # It reformats the "scores" from each individual in the mixture(s) into a two column format used by `rubias`.
  # Each silly is treated as its own mixture.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   sillyvec - character vector of mixture sillys, each silly is treated as its own mixture
  #   loci - character vector of the loci you wish to include
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns a dataframe in `rubias` mixture format
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # load("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_rubias.RData")
  # chignik_2017.rubias_mix <- create_rubias_mixture(sillyvec = paste0("SCHIG17_Strata", 1:6), loci = loci22)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  silly_mix.lst <- lapply(sillyvec, function(silly) {
    my.gcl <- get(paste0(silly, ".gcl"))
    scores.mat <- t(apply(my.gcl$scores[, loci, ], 1, function(ind) {c(t(ind))} ))
    colnames(scores.mat) <- as.vector(sapply(loci, function(locus) {c(locus, paste(locus, 1, sep = "."))} ))
    scores.df <- data.frame(scores.mat, stringsAsFactors = FALSE)
    scores.df$sample_type <- "mixture"
    scores.df$repunit <- NA
    mode(scores.df$repunit) <- "character"
    scores.df$collection <- silly
    scores.df$indiv <- as.character(my.gcl$attributes$SillySource)
    silly_mix.df <- scores.df[, c("sample_type", "repunit", "collection", "indiv", gsub(pattern = "-", replacement = ".", x = colnames(scores.mat)))] } #silly
  )
  return(do.call("rbind", silly_mix.lst))
}