

  
if(FALSE){
  

  while(!require(reshape)){install.packages("reshape")}

  while(!require(abind)){install.packages("abind")}

# nchains <- 2

# NSIMS <- 10

# burn <- NSIMS/2

# thin <- 1

# wd <- "C:/Users/jjasper/Documents" ### "V:/Analysis/Staff/Jim Jasper/Sockeye/SEAK/MarkEnhancedGSI/MarkEnhancedGSIwithAge/All 2011-2014" ### 



# load("V:\\Analysis\\1_SEAK\\Sockeye\\Mixture\\Lynn Canal Inseason\\2015\\LynnCanalInseason2015.NewChilkatBaseline.RData")

# attach("V:\\Analysis\\1_SEAK\\Sockeye\\Mixture\\Lynn Canal Inseason\\2015\\PostSeason/LynnCanal_PostSeason_2015.RData")
# attach("V:/Analysis/1_SEAK/Sockeye/Baseline/2016/SEAK update for Lynn Canal Inseason/SEAK Baseline Chilkat update.RData")
# attach("V:\\Analysis\\1_SEAK\\Sockeye\\Baseline\\2013_2014\\SRO\\SEAKbaseline2014.RData")


# sillyvecMix <- "SGILL15D15"

# loci <- loci48

# metadata_file <- "2015.LynnCanal.Metadata.txt"

# age_classes <- sort(apply(expand.grid(seq(0, 5),seq(1, 5)), 1 ,paste, collapse = ""))

# age_groups <- c(1, 1, 1, 1, 1, 6, 2, 3, 6, 6, 6, 6, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)

# attach("V:/Analysis/1_SEAK/Sockeye/Baseline/2011/Nov2011/BaselineUpdate_Nov2011.RData")

#  sillyvecBase <- dget("V:/Analysis/1_SEAK/Sockeye/Baseline/2016/SEAK update for Lynn Canal Inseason/Output/LynnCanalBasepops173.txt")


  setwd(wd)

  metadat <- read.table(metadata_file, header=TRUE, sep="\t")

  metadat <- cbind(silly_vial = paste(metadat$SILLY_CODE, metadat$FISH_ID, sep = "_"), metadat[, -seq(2)])

  names(loci) <- loci

  nloci <- length(loci)

  nalleles <- LocusControl$nalleles[loci]

  PoolCollections.GCL(collections = sillyvecMix, loci = loci)

  all.gcl <- get(paste(c(sillyvecMix, "gcl"), collapse="."))

  attributes0 <- all.gcl$attributes

  x0 <- all.gcl$counts[,  loci, ]

  rownames(x0) <- attributes0$SillySource
  
# Meta Data Massaging  ########################################################################################################################################################################################################################################################################################################################################################################################################################################################


  metadat$Ages <- sapply(metadat$Ages, function(age){paste0(rep(0, 2-nchar(age)), age)})

  AgeNames <- paste0("Age", age_classes)

  age_groups <- setNames(age_groups, AgeNames)

  age_groups_names <- c("ZeroX", "Age1.2", "Age1.3", "Age2.2", "Age2.3", "AgeOther")

  A <- length(age_groups_names)

  C <- length(age_classes)

  stat_weeks <- sort(unique(metadat$SW))

  metadat <- cbind(metadat, cbind(w = as.numeric(factor(metadat$SW,  levels=stat_weeks,  ordered=TRUE)), a = as.numeric(factor(metadat$Ages,  levels=age_classes,  ordered=TRUE))))

  W <- length(stat_weeks)
 
# Baseline Wild  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################


  y0 <- FreqPop.GCL(sillyvecBase, loci)

  y <- Reduce(cbind, lapply(loci, function(locus){y0[, locus, seq(nalleles[locus])]}))

  K <- length(sillyvecBase)

  rownames(metadat) <- metadat$silly_vial

# Analysis  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  x <- Reduce(cbind, lapply(loci, function(locus){x0[metadat$silly_vial, locus, seq(nalleles[locus])]}))

  x[is.na(x)] <- 0  

  chains <- paste0("Chain", seq(nchains))

  require(foreach)

  require(doParallel)

  cl <- makePSOCKcluster(nchains)

  registerDoParallel(cl, cores = nchains)

  chain_times=foreach(chain = chains)%dopar%{

    rdirich <- function(alpha0){ vec <- rgamma(length(alpha0), alpha0, 1) ; vec / sum(vec) }
  
    beta <- matrix(rep(1/nalleles, nalleles), nrow = nrow(y), ncol = ncol(y), byrow = TRUE, dimnames = dimnames(y))
  
    beta_prm <- y + beta
  
    lnq <- log(t(apply(beta_prm, 1, function(rw){Reduce(cbind, tapply(rw, list(rep(loci, nalleles)), rdirich, simplify = FALSE))})))
  
    genofreq <- exp(x%*%t(lnq))

    piPrior <- array(1/C, c(K, C), list(seq(K), seq(C)))
  
    pPrior <- array(1/(K+H), c(W, K), list(seq(W), seq(K)))
  
    a <- metadat$a
  
    i <- apply(genofreq, 1, function(freqvec){sample(K, 1, TRUE, freqvec)})    
  
    pi0 <- table(i, a)
  
    dmns <- sapply(dimnames(pi0), as.numeric)
  
    pi <- piPrior
  
    pi[dmns$i, dmns$a] <- pi0 + piPrior[dmns$i, dmns$a]
  
    pi <- t(apply(pi, 1, rdirich))
  
    p0 <- table(metadat$w, i = i)
  
    pdmns <- sapply(dimnames(p0), as.numeric) 
  
    p <- pPrior
  
    p[, pdmns$i] <- p0 + pPrior[, pdmns$i]
  
    p <- t(apply(p, 1, rdirich))
  
    a[is.na(metadat$a)] <- sapply(i[is.na(metadat$a)], function(ii){sample(C, 1, TRUE, pi[ii, ])})
  
    beg_time <- Sys.time()
  
    for(sim in seq(NSIMS)){
  
      pi0 <- table(i, a)
  
      dmns <- sapply(dimnames(pi0), as.numeric)
  
      pi <- piPrior
  
      pi[dmns$i, dmns$a] <- pi0 + piPrior[dmns$i, dmns$a]
  
      pi <- t(apply(pi, 1, rdirich))
  
      p0 <- table(metadat$w, i)
  
      pdmns <- sapply(dimnames(p0), as.numeric) 
  
      p <- pPrior
  
      p[, pdmns$i] <- p0 + pPrior[, pdmns$i] 
  
      p <- t(apply(p, 1, rdirich))
  
      a[is.na(metadat$a)] <- sapply(i[is.na(metadat$a)], function(Z){sample(C, 1, TRUE, pi[Z, ])})
  
      a.df <- cbind(metadat$w, a)
  
      i <- apply(cbind(a.df, genofreq), 1, function(wa){sample(K, 1, TRUE, p[wa[1], seq(K)]*pi[seq(K), wa[2]]*wa[-seq(2)])}) 
  
      x_sum <- rowsum(x, group = i, reorder = TRUE)
  
      x_sum_nms <- as.numeric(rownames(x_sum))
  
      beta_prm_prm <- beta_prm 
  
      beta_prm_prm[x_sum_nms,] <- beta_prm[x_sum_nms, ] + x_sum
  
      lnq <- log(t(apply(beta_prm_prm, 1, function(rw){unlist(tapply(rw, INDEX = list(rep(loci, nalleles[loci])), FUN = rdirich))})))
  
      genofreq <- exp(x%*%t(lnq))
  
      if( sim > burn & ! sim %% thin){
  
        empty <- lapply(seq(W), function(w){ cat(c(format(p[w, ], trim = TRUE, digits = 16, scientific = TRUE), "\n"), file = paste0("P_Week", w, "_Chain", chain, ".txt")) })

        Pi <- t(rowsum(t(pi), group = age_groups))
  
        empty <- lapply(seq(A), function(aa){ cat(format(Pi[, aa], trim = TRUE, digits = 16, scientific = TRUE), "\n", file = paste0("Pi_", age_groups_names[aa], "_Chain", chain, ".txt")) })
  
      }# end if  
      
    }#sim
  
    Sys.time() - beg_time
  
  }#chain

}