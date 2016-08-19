

  
if(FALSE){
  
# nchains <- 6

# NSIMS <- 40000

# burn <- NSIMS/2

# thin <- 10

# wd <- "C:/Users/jjasper/Documents/LynnCannal" 



# load("V:/Analysis/1_SEAK/Sockeye\\Mixture/Lynn Canal Inseason\\2015/PostSeason/LynnCanal_PostSeason_2015.RData")

# attach("V:\\Analysis\\1_SEAK\\Sockeye\\Baseline\\2013_2014\\SRO\\SEAKbaseline2014.RData")
# attach("V:\\Analysis\\1_SEAK\\Sockeye\\Baseline\\2016\\SEAK update for Lynn Canal Inseason\\SEAK Baseline Chilkat update.RData")

# 

# sillyvecMix <- "SGILL15D15"

# loci <- loci48

# metadata_file <- "2015.LynnCanal.Metadata.txt"

# age_classes <- sort(apply(expand.grid(seq(0, 5),seq(1, 5)), 1 ,paste, collapse = ""))

# age_groups <- c(1, 1, 1, 1, 1, 6, 2, 3, 6, 6, 6, 6, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)

# age_groups_names <- c("ZeroX", "Age1.2", "Age1.3", "Age2.2", "Age2.3", "AgeOther")

# sillyvecBase <- dget("V:/Analysis/1_SEAK/Sockeye/Baseline/2016/SEAK update for Lynn Canal Inseason/Output/LynnCanalBasepops173.txt")

# groups <- as.numeric(dget("V:\\Analysis\\1_SEAK\\Sockeye\\Mixture\\Lynn Canal Inseason\\2016\\PostSeason\\NewPostSeasonGroupvec.txt"))

# group_names <- c("ChilkatLake","ChilkatMain","Chilkoot","NSEAK","Snett","Taku/StikMain","TakuLakes","Other")



  while(!require(coda)){install.packages("coda")}

  setwd(wd)

  metadat <- read.table(metadata_file, header=TRUE, sep="\t")

  metadat <- cbind(silly_vial = paste(metadat$SILLY_CODE, metadat$FISH_ID, sep = "_"), metadat[, -seq(2)])

  names(loci) <- loci

  nloci <- length(loci)

  nalleles <- LocusControl$nalleles[loci]

  invisible(PoolCollections.GCL(collections = sillyvecMix, loci = loci))

  all.gcl <- get(paste(c(sillyvecMix, "gcl"), collapse="."))

  x0 <- all.gcl$counts[,  loci, ]

  rownames(x0) <- all.gcl$attributes$SillySource
  
# Meta Data Massaging  ########################################################################################################################################################################################################################################################################################################################################################################################################################################################


  metadat$Ages <- sapply(metadat$Ages, function(age){paste0(rep(0, 2-nchar(age)), age)})

  AgeNames <- paste0("Age", age_classes)

  age_groups <- setNames(age_groups, AgeNames)

  A <- length(age_groups_names)

  C <- length(age_classes)

  stat_weeks <- sort(unique(metadat$SW))

  metadat <- cbind(metadat, cbind(w = as.numeric(factor(metadat$SW,  levels=stat_weeks,  ordered=TRUE)), a = as.numeric(factor(metadat$Ages,  levels=age_classes,  ordered=TRUE))))

  W <- length(stat_weeks)
 
# Baseline Wild  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  y <- FreqPop.GCL(sillyvecBase, loci)

  y <- Reduce(cbind, lapply(loci, function(locus){y[, locus, seq(nalleles[locus])]}))

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

  beg_time <- Sys.time()

  chain_times <- foreach(chain = chains)%dopar%{

    rdirich <- function(alpha0){ vec <- rgamma(length(alpha0), alpha0, 1) ; vec / sum(vec) }
  
    beta <- matrix(rep(1/nalleles, nalleles), nrow = nrow(y), ncol = ncol(y), byrow = TRUE, dimnames = dimnames(y))
  
    beta_prm <- y + beta
  
    lnq <- log(t(apply(beta_prm, 1, function(rw){Reduce(cbind, tapply(rw, list(rep(loci, nalleles)), rdirich, simplify = FALSE))})))
  
    genofreq <- exp(x%*%t(lnq))

    piPrior <- matrix((1 / A / table(age_groups))[age_groups], nrow = K, ncol = C, byrow = TRUE)
  
    pPrior <- matrix((1 / max(groups) / table(groups))[groups], nrow = W, ncol = K, dimnames = list(seq(W), seq(K)))
  
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
  
        invisible(lapply(seq(W), function(ww){ cat(format(rowsum(p[ww, ], group = groups), trim = TRUE, digits = 16, scientific = TRUE), "\n", file = paste0("R_Week", ww, "_", chain, ".txt"), append = sim - thin > burn) }))

        Pi <- t(rowsum(t(pi), group = age_groups))

        PiR <- lapply(seq(W), function(ww){ rowsum(diag(p[ww, ]) %*% Pi, group = groups) })
        
        invisible(lapply(seq(W), function(ww){ lapply(seq(A), function(aa){ cat(format(PiR[[ww]][, aa], trim = TRUE, digits = 16, scientific = TRUE), "\n", file = paste0("PiR_Week", ww, "_", age_groups_names[aa], "_", chain, ".txt"), append = sim - thin > burn) }) }))
  
      }# end if  
      
    }#sim
     
  }#chain

  stopCluster(cl)

  tot_time <- Sys.time() - beg_time

  print(tot_time)

}



  R <- lapply(seq(W), function(ww){  lapply(chains, function(chain){ read.table(paste0("R_Week", ww, "_", chain, ".txt")) }) })

  R <- lapply(seq(W), function(ww){  as.mcmc.list(lapply(seq(chains), function(chain){ mcmc(R[[ww]][[chain]]) })) })

  GR_R <- lapply(lapply(lapply(R, gelman.diag, transform = TRUE, autoburnin = FALSE, multivariate = FALSE), "[[", 1), function(gr){ gr[, 1] })

  R <- lapply(seq(W), function(ww){ Reduce(rbind, R[[ww]]) }) 

  R <- lapply(seq(W), function(ww){ t(apply(R[[ww]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), GR = NA) })) }) 

  R <- setNames(lapply(seq(W), function(ww){ rownames(R[[ww]]) <- group_names ; R[[ww]][,"GR"] <- GR_R[[ww]] ; R[[ww]] }), stat_weeks)

  sink("R_all_weeks.txt")

  print(R)

  sink()


  PiR <- lapply(seq(W), function(ww){ lapply(seq(A), function(aa){ lapply(chains, function(chain){ read.table(paste0("PiR_Week", ww, "_", age_groups_names[aa], "_", chain, ".txt")) }) }) })

  GR_PiR <- lapply(seq(W), function(ww){ lapply(seq(A), function(aa){ gelman.diag(as.mcmc.list(lapply(PiR[[ww]][[aa]], mcmc)), transform = TRUE, autoburnin = FALSE, multivariate = FALSE)[[1]][, 1] }) }) 

  PiR <- lapply(seq(W), function(ww){ lapply(seq(A), function(aa){ mat <- Reduce(rbind, PiR[[ww]][[aa]])  ; colnames(mat) <- group_names ; mat }) })

  PiR <- lapply(seq(W), function(ww){ lapply(seq(A), function(aa){ t(apply(PiR[[ww]][[aa]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), GR = NA)} )) }) })

  PiR <- setNames(lapply(seq(W), function(ww){ setNames(lapply(seq(A), function(aa){ PiR[[ww]][[aa]][, "GR"] <- GR_PiR[[ww]][[aa]] ; PiR[[ww]][[aa]] }), age_groups_names) }), stat_weeks)
  
  PiR <- Reduce(rbind, lapply(as.character(stat_weeks), function(ww){ Reduce(rbind, lapply(age_groups_names, function(aa){ data.frame(Group = group_names, StatWeek = ww, AgeClass = aa, PiR[[ww]][[aa]])  })) }))

  write.table(PiR, "PiR.txt", col.names = NA, sep = "\t")

