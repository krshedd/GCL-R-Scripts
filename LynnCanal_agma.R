

  
if(FALSE){
  
# nchains <- 6

# NSIMS <- 40000

# burn <- NSIMS/2

# thin <- 10

# wd <- "C:/Users/jjasper/Documents/LynnCanal" 



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

  while(!require(foreach)){install.packages("foreach")}

  while(!require(doParallel)){install.packages("doParallel")}

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

  metadat <- metadat[as.character(metadat$silly_vial) %in% rownames(x0), ]

  x <- Reduce(cbind, lapply(loci, function(locus){x0[as.character(metadat$silly_vial), locus, seq(nalleles[locus])]}))

  x[is.na(x)] <- 0  

  chains <- paste0("Chain", seq(nchains))

  cl <- makePSOCKcluster(nchains)

  registerDoParallel(cl, cores = nchains)

  beg_time <- Sys.time()

  invisible(foreach(chain = chains)%dopar%{ 

    rdirich <- function(alpha0){ vec <- rgamma(length(alpha0), alpha0, 1) ; vec / sum(vec) }
  
    beta <- matrix(rep(1/nalleles, nalleles), nrow = nrow(y), ncol = ncol(y), byrow = TRUE, dimnames = dimnames(y))
  
    beta_prm <- y + beta
  
    lnq <- log(t(apply(beta_prm, 1, function(rw){unlist(tapply(rw, INDEX = list(rep(loci, nalleles)), FUN = function(g){ g / sum(g) }))})))
  
    genofreq <- exp(x%*%t(lnq))
  
    pPrior <- matrix((1 / max(groups) / table(groups))[groups], nrow = W, ncol = K, byrow = TRUE, dimnames = list(seq(W), seq(K)))
  
    piPrior <- matrix((1 / A / table(age_groups))[age_groups], nrow = K, ncol = C, byrow = TRUE, dimnames = list(seq(K), seq(C)))

    a <- metadat$a
  
    p <- pPrior
  
    i <- apply(cbind(metadat$w, genofreq), 1, function(wfreq){sample(K, 1, TRUE, p[wfreq[1], seq(K)]*wfreq[-1])}) 
  
    for(sim in seq(NSIMS)){
  
      p0 <- table(metadat$w, i)
  
      nms <- colnames(p0) 
  
      p <- pPrior
  
      p[, nms] <- p0 + pPrior[, nms] 
  
      p <- t(apply(p, 1, rdirich))
  
      pi0 <- table(i, a)
  
      nms <- dimnames(pi0)
  
      pi <- piPrior
  
      pi[nms$i, nms$a] <- pi0 + piPrior[nms$i, nms$a]
  
      pi <- t(apply(pi, 1, rdirich))
  
      a[is.na(metadat$a)] <- sapply(i[is.na(metadat$a)], function(ii){sample(C, 1, TRUE, pi[ii, ])})
  
      i <- apply(cbind(cbind(metadat$w, a), genofreq), 1, function(wafreq){sample(K, 1, TRUE, p[wafreq[1], seq(K)]*pi[seq(K), wafreq[2]]*wafreq[-seq(2)])}) 
  
      x_sum <- rowsum(x, group = sillyvecBase[i], reorder = TRUE)
  
      nms <- rownames(x_sum)
  
      beta_prm_prm <- beta_prm 
  
      beta_prm_prm[nms,] <- beta_prm[nms, ] + x_sum
  
      lnq <- log(t(apply(beta_prm_prm, 1, function(rw){unlist(tapply(rw, INDEX = list(rep(loci, nalleles[loci])), FUN = rdirich))})))
 
      genofreq <- exp(x%*%t(lnq))
  
      if( sim > burn & ! sim %% thin){

        Pi <- t(rowsum(t(pi), group = age_groups))

        for(ww in seq(W)){

          write.table(rbind(format(rowsum(p[ww, ], group = groups)[,1], trim = TRUE, digits = 16, scientific = TRUE)), file = paste0("R_Week", ww, "_", chain, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, append = sim - thin > burn) 

          PiR <- rowsum(diag(p[ww, ]) %*% Pi, group = groups)

          for(aa in seq(A)){

            write.table(rbind(format(PiR[, aa], trim = TRUE, digits = 16, scientific = TRUE)), file = paste0("PiR_Week", ww, "_", age_groups_names[aa], "_", chain, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, append = sim - thin > burn)

          }#aa

        }#ww
  
      }# end if  
      
    }#sim

  })#chain

  stopCluster(cl)

  tot_time <- Sys.time() - beg_time

  print(tot_time)




  R <- lapply(seq(W), function(ww){  as.mcmc.list(lapply(chains, function(chain){ mcmc(read.table(paste0("R_Week", ww, "_", chain, ".txt"))) })) })

  GR_R <- lapply(lapply(lapply(R, gelman.diag, transform = TRUE, autoburnin = FALSE, multivariate = FALSE), "[[", 1), function(gr){ gr[, 1] })

  R <- lapply(seq(W), function(ww){ Reduce(rbind, R[[ww]]) }) 

  R <- lapply(seq(W), function(ww){ t(apply(R[[ww]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), GR = NA) })) }) 

  R <- setNames(lapply(seq(W), function(ww){ rownames(R[[ww]]) <- group_names ; R[[ww]][,"GR"] <- GR_R[[ww]] ; R[[ww]] }), stat_weeks)

  sink("R_all_weeks_with_update_and_age.txt")

  print(R)

  sink()


  PiR <- lapply(seq(W), function(ww){ lapply(seq(A), function(aa){ lapply(chains, function(chain){ read.table(paste0("PiR_Week", ww, "_", age_groups_names[aa], "_", chain, ".txt")) }) }) })

  GR_PiR <- lapply(seq(W), function(ww){ lapply(seq(A), function(aa){ gelman.diag(as.mcmc.list(lapply(PiR[[ww]][[aa]], mcmc)), transform = TRUE, autoburnin = FALSE, multivariate = FALSE)[[1]][, 1] }) }) 

  PiR <- lapply(seq(W), function(ww){ lapply(seq(A), function(aa){ mat <- Reduce(rbind, PiR[[ww]][[aa]])  ; colnames(mat) <- group_names ; mat }) })

  PiR <- lapply(seq(W), function(ww){ lapply(seq(A), function(aa){ t(apply(PiR[[ww]][[aa]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), GR = NA)} )) }) })

  PiR <- setNames(lapply(seq(W), function(ww){ setNames(lapply(seq(A), function(aa){ PiR[[ww]][[aa]][, "GR"] <- GR_PiR[[ww]][[aa]] ; PiR[[ww]][[aa]] }), age_groups_names) }), stat_weeks)
  
  PiR <- Reduce(rbind, lapply(as.character(stat_weeks), function(ww){ Reduce(rbind, lapply(age_groups_names, function(aa){ data.frame(Group = group_names, StatWeek = ww, AgeClass = aa, PiR[[ww]][[aa]])  })) }))

  write.table(PiR, "PiR2.txt", col.names = NA, sep = "\t")


}
