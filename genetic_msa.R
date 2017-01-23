genetic_msa <- function(collection_mix, collections_base, loci, groupnames, groups, nchains = 8, nits = 2.5e4, nburn = nits / 2, thin = 1, group_prior, group_inits, out_dir, level = 0.1){

#
# load("V:/Analysis/2_Central/Chinook/Lower Cook Inlet/2016/Mixture/Sport Harvest 2014to2016/2014to2016LCISportFisheryMixtureAnalysis.RData")
#
# attach("V:/Analysis/2_Central/Chinook/Lower Cook Inlet/2015/Baseline/LowerCIChinook2015Baseline.RData")
#
# loci <- LocusControl$locusnames  
#
# collection_mix <- sillyvec_14_15[1] 
#
# collections_base <- PooledNames211 
#
# gcl2Nexus.GCL(sillyvec = collection_mix, loci = LocusControl$locusnames, path = getwd()) # NOT neccessary yet
#
# groupnames <- LCIgroups
#
# groups <- LCIgroupvec
#
# nchains <- max(groups)
#
# group_inits <- diag(rep(0.9, max(groups))) ; group_inits[group_inits == 0] <- 0.1 / max(groups)
#
# group_prior <- rep(1 / max(groups), max(groups)) 
#
# out_dir <- "C:/Users/jjasper/Documents"
#

  while(!require(coda)){ install.packages("coda") }

  while(!require(foreach)){ install.packages("foreach") }

  while(!require(doParallel)){ install.packages("doParallel") }

  setwd(out_dir)

  nloci <- length(loci)

  nalleles <- LocusControl$nalleles[loci]

  all.gcl <- get(paste0(collection_mix, ".gcl"))

  attributes0 <- all.gcl$attributes

  x0 <- all.gcl$counts[, loci, ]

  rownames(x0) <- attributes0$SillySource
  
  y <- FreqPop.GCL(collections_base, loci)

  y <- Reduce(cbind, lapply(loci,function(locus){y[,locus, seq(nalleles[locus])]}))

  K <- length(collections_base)

# Analysis  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  x <- Reduce(cbind, lapply(loci, function(locus){x0[, locus, seq(nalleles[locus])]}))

  x[is.na(x)] <- 0  

  rdirich <- function(alpha0){ vec <- rgamma(length(alpha0), alpha0, 1) ;  vec / sum(vec) }

  beta <- matrix(rep(1 / nalleles, nalleles), nrow = nrow(y), ncol = ncol(y), byrow = TRUE, dimnames = dimnames(y))

  beta_prm <- y + beta

  lnq <- log(t(apply(beta_prm, 1, function(rw){unlist(tapply(rw, INDEX = list(rep(loci, nalleles)), FUN = rdirich))})))

  chains <- paste0("Chain", seq(nchains))

  colnames(group_inits) <- chains 

  alpha <- (group_prior / table(groups))[groups]

  buff_size <- 100

  p_buff <- array(NA, c(buff_size, K))

  cl <- makePSOCKcluster(nchains)

  registerDoParallel(cl, cores = nchains)
  
  beg_time <- Sys.time()

  invisible(foreach(chain=chains)%dopar%{

    genofreq <- exp(x %*% t(lnq))

    p <- (group_inits[, chain] / table(groups))[groups]

    i <- apply(genofreq %*% diag(p), 1, function(freqvec){sample(K, 1, TRUE, freqvec)})

    for(sim in seq(nits)){

      p <- rdirich(tabulate(i, K) + alpha)
           
      i <- apply(genofreq %*% diag(p), 1, function(frq){sample(K, 1, TRUE, frq)})

      x_sum <- rowsum(x, group = collections_base[i], reorder = TRUE)

      x_sum_nms <- rownames(x_sum)

      beta_prm_prm <- beta_prm 

      beta_prm_prm[x_sum_nms,] <- beta_prm[x_sum_nms, ] + x_sum
 
      lnq <- log(t(apply(beta_prm_prm, 1, function(rw){unlist(tapply(rw, INDEX = list(rep(loci, nalleles)), FUN = rdirich))})))

      genofreq <- exp(x%*%t(lnq))

      if( sim > nburn & ! sim %% thin ){

        it <- ifelse(((sim - nburn) / thin) %% buff_size, ((sim - nburn) / thin) %% buff_size, buff_size) 

        p_buff[it, ] <- p 

        if(! ((sim - nburn) / thin) %% buff_size){ write.table(format(p_buff, digits = 16, scientific = TRUE), file = paste0(chain, "_P.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, append = ((sim - nburn) / thin) == buff_size) }

      }

    }#sim

  })#chain

  stopCluster(cl)
  
  assign("run_time", Sys.time() - beg_time, pos = 1)

  print(run_time)

  R <- as.mcmc.list(lapply(paste0(chains, "_P.txt"), function(Pfile){mcmc(t(rowsum(t(read.table(Pfile)), groups)))}))

  GR <- rep(NA, max(groups))

  if(nchains > 1){

    GR <- gelman.diag(R, transform = TRUE, autoburnin = FALSE, multivariate = FALSE)$psrf[,"Point est."]

  }

  R <- Reduce(rbind, lapply(R, cbind))

  R <- setNames(data.frame(setNames(apply(R, 2, mean), groupnames), apply(R, 2, sd), t(apply(R, 2, quantile, probs = c(0.5, level / 2, 1 - level / 2))), GR = GR), c("mean", "sd", "median", paste0("ci", c(level / 2, 1 - level / 2)))) 

  return(R)

}



