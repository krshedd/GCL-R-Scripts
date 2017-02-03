genetic_msa <- function(nalleles, groupnames, groups, group_prior, group_inits, nchains = 8, nits = 25000, nburn = 12500, thin = 1, q_out = FALSE, level = 0.1){
#
# wd <- "C:/Users/jjasper/Documents"
#
# setwd(wd)
#
# attach("V:/Analysis/2_Central/Chinook/Cook Inlet/2014/Baseline/CIChinook2014Baseline.RData")
#
# nalleles <- as.vector(LocusControl$nalleles[loci39])
#
# groups <- groupvec
#
# detach(search()[2], character.only = TRUE)
#
# groupnames <- c("West","Yentna","Susitna","Knik","Turnagain","Kenai","Kasilof")
#
# group_inits <- diag(0.9, nrow = max(groups), ncol = nchains) ; group_inits[group_inits == 0] <- 0.1 / (max(groups) - 1)
#
# group_prior <- rep(1 / max(groups), max(groups)) 
#

  while(!require(coda)){ install.packages("coda") }

  while(!require(foreach)){ install.packages("foreach") }

  while(!require(doParallel)){ install.packages("doParallel") }

  loc_fact <- factor(rep(paste0("Loc", seq(nalleles)), nalleles), levels = paste0("Loc", seq(nalleles)))

  x <- as.matrix(read.table("mix.txt"))
  
  y <- as.matrix(read.table("base.txt"))

  K <- nrow(y)

# Analysis  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  rdirich <- function(alpha0){ vec <- rgamma(length(alpha0), alpha0, 1) ;  vec / sum(vec) }

  beta <- matrix(rep(1 / nalleles, times = nalleles), nrow = nrow(y), ncol = ncol(y), byrow = TRUE)

  beta_prm <- y + beta

  t_lnq <- log(apply(beta_prm, 1, function(rw){unlist(tapply(rw, INDEX = loc_fact, FUN = rdirich))}))

  chains <- paste0("Chain", seq(nchains))

  colnames(group_inits) <- chains 

  alpha <- (group_prior / table(groups))[groups]

  buff_size <- 100

  cl <- makePSOCKcluster(nchains)

  registerDoParallel(cl, cores = nchains)
  
  beg_time <- Sys.time()

  invisible(foreach(chain = chains) %dopar% {

    genofreq <- exp(x %*% t_lnq)

    p_buff <- array(NA, c(buff_size, K))

    i_buff <- array(NA, c(buff_size, nrow(x)))

    if(q_out){ q_buff <- array(NA, c(K * buff_size, ncol(y))) }

    p <- (group_inits[, chain] / table(groups))[groups]

    i <- apply(genofreq %*% diag(p), 1, function(freqvec){sample(K, 1, TRUE, freqvec)})

    for(sim in seq(nits)){

      p <- rdirich(tabulate(i, K) + alpha)
           
      i <- apply(genofreq %*% diag(p), 1, function(frq){sample(K, 1, TRUE, frq)})

      x_sum <- rowsum(x, group = i, reorder = TRUE)

      x_sum_nms <- as.integer(rownames(x_sum))

      beta_prm_prm <- beta_prm 

      beta_prm_prm[x_sum_nms,] <- beta_prm[x_sum_nms, ] + x_sum
 
      t_lnq <- log(apply(beta_prm_prm, 1, function(rw){unlist(tapply(rw, INDEX = loc_fact, FUN = rdirich))}))

      genofreq <- exp(x %*% t_lnq)

      if( sim > nburn & ! sim %% thin ){

        it <- ((sim - nburn) / thin - 1) %% buff_size + 1 

        p_buff[it, ] <- p 

        i_buff[it, ] <- i

        if(q_out){ q_buff[seq(K * (it-1) + 1, K * it), ] <- t(exp(t_lnq)) }

        if(! ((sim - nburn) / thin) %% buff_size){ 

          write.table(format(p_buff, digits = 16, scientific = TRUE), file = paste0(chain, "_P.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, append = it > 1) 

          write.table(format(i_buff), file = paste0(chain, "_i.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, append = it > 1) 

          if(q_out){ write.table(format(q_buff, digits = 16, scientific = FALSE), file = paste0(chain, "_q.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, append = it > 1)  }

        }

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

  R <- setNames(data.frame(setNames(apply(R, 2, mean), groupnames), apply(R, 2, sd), t(apply(R, 2, quantile, probs = c(0.5, level / 2, 1 - level / 2))), GR = GR), c("mean", "sd", "median", paste0("ci", c(level / 2, 1 - level / 2)), "GR")) 

  return(R)

}



