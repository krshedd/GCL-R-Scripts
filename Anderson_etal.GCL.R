Anderson_etal.GCL <- function(popvec, loci, groups, group_names, group_comps = NULL) {
############################################################################################################################################################################################################################################
#
#  rm(list=ls(all=TRUE))
#
#  load("V:/Analysis/2_Central/Chinook/Kenai River/2011 Baseline/KenaiChinook2011.RData")
#
#  popvec <- KenaiNames10 ; loci <- loci40 ; group_names <- group_names ; groups <- groupvec3 ;  group_comps <- NULL
#
############################################################################################################################################################################################################################################

  names(groups) <- popvec

  if(is.null(group_comps)){ group_comps <- rep(1 / length(group_names), length(group_names)) }

  names(group_comps) <- group_names

  P <- setNames((group_comps / table(groups))[groups], popvec)

  nalleles <- LocusControl$nalleles[loci]

  names(popvec) <- popvec

  gcls <- setNames(lapply(paste0(popvec, ".gcl"), get), popvec)

  M <- sapply(gcls, "[[", "n")

  locus_fac <- factor(rep(loci, nalleles), levels = loci)

  alleles <- paste(locus_fac, unlist(lapply(nalleles, seq)), sep = "_")

  vials <- lapply(gcls, function(gcl){ as.character(gcl[["attributes"]]$SillySource) })

  x <- Reduce(rbind, lapply(popvec, function(i){ mat = Reduce(cbind, lapply(loci, function(locus){ gcls[[i]][["counts"]][, locus, seq(nalleles[locus])] })) ; dimnames(mat) = list(vials[[i]], alleles) ; mat[is.na(mat)] = 0 ; mat  }))

  y <- Reduce(cbind, lapply(loci, function(locus){ FreqPop.GCL(popvec, loci = locus)[,1,] }))

  colnames(y) <- alleles

  genefreq <- exp(Reduce(rbind, lapply(popvec, function(i){ Reduce(rbind, lapply(vials[[i]], function(m){ mat = (y + rep(1, length(popvec)) %*% rbind(1 / rep(nalleles, nalleles)) - (popvec == i) %*% rbind(x[m, ])) ; t_lnq = log(t(mat) / apply(rowsum(t(mat), locus_fac, reorder = FALSE), 2, rep, each = 2)) ; x[m, ] %*% t_lnq })) })))

  dimnames(genefreq) <- list(unlist(vials), popvec)

  pop_fac <- factor(rep(popvec, M), levels = popvec)

  group_fac <- factor(group_names[groups[pop_fac]], levels = group_names)

  rel_like <- genefreq / rowSums(genefreq)  # pop flat
  
  colnames(rel_like) <- paste0("Into_", popvec)
  
  pop_mean_rel_like <- sapply(paste0("Into_", popvec), function(p){ tapply(rel_like[,p], factor(paste0("From_", pop_fac), levels = paste0("From_", popvec)), mean) })
  
  rel_like <- genefreq %*% diag(P) / apply(genefreq %*% diag(P), 1, sum)  # group flat (if group_comps = NULL)
  
  group_rel_like <- t(rowsum(t(rel_like), factor(paste0("Into_", group_names[groups]), levels = paste0("Into_", group_names))))

  group_mean_rel_like <- sapply(paste0("Into_", group_names), function(g){ tapply(group_rel_like[,g], factor(paste0("From_", group_fac), levels = paste0("From_", group_names)), mean) })

  return(list(group_mean_rel_like = group_mean_rel_like, pop_mean_rel_like = pop_mean_rel_like, genefreq = data.frame(from_pop = pop_fac, from_group = group_fac, genefreq)))

}