ConfusionMatricesLikeProfile.GCL <- function(likeprof){
  ####################################################################################################################################################################################################################################################################
  #
  # This function takes the genetic likelihood profile from LeaveOneOutLikeProfile.GCL
  # and summarizes into matrices at the group, group by pop, and pop by pop level.
  # This genetic likelihood profile is intended to show how well populations / 
  # reporting groups are differentiated by the markers at hand. These probabilities
  # can give an idea of potential misallocation and where.
  # This function can aid in marker selection.
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # likeprof = 
  #   ~ Unmodified list output from LeaveOneOutLikeProfile.GCL.
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Named List of 3
  #  GroupByGroup matrix (number of RGs x number of RGs)
  #   average probability of each individual from a group originating from a pop in that group
  #  GroupByPop matrix (number of pops x number of RGs)
  #   average probability of each individual from a pop originating from a pop in that group
  #  PopByPop matrix (number of RGs x number of RGs)
  #   average probability of each individual from a pop originating from a pop
  #
  # Originally created by Jim Jasper on to handle output from LeaveOneOutDist.GCL
  # Modified by Kyle Shedd to handle output from LeaveOneOutLikeProfile.GCL
  #
  ####################################################################################################################################################################################################################################################################
  
  sillyvec <- colnames(likeprof$IndividualByPop)
  
  groupnames <- colnames(likeprof$IndividualByGroup)
  
  groupvec <- as.numeric(sapply(sillyvec, function(silly) {which(groupnames == unique(subset(x = likeprof$Attributes, subset = likeprof$Attributes$FromPop == silly)$FromGroup))}))
  
  C <- length(groupvec)
  
  G <- max(groupvec)
  
  Cg <- as.vector(table(groupvec))
  
  cmCg <- c(0,cumsum(Cg)) 
  
  
  ind2popvec <- factor(x = likeprof$Attributes$FromPop, levels = unique(likeprof$Attributes$FromPop))
  ConfusionMatrixPopByPop <- t(sapply(seq(C), function(pop) {colMeans(likeprof$IndividualByPop[which(as.numeric(ind2popvec) == pop), ], na.rm = TRUE)}))
  dimnames(ConfusionMatrixPopByPop) <- list(paste(groupnames[groupvec], sillyvec, sep = "."),
                                            paste(groupnames[groupvec], sillyvec, sep = "."))
  
  ConfusionMatrixGroupByPop <- t(sapply(seq(C), function(pop) {colMeans(likeprof$IndividualByGroup[which(as.numeric(ind2popvec) == pop), ], na.rm = TRUE)}))
  dimnames(ConfusionMatrixGroupByPop) <- list(paste(groupnames[groupvec], sillyvec, sep = "."),
                                              groupnames)

  ConfusionMatrixGroupByGroup <- apply(ConfusionMatrixGroupByPop, 2, function(clm) {sapply(seq(G), function(g) {mean(clm[groupvec == g]) }) })
  dimnames(ConfusionMatrixGroupByGroup) <- list(groupnames, 
                                                groupnames)
  
  return(list(GroupByGroup = ConfusionMatrixGroupByGroup, GroupByPop = ConfusionMatrixGroupByPop, PopByPop = ConfusionMatrixPopByPop))
  
}