PairwiseFstTree.GCL = function(sillyvec, loci, dir, nboots=1000, ncores = 4, returnbootstrapFst = FALSE){
  #############################################################################
  #
  # This function generates a matrix of pairwise Fst values, a neighbor joining
  # tree, bootstrap values for tree nodes, and variance components. This
  # function is able to multicore in order to speed up the calculation of 
  # variance components, with the idea that this function can be run on the
  # servers to take advantage of 16 cores. As this function requires hierfstat
  # it generates a .dat file.
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # sillyvec = KMA211Pops
  #   ~ A character vector of sillys.
  # loci = loci42
  #   ~ A character vector of loci. Can include haploid markers and/or combined loci.
  # dir = "Trees"
  #   ~ Directory where the fstat .dat file is dumped.
  # nboots = 1000
  #   ~ A numeric vector of length one indicating the number of bootstraps.
  # ncores = 4
  #   ~ A numeric vector of length one indicating the number of cores to use.
  # returnbootstrapFst = FALSE
  #   ~ A logical vector of length one indicating whether an Fst matrix should be
  #     saved for each bootstrap iteration.
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PairwiseFstTree: List of 4
  #  tree: List of 4
  #   edge: Numeric matrix
  #   edge.length: Numeric vector
  #   tip.label: Character vector of pop names
  #   Nnode: Integer
  #  bootstrap: Numeric vector of node bootstrap values
  #  PairwiseFst: Numeric matrix (length(sillyvec) x length(sillyvec)) of 
  #               pairwise Fst values
  #  vc: List of choose(n = length(sillyvec), k = 2) with variance components
  #      for each pair of sillys
  #
  # "fstatfile.dat" is put into "dir"
  #
  # PairwiseFstTree is 'dput' into "dir", named 
  # "paste(dir,"\\", length(sillyvec), "Pops", length(loci), "Loci_","PairwiseFstTree.txt",sep="")"
  #
  # Created by Jim Jasper on unknown???
  # Modified by Kyle Shedd to multicore
  #
  #############################################################################
  
  if(!exists("LocusControl")) {
    stop("'LocusControl' is required and not found, please create.")
  }
  
  while(!require("hierfstat")){install.packages("hierfstat")}
  
  while(!require("ape")){install.packages("ape")}

  while(!require("foreach")){install.packages("foreach")}
  
  while(!require("doParallel")){install.packages("doParallel")}
  
  while(!require("parallel")){install.packages("parallel")}
  
  if(ncores > detectCores()) {
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
  }
  
  message("\n4 Main function tasks: \n1) Create fstat .dat file\n2) Calculate variance componenets for each pair of sillys\n3) Calculate bootstrap Fst values\n4) Bootstrap tree nodes\n")
  if (.Platform$OS.type == "windows") flush.console()
  
  start.time <- Sys.time() 
  
  fstatdir=paste(dir,"\\fstatfile.dat",sep="")  
  
  cl=makePSOCKcluster(ncores)
  
  registerDoParallel(cl,cores=ncores)  
  
  nsillys=length(sillyvec)

  maxsillychar=nchar(nsillys)+1

  nloci=length(loci)

  ploidy=LocusControl$ploidy[loci]

  nalleles=LocusControl$nalleles[loci]

  maxchar=nchar(nalleles)+1
  names(maxchar)=loci

  alleles=LocusControl$alleles[loci]

  my.gcl=lapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)})
  names(my.gcl)=sillyvec

  n=sapply(sillyvec,function(silly){my.gcl[[silly]]$n})
  names(n)=sillyvec

  scores=lapply(sillyvec,function(silly){scores0 <- my.gcl[[silly]]$scores[,loci,] ; scores0[scores0 %in% c("Unk", "XXX", "0")] <- NA ; scores0})
  names(scores)=sillyvec

  counts <- foreach(silly=sillyvec) %dopar% {
    sapply(1:n[silly],function(i){
      paste(c(match(silly,sillyvec),sapply(loci,function(locus){
        ifelse(is.na(scores[[silly]][i,locus,1]),paste(rep(0,ploidy[locus]*maxchar[locus]),collapse=""),paste(sapply(1:ploidy[locus],function(allele){
          paste(c(rep(0,maxchar[locus]-nchar(match(scores[[silly]][i,locus,allele],alleles[[locus]]))),match(scores[[silly]][i,locus,allele],alleles[[locus]])),collapse="")
        }),collapse=""))
      })),collapse=" ")
    })
  }
  names(counts)=sillyvec

  fstat=paste(nsillys,nloci,max(nalleles),max(maxchar),sep=" ")

  fstat=rbind(fstat,cbind(loci))

  fstat=rbind(fstat,cbind(as.vector(unlist(counts))))

  write.table(fstat,fstatdir,row.names=FALSE,col.names=FALSE,quote=FALSE) 
  
  message("\nfstat file created\n", sep = '')
  
  dat=read.fstat.data(fstatdir)
  
  pairs=combn(sillyvec,2)

  pairnames=apply(pairs,2,function(col){paste(col,collapse=".")})
  
  dimnames(pairs)[[2]]=pairnames
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Calculate variance components for each pair of sillys
  message("\nCalculate variance components for each pair of sillys\n", sep = '')

  vc <- foreach(pair=pairnames, .packages = "hierfstat") %dopar% {
    sillys=pairs[,pair]
    sillyIND=sapply(sillys,function(silly){match(silly,sillyvec)})
    mydat=dat[dat[,1]==sillyIND[1] | dat[,1]==sillyIND[2],loci]
    myn=n[sillys]
    levels=rep(1:2,myn)
    
    mito.loci <- which(ploidy[loci] == 1)
    
    # Create matrix of variance components with order loci
    if(length(mito.loci) == 0) {
      vc.pair <- 
        # Diploid loci
        t(sapply(loci, function(locus) {
          mydata = data.frame(levels,mydat[,locus])
          hierfstat::varcomp(data = data.matrix(mydata), diploid=TRUE)$overall
        } ))
    } else {
      vc.pair <- rbind(
      # Diploid loci
      t(sapply(loci[-mito.loci], function(locus) {
        mydata = data.frame(levels,mydat[,locus])
        hierfstat::varcomp(data = data.matrix(mydata), diploid=TRUE)$overall
      } )),
      
      # Haploid loci
      t(sapply(loci[mito.loci], function(locus) {
        mydata = data.frame(levels,mydat[,locus])
        c(hierfstat::varcomp(data = data.matrix(mydata), diploid=FALSE)$overall,0)
      } ))
    )[loci, ]
    }
    
    # Replace NAs with 0s
    vc.pair[is.na(vc.pair)] <- 0
    vc.pair
  }
  names(vc) <- pairnames
  
  # Calculate Fst from VC for each pair
  Fst.ls <- lapply(vc, function(pair) {
    sum(pair[, 1]) / sum(pair)
  })
  
  
  # Create matrix of pairwise Fst
  Fst <- array(0,c(nsillys,nsillys),dimnames=list(sillyvec,sillyvec))
  
  Fst[lower.tri(Fst)] <- unlist(Fst.ls)
  Fst[upper.tri(Fst)] <- t(Fst)[upper.tri(Fst)]
  
  # Create tree
  tree=ape::nj(Fst)

  # Bootstraps
  message("\nCalculate bootstrap Fst values\n", sep = '')
  
  trees.bootstrapFst <- foreach(i=seq(nboots), .packages = "ape") %dopar% {
    
    temploci=sample(loci,nloci,replace=TRUE)
    tempFst.ls <- lapply(vc, function(pair) {
      sum(pair[temploci, 1]) / sum(pair[temploci, 1:3])
    })
    tempFst <- array(0,c(nsillys,nsillys),dimnames=list(sillyvec,sillyvec))
    tempFst[lower.tri(tempFst)] <- unlist(tempFst.ls)
    tempFst[upper.tri(tempFst)] <- t(Fst)[upper.tri(tempFst)]
    
    list(trees = ape::nj(tempFst), bootstrapFst = tempFst)
  }

  bootstrapFst <- lapply(trees.bootstrapFst, function(i) {
    i[["bootstrapFst"]]
  })
  
  trees <- lapply(trees.bootstrapFst, function(i) {
    i[["trees"]]
  })
  
  message("\nBootstrap tree nodes\n", sep = '')
  
  bootstrap=prop.clades(tree,trees)

  # Create final PairwiseFstTree list object
  if(returnbootstrapFst) {
    PairwiseFstTree <- list(tree=tree,bootstrap=bootstrap,PairwiseFst=Fst,vc=vc,BootstrapFst=bootstrapFst)
  } else {
    PairwiseFstTree <- list(tree=tree,bootstrap=bootstrap,PairwiseFst=Fst,vc=vc)
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Save tree before exiting
  dput(x = PairwiseFstTree, file = paste(dir,"\\", length(sillyvec), "Pops", length(loci), "Loci_","PairwiseFstTree.txt",sep="")  )
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime) 
  
  return(PairwiseFstTree)

}