
if(FALSE){

# ARGS  ########################################################################################################################################################################################################################################################################################################################################################################################################################################################

  wd <- ""

  nchains <- 6

  NSIMS <- 40000

  burn <- NSIMS/2

  thin <- 10

# load("Mixture_Workspace")

  sillyvecMix <- c("SGILL12D11", "SGILL12D6", "SGILL12D8", "SGILL13D11", "SGILL13D6", "SGILL13D8")

  meta.data.file <- "MetaDataAll.txt"

  agevec <- c(1, 1, 1, 1, 1, 6, 2, 3, 6, 6, 6, 6, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)

  AgeGroups <- c("ZeroX", "Age1.2", "Age1.3", "Age2.2", "Age2.3", "AgeOther")

# attach("Baseline_Workspace")

  sillyvecBase <- SEAK151pops

  loci <- LocusControl$locusnames[ ! LocusControl$locusnames %in% c("One_Cytb_26", "One_CO1", "One_Cytb_17", "One_c3-98", "One_GPDH", "One_MHC2_251")]

  # I am assuming that I need to read in the groups names sometime before this next line of code? Something like this:
  
  D106groups <- dget("V/...D106groups")
  D108groups <- dget("V/...D108groups")
  D111groups <- dget("V/...D111groups")
  
  # I brought in the groups and then ran this line and it bombed. It appeared to be missing parenthesis, so I added them and I got error
  groups_list <- lapply(list(D106groups, D108groups, D111groups), function(groups){ c(groups, rep(max(groups) + 1, length(Hatcheries)))})

  group_names_list <- list(c(D106group_names, "Hatcheries"), c(D108group_names, "Hatcheries"), c(D111group_names, "Hatcheries"))

  Hatcheries <- c("BURN", "MAIN", "MCDO", "NECK", "SPEE", "SWEE", "TAHL", "TATS", "TRAP", "TUYA")

##################################################################################################################################################################################################################################################################################################################################################################################################################################################################




  while(!require(reshape)){install.packages("reshape")}

  while(!require(coda)){install.packages("coda")}

  while(!require(foreach)){ install.packages("foreach") }

  while(!require(doParallel)){ install.packages("doParallel") }

  setwd(wd)

  metadat0 <- read.table(meta.data.file, header=TRUE, sep="\t")

  names(loci) <- loci

  nloci <- length(loci)

  nalleles <- LocusControl$nalleles[loci]

  invisible(PoolCollections.GCL(collections = sillyvecMix, loci = loci))

  x0 <- get(paste(c(sillyvecMix, "gcl"), collapse="."))$counts[,  loci, ]

  attributes0 <- get(paste(c(sillyvecMix, "gcl"), collapse="."))$attributes

  rownames(x0) <- attributes0$SillySource
  
# Meta Data Massaging  ########################################################################################################################################################################################################################################################################################################################################################################################################################################################

  years  <- sort(unique(metadat0$YEAR))

  T  <- length(years)

  YearData  <- data.frame(YEAR=years, t=seq(T))

  EuropeanAge <- sapply(metadat0$AGE_EUROPEAN, function(age){paste0(rep(0, 2-nchar(age)), age)})

  EuropeanAge[is.na(metadat0$AGE_EUROPEAN)] <- NA

  EuropeanAges <- sort(apply(expand.grid(seq(0, 5),seq(1, 5)), 1 ,paste, collapse = ""))
  
  AgeNames <- paste0("Age", EuropeanAges)

  agevec <- setNames(agevec, AgeNames)

  A <- length(AgeGroups)

  C <- length(EuropeanAges)

  AgeClass <- as.numeric(factor(EuropeanAge,  levels=EuropeanAges,  ordered=TRUE))

  metadatAux <- data.frame(

                  SillySource = metadat0$SillySource,

                  AgeClass = AgeClass,

                  WildGenetics = metadat0$SillySource %in% attributes0$SillySource & metadat0$ThermalGroupsForSEAK %in% "WILD",

                  WildNoGenetics = ! metadat0$SillySource %in% attributes0$SillySource & metadat0$ThermalGroupsForSEAK %in% "WILD",

                  NonWild = ! metadat0$ThermalGroupsForSEAK %in% "WILD"

               )

  metadat1 <- merge(metadat0, metadatAux)
 
  Districts <- sort(unique(metadat1$DISTRICT))

  D <- length(Districts)

  SubDistricts=sapply(Districts, function(d){as.vector(sort(unique(metadat1$SUBDISTRICT[metadat1$DISTRICT%in%d])))}, simplify=FALSE)

  S=unique(sapply(SubDistricts, length))

  index.spatial <- NULL 

  for(d in Districts){

    district <- match(d, Districts)

    index.spatial <- rbind(index.spatial, data.frame(DISTRICT=d, SUBDISTRICT=SubDistricts[[district]], d=district, s=seq(length(SubDistricts[[district]]))))

  }#d
  
  StatWeeks=sort(unique(metadat1$STAT_WEEK))

  W=length(StatWeeks)
 
  stat.week.dat=data.frame(STAT_WEEK=StatWeeks, w=seq(W))

  metadat2=merge(x = merge(x = merge(x = metadat1, y = YearData), y =  index.spatial), y = stat.week.dat)

# Baseline Wild  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  y0 <- FreqPop.GCL(sillyvecBase, loci)

  detach(search()[2], character.only = TRUE)

  y <- Reduce(cbind, lapply(loci, function(locus){y0[, locus, seq(nalleles[locus])]}))

  K <- length(sillyvecBase)

# Baseline Hatheries  ################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  H <- length(Hatcheries)

  Originmetadat <- data.frame(ThermalGroupsForSEAK = c(Hatcheries, "WILD"), z = c(seq(H) + K, NA))

  metadat3 <- merge(metadat2, Originmetadat)

  metadat <- sort_df(metadat3, vars=c("t", "NonWild", "WildNoGenetics", "WildGenetics", "d", "s", "w"))

  n <- tapply(metadat$WildGenetics, metadat[, c("t", "d", "s", "w")], length)

  n[is.na(n)] <- 0

  n0 <- as.data.frame.table(tapply(metadat$WildGenetics, metadat[, c("t", "d", "s", "w")], sum), responseName="WildGenetics")

  WildNoGenetics0 <- as.data.frame.table(tapply(metadat$WildNoGenetics, metadat[, c("t", "d", "s", "w")], sum), responseName="WildNoGenetics")$WildNoGenetics

  NonWild0 <- as.data.frame.table(tapply(metadat$NonWild, metadat[, c("t", "d", "s", "w")], sum), responseName="NonWild")$NonWild

  nn <- data.frame(n0, WildNoGenetics=WildNoGenetics0, NonWild=NonWild0)

  rownames(metadat) <- metadat$SillySource

  sillyvecBase <- c(sillyvecBase, Hatcheries)




# Analysis  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  metadat <- metadat[as.character(metadat$SillySource[metadat$WildGenetics]) %in% rownames(x0), ]

  x <- Reduce(cbind, lapply(loci, function(locus){x0[as.character(metadat$SillySource[metadat$WildGenetics]), locus, seq(nalleles[locus])]}))

  x[is.na(x)] <- 0  
  
  chains <- paste0("Chain", seq(nchains))

  cl <- makePSOCKcluster(nchains)

  registerDoParallel(cl, cores = nchains)

  beg_time <- Sys.time()
  
  invisible(foreach(chain = chains)%dopar%{

    rdirich <- function(alpha0){ vec <- rgamma(length(alpha0), alpha0, 1) ; vec / sum(vec) }
   
    beta <- matrix(rep(1/nalleles, nalleles), nrow = nrow(y), ncol = ncol(y), byrow = TRUE, dimnames = dimnames(y))
   
    beta_prm <- y + beta
   
    lnq <- log(t(apply(beta_prm, 1, function(rw){unlist(tapply(rw, list(rep(loci, nalleles)), rdirich, simplify = FALSE))})))
   
    genofreq <- exp(x%*%t(lnq))
   
    piPrior <- array(matrix((1 / A / table(agevec))[agevec], nrow = K + H, ncol = C, byrow = TRUE), dim = c(T, K+H, C), dimnames = list(seq(T), seq(K+H), seq(C)))
   
    pPrior <- array(1/(K+H), c(T, D, S, W, K+H), list(seq(T), seq(D), seq(S), seq(W), seq(K+H)))
   
    a <- metadat$AgeClass
    
    i <- metadat$z
    
    p <- pPrior
    
    i[metadat$WildGenetics]=apply(cbind(metadat[metadat$WildGenetics,c("t", "d", "s", "w")], genofreq), 1, function(tdswa){sample(K, 1, TRUE, p[tdswa[1], tdswa[2], tdswa[3], tdswa[4], seq(K)]*tdswa[-seq(4)])}) 
           
    p0 <- table(cbind(metadat[,c("t", "d", "s", "w")], i=i))
    
    pdmns <- dimnames(p0) 
        
    p[, , , , pdmns$i] <- p0 + pPrior[, , , , pdmns$i] 
    
    p <- aperm(apply(p, seq(4), rdirich), c(seq(2, 5), 1))

    i[metadat$WildNoGenetics]=apply(metadat[metadat$WildNoGenetics, c("t", "d", "s", "w")], 1, function(tdswa){sample(K, 1, TRUE, p[tdswa[1], tdswa[2], tdswa[3], tdswa[4], seq(K)])})   
    
    for(sim in seq(NSIMS)){
    
      p0 <- table(cbind(metadat[,c("t", "d", "s", "w")], i=i))
    
      pdmns <- dimnames(p0) 
    
      p <- pPrior
    
      p[, , , , pdmns$i] <- p0 + pPrior[, , , , pdmns$i] 
    
      p <- aperm(apply(p, seq(4), rdirich), c(seq(2, 5), 1))
    
      pi0 <- table(t = metadat$t, i, a)
    
      dmns <- dimnames(pi0)
    
      pi <- piPrior
    
      pi[, dmns$i, dmns$a] <- pi0 + piPrior[, dmns$i, dmns$a]
    
      pi <- aperm(apply(pi, 1:2, rdirich), c(2, 3, 1))
        
      a[is.na(metadat$AGE_EUROPEAN)] <- unlist(lapply(seq(T), function(t){sapply(i[is.na(metadat$AGE_EUROPEAN) & metadat$t == t], function(ii){sample(C, 1, TRUE, pi[t, ii, ])})}))

      a.df <- cbind(metadat[,c("t", "d", "s", "w")], a)
 
      i[metadat$WildNoGenetics]=apply(a.df[metadat$WildNoGenetics, ], 1, function(tdswa){sample(K, 1, TRUE, p[tdswa[1], tdswa[2], tdswa[3], tdswa[4], seq(K)]*pi[tdswa[1], seq(K), tdswa[5]])})   
        
      i[metadat$WildGenetics]=apply(cbind(a.df[metadat$WildGenetics, ], genofreq), 1, function(tdswa){sample(K, 1, TRUE, p[tdswa[1], tdswa[2], tdswa[3], tdswa[4], seq(K)]*pi[tdswa[1], seq(K), tdswa[5]]*tdswa[-seq(5)])}) 

      x_sum <- rowsum(x, group = sillyvecBase[i[metadat$WildGenetics]], reorder = TRUE)

      nms <- rownames(x_sum)
  
      beta_prm_prm <- beta_prm 
  
      beta_prm_prm[nms,] <- beta_prm[nms, ] + x_sum
  
      lnq <- log(t(apply(beta_prm_prm, 1, function(rw){unlist(tapply(rw, INDEX = list(rep(loci, nalleles[loci])), FUN = rdirich))})))
 
      genofreq <- exp(x%*%t(lnq))
      
      if( sim > burn & ! sim %% thin){
    
        for(t in seq(T)){

          Pi <- t(rowsum(t(pi[t, , ]), group = agevec))

          for(d in seq(D)){

            for(s in seq(S)){

              for(w in seq(W)){

                PiR <- rowsum(diag(p[t, d, s, w, ]) %*% Pi, group = groups_list[[d]])

                for(aa in seq(A)){

                   file_name <- paste0("PiR_", years[t], "Dist", Districts[d], "Sub", SubDistricts[[d]][s], "_Week", w, "_", AgeGroups[aa], chain, ".txt")

                   write.table(rbind(format(PiR[, aa], trim = TRUE, digits = 16, scientific = TRUE)), file_name, row.names = FALSE, col.names = FALSE, quote = FALSE, append = sim - thin > burn)  

                }#aa

              }#w

            }#s

          }#d

        }#t
   
      }# end if  
      
    }#sim
    
   })#chain

  end_time <- Sys.time() - beg_time

  print(end_time)

  stopCluster(cl);

  

# Harvest    #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  Harvest00 <- read.table("Harvest2011-14.txt", header=TRUE, sep="\t")

  metamerge <- sort_df(merge(merge(YearData, index.spatial), stat.week.dat), vars=c("t", "d", "s", "w"))

  Harvest <- merge( x = merge( x = nn, y = metamerge, by = c("t", "d", "s", "w")), y = Harvest00, all.x = TRUE)

  Harvest$HARVEST[is.na(Harvest$HARVEST)] <- 0

  Harvestvec <- Harvest$HARVEST

  h <- tapply(Harvest$HARVEST, Harvest[, c("t", "d", "s", "w")], c)

# Summary  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  threshold=5e-7

  
  PiR0 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(S), function(s){ lapply(seq(W), function(w){ lapply(chains, function(chain){ lapply(seq(A), function(aa){ read.table(paste0("PiR_", years[t], "Dist", Districts[d], "Sub", SubDistricts[[d]][s], "_Week", w, "_", AgeGroups[aa], chain, ".txt")) }) }) }) }) }) })
  
  R0 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(S), function(s){ lapply(seq(W), function(w){ as.mcmc.list(lapply(seq(chains), function(chain){ mcmc(Reduce("+", lapply(seq(A), function(aa){ PiR0[[t]][[d]][[s]][[w]][[chain]][[aa]] }))) })) }) }) }) })

 
####   By Subdistricts By Stat Week

  GR_R1 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(S), function(s){ lapply(seq(W), function(w){ gelman.diag(R0[[t]][[d]][[s]][[w]], autoburnin = FALSE, transform = TRUE, multivariate = FALSE)[[1]][, 1] }) }) }) })

  R1 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(S), function(s){ lapply(seq(W), function(w){ Reduce(rbind, R0[[t]][[d]][[s]][[w]]) }) }) }) })
 
  R1 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(S), function(s){ lapply(seq(W), function(w){ cbind(t(apply(R1[[t]][[d]][[s]][[w]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), P0 = sum(cl > threshold) / length(cl)) })), GR = GR_R1[[t]][[d]][[s]][[w]]) }) }) }) })

####  Subdistricts Combined By Stat Week

  R2 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(W), function(w){ as.mcmc.list(lapply(seq(chains), function(chain){ mcmc(Reduce("+", lapply(seq(S), function(s){ (n[t, d, s, w] * h[t, d, s, w] > 0) * h[t, d, s, w] * R0[[t]][[d]][[s]][[w]][[chain]] })) / ifelse(sum((n[t, d, , w] * h[t, d, , w] > 0) * h[t, d, , w]), sum((n[t, d, , w] * h[t, d, , w] > 0) * h[t, d, , w]), 1)) })) }) }) })

  GR_R2 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(W), function(w){ gelman.diag(R2[[t]][[d]][[w]], autoburnin = FALSE, transform = TRUE, multivariate = FALSE)[[1]][, 1] }) }) })

  R2 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(W), function(w){ Reduce(rbind, R2[[t]][[d]][[w]]) }) }) })

  R2 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(W), function(w){ cbind(t(apply(R2[[t]][[d]][[w]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), P0 = sum(cl > threshold) / length(cl)) })), GR = GR_R2[[t]][[d]][[w]]) }) }) })



####  Subdistricts Combined Full Season By Age


  PiR1 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(A), function(aa){ as.mcmc.list(lapply(seq(chains), function(chain){ mcmc(Reduce("+", lapply(seq(S), function(s){ Reduce("+", lapply(seq(W), function(w){(n[t, d, s, w] * h[t, d, s, w] > 0) * h[t, d, s, w] * PiR0[[t]][[d]][[s]][[w]][[chain]][[aa]] })) })) / sum((n[t, d, , ] * h[t, d, , ] > 0) * h[t, d, , ])) })) }) }) })

  GR_PiR1 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(A), function(aa){ gelman.diag(PiR1[[t]][[d]][[aa]], autoburnin = FALSE, transform = TRUE, multivariate = FALSE)[[1]][, 1]  }) }) })

  PiR1 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(A), function(aa){ Reduce(rbind, PiR1[[t]][[d]][[aa]]) }) }) })

  PiR1 <- lapply(seq(T), function(t){ lapply(seq(D), function(d){ lapply(seq(A), function(aa){ cbind(t(apply(PiR1[[t]][[d]][[aa]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), P0 = sum(cl > threshold) / length(cl)) })), GR = GR_R2[[t]][[d]][[aa]]) }) }) })

  




}




