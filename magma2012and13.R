if(FALSE){### This line prevents any commands in this script from being accidently executed upon being sourced. 


##### Mixture read-in and Prep ####################################

  source("H:/Desktop/R/Functions.GCL.R")

  username <- ""

  password <- ""

  sillyvecMix <- c("SGILL16D11","SGILL16D6","SGILL16D8")

  CreateLocusControl.GCL("Sockeye2011_96SNPs",username,password)

  read_time <- LOKI2R.GCL(sillyvecMix, username, password) 

  rm(username, password)

#### Data Massage
  ColSize <- sapply(paste(sillyvecMix,".gcl",sep=''), function(x) get(x)$n)
  
  RemoveMissing <- RemoveIndMissLoci.GCL(sillyvecMix, proportion=0.8)

  MyDupCheck <- CheckDupWithinSilly.GCL(sillyvecMix, loci = LocusControl$locusnames, quantile = NULL)

  RemoveDups <- RemoveDups.GCL(MyDupCheck)

  invisible(CombineLoci.GCL(sillyvecMix, markerset = c("One_CO1", "One_Cytb_17", "One_Cytb_26"), update = TRUE))

  loci2remove <- c("One_Cytb_26", "One_CO1", "One_Cytb_17", "One_c3-98", "One_GPDH", "One_MHC2_251")
 
  loci <- LocusControl$locusnames[ ! LocusControl$locusnames %in% loci2remove]

  names(loci) <- loci

  nloci <- length(loci)

  nalleles <- LocusControl$nalleles[loci]

  invisible(PoolCollections.GCL(sillyvecMix, loci = loci))

  all.gcl <- get(paste(c(sillyvecMix, "gcl"), collapse = "."))

 attributes0 <- all.gcl$attributes

  x0 <- all.gcl$counts[,loci,]

  rownames(x0) <- attributes0$SillySource
  
  require(xlsx)
  
  write.xlsx(c(RemoveDups[1],RemoveDups[2],RemoveDups[3]),file="C:\\Analysis\\MAGMA\\TBR 2016\\DupsRemoved.xlsx")
  
  FinalColSize <- sapply(paste(sillyvecMix,".gcl",sep=''), function(x) get(x)$n)

############## magma2012and2013 function starts here #############################################################################################################
  
  # ARGS  ########################################################################################################################################################################################################################################################################################################################################################################################################################################################
  
  wd <- "C:\\Analysis\\MAGMA\\TBR 2016"
  
  nchains <- 6
  
  NSIMS <- 80000
  
  burn <- NSIMS / 2
  
  thin <- 10
    
  meta.data.file <- "TBR_metadata_2016.txt"  

  agevec <- c(1, 1, 1, 1, 1, 6, 2, 3, 6, 6, 6, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)
  
  AgeGroups <- c("ZeroX", "Age1.2", "Age1.3", "Age2.2", "Age2.3", "AgeOther")
  
  attach("V:\\Analysis\\1_SEAK\\Sockeye\\Baseline\\2013_2014\\SRO\\SEAKbaseline2014.RData")

  sillyvecBase <- OrderedPops

  y0 <- FreqPop.GCL(sillyvecBase, loci)

  detach("file:V:\\Analysis\\1_SEAK\\Sockeye\\Baseline\\2013_2014\\SRO\\SEAKbaseline2014.RData", character.only = TRUE)
   
  Hatcheries <- c("BURN", "MAIN", "SPEE", "SWEE", "TAHL", "TATS", "TUYA", "KING","TRAP")

  setwd(wd)

  groups_files <- c("groupvecD6.txt", "groupvecD8.txt", "groupvecD11.txt")

  group_name_files <- c("GroupsD6.txt", "GroupsD8.txt", "GroupsD11.txt")
 
  groups_list <- lapply(lapply(groups_files, dget), setNames, nm = c(sillyvecBase, Hatcheries))
  write.xlsx(groups_list,file="C:\\Analysis\\MAGMA\\TBR 2016\\GroupvecCheck.xlsx")

  group_names_list <- list(lapply(group_name_files, dget))

##################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  while(!require(reshape)){install.packages("reshape")}

  while(!require(coda)){install.packages("coda")}

  while(!require(foreach)){ install.packages("foreach") }

  while(!require(doParallel)){ install.packages("doParallel") }

  metadat0 <- subset(read.table(meta.data.file, header=TRUE, sep="\t", stringsAsFactors = FALSE),SILLY %in% sillyvecMix)

## Meta Data Massaging  ########################################################################################################################################################################################################################################################################################################################################################################################################################################################

  years  <- sort(unique(metadat0$YEAR))

  T  <- length(years)

  YearData  <- data.frame(YEAR = years, t = seq(T))
  
  metadat0$AGE_EUROPEAN[is.na(metadat0$AGE_EUROPEAN)] = "NA"
  
  EuropeanAge <- sapply(metadat0$AGE_EUROPEAN, function(age){paste0(rep(0, 2-nchar(age)), age)})
  
  metadat0$AGE_EUROPEAN[   metadat0$AGE_EUROPEAN == "NA" ] = NA

  #EuropeanAge[is.na(metadat0$AGE_EUROPEAN)] <- NA

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

  SubDistricts <- sapply(Districts, function(d){as.vector(sort(unique(metadat1$SUBDISTRICT[metadat1$DISTRICT%in%d])))}, simplify=FALSE)

  S  <- unique(sapply(SubDistricts, length))

  index.spatial <- NULL 

  for(d in Districts){
  
    district <- match(d, Districts)
  
    index.spatial <- rbind(index.spatial, data.frame(DISTRICT=d, SUBDISTRICT=SubDistricts[[district]], d=district, s=seq(length(SubDistricts[[district]]))))
  
  }#d

  StatWeeks=sort(unique(metadat1$STAT_WEEK))

  W <- length(StatWeeks)

  stat.week.dat <- data.frame(STAT_WEEK=StatWeeks, w=seq(W))

  metadat2 <- merge(x = merge(x = merge(x = metadat1, y = YearData), y =  index.spatial), y = stat.week.dat)

# Baseline Wild  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  y <- Reduce(cbind, lapply(loci, function(locus){y0[, locus, seq(nalleles[locus])]}))

  colnames(y) <- unlist(lapply(loci, function(locus){ as.vector(paste(locus, seq(nalleles[locus]), sep = "_")) }))

  K <- nrow(y)

# Baseline Hatcheries  ################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  H <- length(Hatcheries)

  Originmetadat <- data.frame(ThermalGroupsForSEAK = c(Hatcheries, "WILD"), i = c(seq(H) + K, NA))

  metadat3 <- merge(metadat2, Originmetadat)

  metadat <- sort_df(metadat3, vars=c("t", "NonWild", "WildNoGenetics", "WildGenetics", "d", "s", "w"))

  n <- tapply(metadat$WildGenetics, metadat[, c("t", "d", "s", "w")], length)

  n[is.na(n)] <- 0

  n0 <- as.data.frame.table(tapply(metadat$WildGenetics, metadat[, c("t", "d", "s", "w")], sum), responseName="WildGenetics")

  WildNoGenetics0 <- as.data.frame.table(tapply(metadat$WildNoGenetics, metadat[, c("t", "d", "s", "w")], sum), responseName="WildNoGenetics")$WildNoGenetics

  NonWild0 <- as.data.frame.table(tapply(metadat$NonWild, metadat[, c("t", "d", "s", "w")], sum), responseName="NonWild")$NonWild

  nn <- data.frame(n0, WildNoGenetics=WildNoGenetics0, NonWild=NonWild0)

  ## Dumping nn to double-check samples sizes
  write.xlsx(nn,file="C:\\Analysis\\MAGMA\\TBR 2016\\SampleSizeCheckAfterD108Update_nn.xlsx")
  
  
  rownames(metadat) <- metadat$SillySource

  sillyvecBase <- c(sillyvecBase, Hatcheries)

# Analysis  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  if(any(! metadat$SillySource[metadat$WildGenetics] %in% rownames(x0))){ print("meta-data is messed up.") }

  x <- Reduce(cbind, lapply(loci, function(locus){ x0[metadat$SillySource[metadat$WildGenetics], locus, seq(nalleles[locus])] }) )

  x[is.na(x)] <- 0  

  colnames(x) <- colnames(y)

  chains <- paste0("Chain", seq(nchains))

  cl <- makePSOCKcluster(nchains)

  registerDoParallel(cl, cores = nchains)

  update_base <- TRUE # Argument for turning on/off baseline updating--useful for troubleshooting.

  nadapt <- 100 # A tuning parameter used to set the number of iterations to run before baseline updating begins. 

  beg_time <- Sys.time()

  invisible(foreach(chain = chains)%dopar%{
  
    rdirich <- function(alpha0){ vec <- rgamma(length(alpha0), alpha0, 1) ; vec / sum(vec) }
  
    beta <- matrix(rep(1/nalleles, nalleles), nrow = nrow(y), ncol = ncol(y), byrow = TRUE, dimnames = dimnames(y))
  
    beta_prm <- y + beta

    loc_all_fac <- factor(rep(loci, nalleles), levels = loci)
  
    t_lnq <- log(apply(beta_prm, 1, function(rw){unlist(tapply(rw, loc_all_fac, function(alls){ alls / sum(alls) }, simplify = FALSE))}))
  
    genofreq <- exp(x %*% t_lnq)
  
    piPrior <- array(matrix((1 / A / table(agevec))[agevec], nrow = K + H, ncol = C, byrow = TRUE), dim = c(T, K+H, C), dimnames = list(seq(T), seq(K+H), seq(C)))
  
    pPrior <- array(1/(K+H), c(T, D, S, W, K+H), list(seq(T), seq(D), seq(S), seq(W), seq(K+H)))
  
    na_a <- is.na(metadat$AgeClass)

    tdsw <- metadat[,c("t", "d", "s", "w")]

    wild_gen <- which(metadat$WildGenetics)

    wild_no_gen <- which(metadat$WildNoGenetics)

    a <- metadat$AgeClass
  
    i <- metadat$i 
  
    p <- pPrior 
  
    i[wild_gen] <- apply(cbind(tdsw[wild_gen, ], genofreq), 1, function(tdswa){sample(K, 1, TRUE, p[tdswa[1], tdswa[2], tdswa[3], tdswa[4], seq(K)]*tdswa[-seq(4)])})

    p0 <- table(cbind(tdsw, i = i))
    
    pdmns <- dimnames(p0) 
    
    p[pdmns$t,,,, pdmns$i] <- p0[pdmns$t,,,, pdmns$i] + pPrior[pdmns$t,,,, pdmns$i]

    i[wild_no_gen] <- apply(tdsw[wild_no_gen, ], 1, function(tdswa){sample(K, 1, TRUE, p[tdswa[1], tdswa[2], tdswa[3], tdswa[4], seq(K)])})   
  
    for(sim in seq(NSIMS)){
    
      p0 <- table(cbind(tdsw, i = i))
    
      pdmns <- dimnames(p0) 
    
      p <- pPrior
    
      p[pdmns$t,,,, pdmns$i] <- p0[pdmns$t,,,, pdmns$i] + pPrior[pdmns$t,,,, pdmns$i]
    
      p <- aperm(apply(p, seq(4), rdirich), c(seq(2, 5), 1))
    
      pi0 <- table(t = tdsw$t, i, a)
    
      dmns <- dimnames(pi0)
    
      pi <- piPrior
 
      pi[dmns$t, dmns$i, dmns$a] <- pi0[dmns$t, dmns$i, dmns$a] + piPrior[dmns$t, dmns$i, dmns$a]   
    
      pi <- aperm(apply(pi, 1:2, rdirich), c(2, 3, 1))
    
      a[na_a] <- unlist(lapply(seq(T), function(t){sapply(i[na_a & tdsw$t == t], function(ii){sample(C, 1, TRUE, pi[t, ii, ])})}))
    
      i[wild_gen] <- apply(cbind(cbind(tdsw[wild_gen, ], a = a[wild_gen]), genofreq), 1, function(tdswa){sample(K, 1, TRUE, p[tdswa[1], tdswa[2], tdswa[3], tdswa[4], seq(K)]*pi[tdswa[1], seq(K), tdswa[5]]*tdswa[-seq(5)])}) 
    
      i[wild_no_gen] <- apply(cbind(tdsw[wild_no_gen, ], a = a[wild_no_gen]), 1, function(tdswa){sample(K, 1, TRUE, p[tdswa[1], tdswa[2], tdswa[3], tdswa[4], seq(K)]*pi[tdswa[1], seq(K), tdswa[5]])})   
    
      if(update_base & sim > nadapt){

        x_sum <- rowsum(x, group = sillyvecBase[i[wild_gen]], reorder = TRUE)
    
        nms <- rownames(x_sum)
    
        beta_prm_prm <- beta_prm 
    
        beta_prm_prm[nms,] <- beta_prm[nms, ] + x_sum
    
        t_lnq <- log(apply(beta_prm_prm, 1, function(rw){ unlist(tapply(rw, INDEX = loc_all_fac, FUN = rdirich)) }) )
    
        genofreq <- exp(x %*% t_lnq)

      }#update_base 
    
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

  stopCluster(cl)

save.image ("C:\\Analysis\\MAGMA\\TBR 2016\\TBR2016.RData")
#load("C:\\Analysis\\MAGMA\\TBR 2016\\TBR2016.RData")
# Harvest ############################################################################################################################################################################################################################

  Harvest00 <- subset(read.table("2016_TBR_Harvest_forR.txt", header=TRUE, sep="\t"), YEAR %in% years)

  metamerge <- sort_df(merge(merge(YearData, index.spatial), stat.week.dat), vars=c("t", "d", "s", "w"))

  Harvest <- merge( x = merge( x = nn, y = metamerge, by = c("t", "d", "s", "w")), y = Harvest00, all.x = TRUE)

  Harvest$HARVEST[is.na(Harvest$HARVEST)] <- 0
  
  write.xlsx(Harvest,file="C:\\Analysis\\MAGMA\\TBR 2016\\HarvestCheck.xlsx")

  h <- tapply(Harvest$HARVEST, Harvest[, c("t", "d", "s", "w")], c)
  
  h[is.na(h)] <- 0

# useful things ###################################################################################################################################################################################################################

  TT <- setNames(seq(T), years)

  AA <- setNames(seq(A), AgeGroups)

  DD <- setNames(seq(D), Districts)

  SS <- sapply(SubDistricts, function(nms){ setNames(seq(nms), nms) }, simplify = FALSE)

  WW <- setNames(seq(W), paste0("StatWeek", StatWeeks))

  CHAINS <- setNames(seq(chains), chains)

  threshold <- 5e-7

# Summary ################################################################################################################################################################################################################################

  PiR0 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ lapply(WW, function(w){ lapply(CHAINS, function(chain){ lapply(AA, function(aa){ read.table(paste0("PiR_", years[t], "Dist", Districts[d], "Sub", SubDistricts[[d]][s], "_Week", w, "_", AgeGroups[aa], "Chain",chain, ".txt")) }) }) }) }) }) })

  R0 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ lapply(WW, function(w){ as.mcmc.list(lapply(CHAINS, function(chain){ mcmc(Reduce("+", lapply(seq(A), function(aa){ PiR0[[t]][[d]][[s]][[w]][[chain]][[aa]] }))) })) }) }) }) })

####   By Subdistricts By Stat Week

  GR_R1 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ lapply(WW, function(w){ gelman.diag(R0[[t]][[d]][[s]][[w]], autoburnin = FALSE, transform = TRUE, multivariate = FALSE)[[1]][, 1] }) }) }) })

  R1 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ lapply(WW, function(w){ Reduce(rbind, R0[[t]][[d]][[s]][[w]]) }) }) }) })

  R1 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ lapply(WW, function(w){  mat <- cbind(t(apply((n[t, d, s, w] > 0) * R1[[t]][[d]][[s]][[w]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), P0 = sum(cl < threshold) / length(cl)) })), GR = GR_R1[[t]][[d]][[s]][[w]]); rownames(mat) <- group_names_list[[1]][[d]] ; mat }) }) }) })


####  Subdistricts Combined By Stat Week

  R2 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(WW, function(w){ as.mcmc.list(lapply(CHAINS, function(chain){ mcmc(Reduce("+", lapply(SS[[d]], function(s){ (n[t, d, s, w] > 0) * h[t, d, s, w] * R0[[t]][[d]][[s]][[w]][[chain]] })) / ifelse(sum((n[t, d, , w] * h[t, d, , w] > 0) * h[t, d, , w]), sum((n[t, d, , w] * h[t, d, , w] > 0) * h[t, d, , w]), 1)) })) }) }) })

  GR_R2 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(WW, function(w){ gelman.diag(R2[[t]][[d]][[w]], autoburnin = FALSE, transform = TRUE, multivariate = FALSE)[[1]][, 1] }) }) })

  R2 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(WW, function(w){ Reduce(rbind, R2[[t]][[d]][[w]]) }) }) })

  R2 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(WW, function(w){ mat <- cbind(t(apply(R2[[t]][[d]][[w]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), P0 = sum(cl < threshold) / length(cl)) })), GR = GR_R2[[t]][[d]][[w]]) ; rownames(mat) <- group_names_list[[1]][[d]] ; mat  }) }) })


  
####  Subdistricts Combined  By Age

  PiR1 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(AA, function(aa){ as.mcmc.list(lapply(CHAINS, function(chain){ mcmc(Reduce("+", lapply(SS[[d]], function(s){ Reduce("+", lapply(WW, function(w){(n[t, d, s, w] > 0) * h[t, d, s, w] * PiR0[[t]][[d]][[s]][[w]][[chain]][[aa]] })) })) / sum((n[t, d, , ] > 0) * h[t, d, , ])) })) }) }) })

  GR_PiR1 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(AA, function(aa){ gelman.diag(PiR1[[t]][[d]][[aa]], autoburnin = FALSE, transform = TRUE, multivariate = FALSE)[[1]][, 1]  }) }) })

  PiR1 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(AA, function(aa){ Reduce(rbind, PiR1[[t]][[d]][[aa]]) }) }) })

  PiR1 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(AA, function(aa){ mat <- cbind(t(apply(PiR1[[t]][[d]][[aa]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), P0 = sum(cl < threshold) / length(cl)) })), GR = GR_PiR1[[t]][[d]][[aa]]) ; rownames(mat) <- group_names_list[[1]][[d]] ; mat }) }) })



  ### Stat Weeks Combined By SubDistrict
  
  R3 <-  lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ as.mcmc.list(lapply(CHAINS, function(chain){ mcmc(Reduce("+", lapply(WW, function(w){ (n[t, d, s, w] > 0) * h[t, d, s, w] * R0[[t]][[d]][[s]][[w]][[chain]] })) / sum((n[t, d, s, ] > 0) * h[t, d, s, ])) })) }) }) })
  
  GR_R3 <-  lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ gelman.diag(R3[[t]][[d]][[s]], autoburnin = FALSE, transform = TRUE, multivariate = FALSE)[[1]][, 1] }) }) })
  
  R3 <-  lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ Reduce(rbind, R3[[t]][[d]][[s]]) }) }) })
  
  R3 <-  lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ mat <- cbind(t(apply(R3[[t]][[d]][[s]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), P0 = sum(cl < threshold) / length(cl)) })), GR = GR_R3[[t]][[d]][[s]]) ; rownames(mat) <- group_names_list[[1]][[d]] ; mat }) }) })
  
## By subdistrict combined
  R4 <-  lapply(TT, function(t){ lapply(DD, function(d){ mcmc.list(lapply(CHAINS, function(chain){ mcmc(Reduce("+", lapply(SS[[d]], function(s){ Reduce("+", lapply(WW, function(w){ (n[t, d, s, w] > 0) * h[t, d, s, w] * R0[[t]][[d]][[s]][[w]][[chain]] })) })) / sum((n[t, d, , ] > 0) * h[t, d, , ])) })) }) })
  
  GR_R4 <-  lapply(TT, function(t){ lapply(DD, function(d){ gelman.diag(R4[[t]][[d]], autoburnin = FALSE, transform = TRUE, multivariate = FALSE)[[1]][, 1] }) })
  
  R4 <-  lapply(TT, function(t){ lapply(DD, function(d){ Reduce(rbind, R4[[t]][[d]]) }) })
  
  R4 <-  lapply(TT, function(t){ lapply(DD, function(d){ mat <- cbind(t(apply(R4[[t]][[d]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), P0 = sum(cl < threshold) / length(cl)) })), GR = GR_R4[[t]][[d]]) ; rownames(mat) <- group_names_list[[1]][[d]] ; mat }) })
  
  
  ####  Full Season By Subdistricts By Age 
  
  PiR2 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ lapply(AA, function(aa){ as.mcmc.list(lapply(CHAINS, function(chain){ mcmc(Reduce("+", lapply(WW, function(w){(n[t, d, s, w] > 0) * h[t, d, s, w] * PiR0[[t]][[d]][[s]][[w]][[chain]][[aa]] })) / sum((n[t, d, s, ] > 0) * h[t, d, s, ])) })) }) }) }) })
  
  GR_PiR2 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ lapply(AA, function(aa){ gelman.diag(PiR2[[t]][[d]][[s]][[aa]], autoburnin = FALSE, transform = TRUE, multivariate = FALSE)[[1]][, 1]  }) }) }) })
  
  PiR2 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ lapply(AA, function(aa){ Reduce(rbind, PiR2[[t]][[d]][[s]][[aa]]) }) }) }) })
  
  PiR2 <- lapply(TT, function(t){ lapply(DD, function(d){ lapply(SS[[d]], function(s){ lapply(AA, function(aa){ mat <- cbind(t(apply(PiR2[[t]][[d]][[s]][[aa]], 2, function(cl){ c(mean = mean(cl), sd = sd(cl), median = quantile(cl, 0.5), ci = quantile(cl, 0.05), ci = quantile(cl, 0.95), P0 = sum(cl < threshold) / length(cl)) })), GR = GR_PiR2[[t]][[d]][[s]][[aa]]) ; rownames(mat) <- group_names_list[[1]][[d]] ; mat }) }) }) })
  
}#end if(FALSE)



# Saving results to file

sink(file="C:\\Analysis\\MAGMA\\TBR 2016\\MAGMA output\\BySDbySW_R1.txt")

print(R1)

sink()

R1.df <- Reduce( R1)

sink(file="C:\\Analysis\\MAGMA\\TBR 2016\\MAGMA output\\DistrictbySW_R2.txt")

print(R2)

sink()


sink(file="C:\\Analysis\\MAGMA\\TBR 2016\\MAGMA output\\FUllSeasonbySD_R3.txt")

print(R3)

sink()


sink(file="C:\\Analysis\\MAGMA\\TBR 2016\\MAGMA output\\FUllSeasonDistricts_R4.txt")

print(R4)

sink()


sink(file="C:\\Analysis\\MAGMA\\TBR 2016\\MAGMA output\\DistrictsbyAge_PiR1.txt")

print(PiR1)

sink()


sink(file="C:\\Analysis\\MAGMA\\TBR 2016\\MAGMA output\\DistrictsbyAgeSD_PiR2.txt")

print(PiR2)

sink()
