
if(FALSE){

# wd <- "C:/Users/jjasper/Documents"#"V:/Analysis/Staff/Jim Jasper/Sockeye/SEAK/MarkEnhancedGSI/MarkEnhancedGSIwithAge/All 2011-2014"

# nchains <- 2

# NSIMS <- 10

# burn <- NSIMS/2

# thin <- 1

# Meta Data Read-In  #################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# load("ReadIn2012and13.RData")

# sillyvecMix <- c("SGILL12D11", "SGILL12D6", "SGILL12D8", "SGILL13D11", "SGILL13D6", "SGILL13D8")

# meta.data.file <- "MetaDataAll.txt"

# agevec <- c(1, 1, 1, 1, 1, 6, 2, 3, 6, 6, 6, 6, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)

# attach("V:/Analysis/1_SEAK/Sockeye/Baseline/2011/Nov2011/BaselineUpdate_Nov2011.RData")

# sillyvecBase <- SEAK151pops

# loci2remove <- c("One_Cytb_26", "One_CO1", "One_Cytb_17", "One_c3-98", "One_GPDH", "One_MHC2_251")
 
# loci <- LocusControl$locusnames[ ! LocusControl$locusnames %in% loci2remove]

  setwd(wd)

  while(!require(reshape)){install.packages("reshape")}

  while(!require(abind)){install.packages("abind")}

  while(!require(foreach)){ install.packages("foreach") }

  while(!require(doParallel)){ install.packages("doParallel") }

  metadat0 <- read.table(meta.data.file, header=TRUE, sep="\t")

  names(loci) <- loci

  nloci <- length(loci)

  nalleles <- LocusControl$nalleles[loci]

  invisible(PoolCollections.GCL(collections = sillyvecMix, loci = loci))

  x0 <- get(paste(c(sillyvecMix, "gcl"), collapse="."))$counts[,  loci, ]

  rownames(x0) <- get(paste(c(sillyvecMix, "gcl"), collapse="."))$attributes$SillySource
  
# Meta Data Massaging  ########################################################################################################################################################################################################################################################################################################################################################################################################################################################

  years  <- sort(unique(metadat0$YEAR))

  T  <- length(years)

  YearData  <- data.frame(YEAR=years, t=seq(T))

  EuropeanAge <- sapply(metadat0$AGE_EUROPEAN, function(age){paste0(rep(0, 2-nchar(age)), age)})

  EuropeanAge[is.na(metadat0$AGE_EUROPEAN)] <- NA

  EuropeanAges <- sort(apply(expand.grid(seq(0, 5),seq(1, 5)), 1 ,paste, collapse = ""))
  
  AgeNames <- paste0("Age", EuropeanAges)

  agevec <- setNames(agevec, AgeNames)

  AgeGroups <- c("ZeroX", "Age1.2", "Age1.3", "Age2.2", "Age2.3", "AgeOther")

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

  Hatcheries <- levels(metadat2$ThermalGroupsForSEAK)

  Hatcheries <- Hatcheries[Hatcheries != "WILD"]

  H <- length(Hatcheries)

  Originmetadat <- data.frame(ThermalGroupsForSEAK = levels(metadat2$ThermalGroupsForSEAK), z = c(seq(H) + K, NA))

  metadat3 <- merge(metadat2, Originmetadat)

  metadat <- sort_df(metadat3, vars=c("t", "NonWild", "WildNoGenetics", "WildGenetics", "d", "s", "w"))

  n <- tapply(metadat$WildGenetics, metadat[, c("t", "d", "s", "w")], length)

  n[is.na(n)] <- 0

  n0 <- as.data.frame.table(tapply(metadat$WildGenetics, metadat[, c("t", "d", "s", "w")], sum), responseName="WildGenetics")

  WildNoGenetics0 <- as.data.frame.table(tapply(metadat$WildNoGenetics, metadat[, c("t", "d", "s", "w")], sum), responseName="WildNoGenetics")$WildNoGenetics

  NonWild0 <- as.data.frame.table(tapply(metadat$NonWild, metadat[, c("t", "d", "s", "w")], sum), responseName="NonWild")$NonWild

  nn <- data.frame(n0, WildNoGenetics=WildNoGenetics0, NonWild=NonWild0)

# head(metadat)

  rownames(metadat) <- metadat$SillySource

  sillyvecBase <- c(sillyvecBase, Hatcheries)


# Analysis  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

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
   
    lnq <- log(t(apply(beta_prm, 1, function(rw){Reduce(cbind, tapply(rw, list(rep(loci, nalleles)), rdirich, simplify = FALSE))})))
   
    genofreq <- exp(x%*%t(lnq))
   
    piPrior <- array(1/C, c(T, K+H, C), list(seq(T), seq(K+H), seq(C)))
   
    pPrior <- array(1/(K+H), c(T, D, S, W, K+H), list(seq(T), seq(D), seq(S), seq(W), seq(K+H)))
   
    a <- metadat$AgeClass
    
    i <- metadat$z
    
    i[metadat$WildGenetics] <- apply(genofreq, 1, function(freqvec){sample(K, 1, TRUE, freqvec)})    
    
    pi0 <- table(t = metadat$t, i, a)
    
    dmns <- sapply(dimnames(pi0), as.numeric)
    
    pi <- piPrior
    
    pi[, dmns$i, dmns$a] <- pi0 + piPrior[, dmns$i, dmns$a]
    
    pi <- aperm(apply(pi, 1:2, rdirich), c(2, 3, 1))
    
    p0 <- table(cbind(metadat[,c("t", "d", "s", "w")], i=i))
    
    pdmns <- sapply(dimnames(p0), as.numeric) 
    
    p <- pPrior
    
    p[, , , , pdmns$i] <- p0 + pPrior[, , , , pdmns$i]
    
    p <- aperm(apply(p, seq(4), rdirich), c(seq(2, 5), 1))
    
    i[metadat$WildNoGenetics]=apply(metadat[metadat$WildNoGenetics,c("t", "d", "s", "w")], 1, function(rwvec){sample(K, 1, TRUE, p[rwvec[1], rwvec[2], rwvec[3], rwvec[4], seq(K)])})
    
    a[is.na(metadat$AGE_EUROPEAN)] <- unlist(lapply(seq(T), function(t){sapply(i[is.na(metadat$AGE_EUROPEAN) & metadat$t == t], function(Z){sample(C, 1, TRUE, pi[t, Z, ])})}))
    
    for(sim in seq(NSIMS)){
    
      pi0 <- table(t = metadat$t, i, a))
    
      dmns <- sapply(dimnames(pi0), as.numeric)
    
      pi <- piPrior
    
      pi[, dmns$i, dmns$a] <- pi0 + piPrior[, dmns$i, dmns$a]
    
      pi <- aperm(apply(pi, 1:2, rdirich), c(2, 3, 1))
    
      p0 <- table(cbind(metadat[,c("t", "d", "s", "w")], i=i))
    
      pdmns <- sapply(dimnames(p0), as.numeric) 
    
      p <- pPrior
    
      p[, , , , pdmns$i] <- p0 + pPrior[, , , , pdmns$i] 
    
      p <- aperm(apply(p, seq(4), rdirich), c(seq(2, 5), 1))
    
      a[is.na(metadat$AGE_EUROPEAN)] <- unlist(lapply(seq(T), function(t){sapply(i[is.na(metadat$AGE_EUROPEAN) & metadat$t == t], function(ii){sample(C, 1, TRUE, pi[t, ii, ])})}))
    
      a.df <- cbind(metadat[,c("t", "d", "s", "w")], a)
    
      i[metadat$WildNoGenetics]=apply(a.df[metadat$WildNoGenetics, ], 1, function(tdswa){sample(K, 1, TRUE, p[tdswa[1], tdswa[2], tdswa[3], tdswa[4], seq(K)]*pi[tdswa[1], seq(K), tdswa[5]])})   
    
      i[metadat$WildGenetics]=apply(cbind(a.df[metadat$WildGenetics, ], genofreq), 1, function(tdswa){sample(K, 1, TRUE, p[tdswa[1], tdswa[2], tdswa[3], tdswa[4], seq(K)]*pi[tdswa[1], seq(K), tdswa[5]]*tdswa[-seq(5)])}) 
    
      if( sim > burn & ! sim %% thin){
    
        invisible(lapply(seq(T), function(t){lapply(seq(D), function(d){lapply(seq(S), function(s){lapply(seq(W), function(w){ cat(format(p[t, d, s, w, ], trim = TRUE, digits = 16, scientific = TRUE), "\n", file = paste0("P_", years[t], "Dist", Districts[d], "Sub", SubDistricts[[d]][s], "_Week", w, chain, ".txt")) }) }) }) }) )
    
        invisible(lapply(seq(T), function(t){lapply(seq(C), function(c){ cat(format(pi[t, , c], trim = TRUE, digits = 16, scientific = TRUE), "\n", file = paste0("Pi_", years[t], AgeNames[c], chain, ".txt")) }) }) )
    
      }# end if  
      
    }#sim
    
   })#chain

  end.time <- Sys.time() - beg.time


  stopCluster(cl);

  

# Harvest    #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  Harvest00 <- read.table("Harvest2011-14.txt", header=TRUE, sep="\t")

  metamerge <- sort_df(merge(merge(YearData, index.spatial), stat.week.dat), vars=c("t", "d", "s", "w"))

  Harvest <- merge( x = merge( x = nn, y = metamerge, by = c("t", "d", "s", "w")), y = Harvest00, all.x = TRUE)

  Harvest$HARVEST[is.na(Harvest$HARVEST)] <- 0

  Harvestvec <- Harvest$HARVEST

# Summary  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  threshold=5e-7

  groups=list(dget(paste0(dir, "/GroupsD106.txt")), dget(paste0(dir, "/GroupsD108.txt")), dget(paste0(dir, "/GroupsD111.txt")))

  groupvecNames=dget(paste0(dir, "/groupvecNames.txt"))

  groupvecs=list(  setNames(dget(paste0(dir, "/GroupvecD106.txt")), groupvecNames)[sillyvecBase], 
                   setNames(dget(paste0(dir, "/GroupvecD108.txt")), groupvecNames)[sillyvecBase], 
                   setNames(dget(paste0(dir, "/GroupvecD111.txt")), groupvecNames)[sillyvecBase]
                 )





  Pi <- lapply(seq(T), function(t){lapply(seq(C), function(c){ read.table(paste0("Pi_", years[t], AgeNames[c], chain, ".txt")) }) }) 

  P <- lapply(seq(T), function(t){lapply(seq(D), function(d){lapply(seq(S), function(s){lapply(seq(W), function(w){ read.table(paste0("P_", years[t], "Dist", Districts[d], "Sub", SubDistricts[[d]][s], "_Week", w, chain, ".txt")) }) }) }) })



}




