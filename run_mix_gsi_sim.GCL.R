run_mix_gsi_sim.GCL <- function(mix.path, base.path, out.dir, groups, group.names, group.prior, n.its, n.burn = floor(n.its / 2)){

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# load("Z:\\Analysis\\2_Central\\Chinook\\Cook Inlet\\2014\\Mixture\\2015 ESSN Chinook postseason analysis\\ESSN2015PostseasonAnalysis.RData")
#
# attach("V:/Analysis/2_Central/Chinook/Cook Inlet/2014/Baseline/CIChinook2014Baseline.RData")
#
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# out.dir <- "C:/Users/jjasper/Documents/R/gsi_sim/test_data"
#
# mix.path <- "C:/Users/jjasper/Documents/R/gsi_sim/test_data/KasilofEarly2015.txt"
#
# base.path <- "C:/Users/jjasper/Documents/R/gsi_sim/test_data/CIbaselineChinook2014.txt"
#
# groups <- groupvec
#
# group.names <- c("KenaiTrib","KenaiMain","KasilofMain","CIOther") 
#
# group.pi.prior <- KasilofEarly_prior
#
# n.its <- 25000
#
# n.burn <- 5000
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  setwd(out.dir)

  base.collections <- readLines(base.path)

  base.collections <- base.collections[substr(base.collections, 1, 3) == "POP"]

  base.collections <- substr(base.collections, 5, nchar(base.collections))

  out.path <- paste0(getwd(), "/output_", basename(mix.path))

  groups.path <- paste0(getwd(), "/groups_", basename(base.path))

  groups.file <- tapply(base.collections, list(group.names[groups]), c, simplify=FALSE)

  groups.file <- Reduce(rbind, lapply(group.names, function(group){cbind(c(" ", paste("REPUNIT", group), groups.file[[group]]))}))

  write.table(groups.file, groups.path, row.names = FALSE, col.names = FALSE, quote = FALSE)

  command <- paste0(path.expand("~\\R"), "\\gsi_sim\\gsi_sim-MINGW64_NT-6.1")

  args <- c(paste("-b", base.path), paste("--rep-unit-file", groups.path), paste("-t", mix.path), paste("--mcmc-sweeps", n.its), paste("--mcmc-burnin", n.burn), paste(c("--bayes-pi-prior", Prior.GCL(groups, group.prior)), collapse = " "))

  err <- system2(command = command, args = args, stdout = out.path, stderr = TRUE) 

  return(err)

} 


