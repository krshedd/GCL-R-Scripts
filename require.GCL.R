

require.GCL <- function(packages, repos = c("http://rstudio.org/_packages", "http://cran.rstudio.com")) 

  invisible(sapply(packages, function(pkg) while(!require(pkg, character.only = TRUE)) install.packages(pkg, repos)) ) 