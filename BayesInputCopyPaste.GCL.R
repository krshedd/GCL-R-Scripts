###################################################################################################################################################################################
# This function will copy the BAYES input files (.ctl and .mix) from the V drive to the server
#
# origindir <- where your files are on the V drive (e.g. "V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt")
# targetdir <- where you want your files on the server (e.g. "C:/Kyle")
#   NOTE: this assumes a file structure of "Bayes A", "Bayes B" within your targetdir (e.g. "C:/Kyle/Bayes A")
# sillyvec <- what your mixtures are called (e.g. "KarlukSmolt2013.1" or "SKARL13s" or dget(file="V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/Objects/KarlukMixtures.txt")[1:5])
#   NOTE: sillyvec can only be up to 5 mixtures long (1st one goes in "Bayes A", 2nd in "Bayes B", etc.)
#
# Created by Kyle Shedd on Fri Mar 27 16:05:57 2015
###################################################################################################################################################################################

BayesInputCopyPaste.GCL <- function(origindir, targetdir, sillyvec){
  
  if(length(sillyvec) > 5) {stop(paste("sillyvec can only be up to 5 mixtures long!!!"))}
  
  controlfiles=list.files(path=paste(origindir, "/BAYES/Control", sep=""), pattern=".ctl", full.names=TRUE, recursive=FALSE) # All control files in origindir
  mixturefiles=list.files(path=paste(origindir, "/BAYES/Mixture", sep=""), pattern=".mix", full.names=TRUE, recursive=FALSE) # All mixture files in origindir

  filestocopy=sapply(sillyvec, function(silly) {c(controlfiles[grep(pattern=silly, x=controlfiles)], mixturefiles[grep(pattern=silly, x=mixturefiles)])}, simplify=FALSE) # Just the control/mixture files for sillyvec
  
  invisible(sapply(seq(sillyvec), function(mix) {file.copy(from=filestocopy[[mix]], to=paste(targetdir, "/Bayes ", LETTERS[mix], sep=""))}))
  
}