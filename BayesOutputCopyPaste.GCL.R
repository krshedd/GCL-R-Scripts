###################################################################################################################################################################################
# This function will remove the BAYES input files (.ctl and .mix) from the server and cut/paste the BAYES output files (.BO1, .BOT, .RGN, .SUM, and .CLS) to the V drive
#
# origindir <- where your files are on the server (e.g. "C:/Kyle")
# targetdir <- where you want your files to go on the V drive (e.g. "V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt")
#   NOTE: this assumes a file structure of "BAYES/Output/" within your targetdir (e.g. ""V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/BAYES/Output/")
# sillyvec <- what your mixtures are called (e.g. "KarlukSmolt2013.1" or "SKARL13s" or dget(file="V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/Objects/KarlukMixtures.txt")[1:5])
#   NOTE: sillyvec can only be up to 5 mixtures long (1st one goes in "Bayes A", 2nd in "Bayes B", etc.)
#
# Created by Kyle Shedd on Fri Mar 27 17:10:57 2015
###################################################################################################################################################################################

BayesOutputCopyPaste.GCL <- function(origindir, targetdir, sillyvec){
  
  if(length(sillyvec) > 5) {stop(paste("sillyvec can only be up to 5 mixtures long!!!"))}
  
  # Removing input files (".ctl" and ".mix")
  controlfiles=list.files(path=origindir, pattern=".ctl", full.names=TRUE, recursive=TRUE) # All control files in origindir
  mixturefiles=list.files(path=origindir, pattern=".mix", full.names=TRUE, recursive=TRUE) # All mixture files in origindir
  
  filestoremove=sapply(sillyvec, function(silly) {c(controlfiles[grep(pattern=silly, x=controlfiles)], mixturefiles[grep(pattern=silly, x=mixturefiles)])}, simplify=FALSE) # Just the control/mixture files for sillyvec 
  
  invisible(lapply(filestoremove, file.remove)) # Remove input files from server
  
  # Copy/paste output files (".BO1", ".BOT", ".RGN", ".SUM", and ".CLS")
  outputfiles=sapply(c(".BO1", ".BOT", ".RGN", ".SUM", ".CLS"), function(outfile) {list.files(path=origindir, pattern=outfile, full.names=TRUE, recursive=TRUE)}) # All output files in origindir
  
  filestocopy=sapply(sillyvec, function(silly) {outputfiles[grep(pattern=silly, x=outputfiles)]}, simplify=FALSE) # Create list of output files by sillyvec
  
  invisible(sapply(sillyvec, function(silly) {file.copy(from=filestocopy[[silly]], to=paste(targetdir, "/BAYES/Output/", silly, sep=""))})) # Move to appropriate V drive directory
  
  invisible(lapply(filestocopy, file.remove)) # Remove output files from server
  
}