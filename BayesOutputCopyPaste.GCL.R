BayesOutputCopyPaste.GCL <- function(origindir, targetdir, sillyvec){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function will copy the BAYES output files (.BO1, .BOT, .RGN, .SUM, .CLS)
  # from the server and paste them onto the V drive. All BAYES files are removed
  # from the server (input and output). The function looks recursively into the
  # orgindir for all BAYES files that match the sillyvec.
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # origindir = original directory of files
  #   ~ A character vector indicating where BAYES output files are on the server.
  # targetdir = target directory structure
  #   ~ A character vector indicating where you want BAYES files to go (i.e. V drive).
  # sillyvec = mixtures you want to move
  #   ~ A character vector indicating which mixtures you want to move.
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # No output, just cuts and pastes BAYES output files from server back to 
  # and deletes the BAYES input files from the server.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   source("https://raw.githubusercontent.com/krshedd/GCL-R-Scripts/master/BayesOutputCopyPaste.GCL.R")
  #   
  #   BayesOutputCopyPaste.GCL(origindir="C:/Users/krshedd/BAYES", 
  #                            targetdir="V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/BAYES/Late August 89loci",
  #                            sillyvec=dget(file="V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects/LateAugustMixtures2014Strata.txt"))
  #    
  # Created by Kyle Shedd Fri Mar 27 17:10:57 2015
  # Updated by Kyle Shedd Tue Apr 05 16:46:39 2016 in order to generalize and better document
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Removing input files (".ctl" and ".mix")
  controlfiles=list.files(path=origindir, pattern=".ctl", full.names=TRUE, recursive=TRUE) # All control files in origindir
  mixturefiles=list.files(path=origindir, pattern=".mix", full.names=TRUE, recursive=TRUE) # All mixture files in origindir
  
  filestoremove=sapply(sillyvec, function(silly) {c(controlfiles[grep(pattern=silly, x=controlfiles)], mixturefiles[grep(pattern=silly, x=mixturefiles)])}, simplify=FALSE) # Just the control/mixture files for sillyvec 
  
  invisible(lapply(filestoremove, file.remove)) # Remove input files from server
  
  # Copy/paste output files (".BO1", ".BOT", ".RGN", ".SUM", and ".CLS")
  outputfiles=sapply(c(".BO1", ".BOT", ".RGN", ".SUM", ".CLS"), function(outfile) {list.files(path=origindir, pattern=outfile, full.names=TRUE, recursive=TRUE)}) # All output files in origindir
  
  filestocopy=sapply(sillyvec, function(silly) {outputfiles[grep(pattern=silly, x=outputfiles)]}, simplify=FALSE) # Create list of output files by sillyvec
  
  invisible(sapply(sillyvec, function(silly) {file.copy(from=filestocopy[[silly]], to=paste(targetdir, "/Output/", silly, sep=""))})) # Move to appropriate V drive directory
  
  invisible(lapply(filestocopy, file.remove)) # Remove output files from server
  
}