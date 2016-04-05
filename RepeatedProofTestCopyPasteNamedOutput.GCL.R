RepeatedProofTestCopyPasteNamedOutput.GCL=function(origindir,targetdir,rep,rg,chains=5,IA=FALSE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function will copy the BAYES output files (.BO1,.BOT,.RGN,.SUM) and move
  # them to wherever you want (i.e on the V drive) and then remove them from the 
  # gen server along with the BAYES input files (.bse, .ctl, .mix)
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # origindir = original directory of files
  #   ~ A character vector indicating where BAYES output files are on the server.
  # targetdir = target directory structure
  #   ~ A character vector indicating where you want BAYES files to go back on 
  #     the V drive. Assumes that each repeat has its own directory.
  # rep = the repeat you would like to move (1:10)
  #   ~ A numberic vector indicated which repeat you want to move.
  # rg = the reporting group you would like to move
  #   ~ A character vector indicating which reporting group to move (i.e. "Afognak").
  # chains = 5
  #   ~ How many chains do you have? Default is 5, but perhaps you should only have 1 chain
  #     and invest more time in looking at more repeats.
  # IA = FALSE
  #   ~ Are you performing Individual Assignment (IA)? If so, then this will also
  #     look for and move .cls files for you.
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # No output, just cuts and pastes BAYES output files from server back to 
  # and deletes the BAYES input files from the server.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Moving 5 chains of 5 repeat to 5 different instances of BAYES with each repeat in its own folder (Bayes A, Bayes B, ..., Bayes E)
  #   RGvec <- rownames(dget(file = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/Objects/MixtureProofTest.SampleSize.txt"))
  #   
  #   RepGroup <- 1 # middle.Uyak
  #   
  #   for(i in 1:5){
  #     RepeatedProofTestCopyPasteNamedOutput.GCL(
  #       origindir = paste("C:/Users/krshedd/BAYES/Bayes ", LETTERS[i], sep = ''),
  #       targetdir = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/BAYES/Mixture Proof Tests/loci89/BAYES.output",
  #       rep = i,
  #       rg = RGvec[RepGroup],
  #       chains = 5,
  #       IA = TRUE)
  #   }
  #   
  # Created by Kyle Shedd September 30, 2014 15:54:11
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  invisible(ifelse(IA==FALSE,
                   assign("output",list("BO1.BO1","BOT.BOT","RGN.RGN","SUM.SUM")),
                   assign("output",list("BO1.BO1","BOT.BOT","RGN.RGN","SUM.SUM","CLS.CLS"))))
  
  x=NULL
  for(chain in 1:chains){
    for(bayesout in output){
    x=c(x,paste(origindir,"/",rg,rep,"Chain",chain,bayesout,sep=''))
    }
  }
  
  invisible(assign("filestocopy",x))
  
  targetdir2=paste(targetdir,"/",rg,rep,sep='')
  
  invisible(file.copy(from=filestocopy,to=targetdir2,overwrite=TRUE,recursive=FALSE,copy.mode=TRUE))
  
  invisible(file.remove(filestocopy))
  
  
  filestoremove=c(paste(origindir,"/",rg,rep,".bse",sep=''),
                  paste(origindir,"/",rg,rep,".mix",sep=''),
                  paste(origindir,"/",rg,rep,"Chain",1:chains,".ctl",sep=''))

  invisible(file.remove(filestoremove))
}