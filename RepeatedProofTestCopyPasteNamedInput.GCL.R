RepeatedProofTestCopyPasteNamedInput.GCL=function(origindir,targetdir,rep,rg,chains=5){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function will copy the BAYES input files (.bse, .ctl, .mix) generated 
  # by the either ProofTest.GCL or RepeatedProofTest.GCL function from the V drive
  # and paste them onto the server to run BAYES.
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # origindir = original directory of files
  #   ~ A character vector indicating where BAYES files are (i.e. BAYES.baseline, etc.)
  # targetdir = target directory structure
  #   ~ A character vector indicating where you want BAYES files to go.
  # rep = the repeat you would like to move (1:10)
  #   ~ A numberic vector indicated which repeat you want to move.
  # rg = the reporting group you would like to move
  #   ~ A character vector indicating which reporting group to move (i.e. "Afognak").
  # chains = 5
  #   ~ How many chains do you have? Default is 5, but perhaps you should only have 1 chain
  #     and invest more time in looking at more repeats.
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # No output, just copies files and pastes them on to the server for you.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Moving 5 chains of 5 repeat to 5 different instances of BAYES with each repeat in its own folder (Bayes A, Bayes B, ..., Bayes E)
  #   RGvec <- rownames(dget(file = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/Objects/MixtureProofTest.SampleSize.txt"))
  #   
  #   RepGroup <- 1 # middle.Uyak
  #   
  #   for(i in 1:5){
  #     RepeatedProofTestCopyPasteNamedInput.GCL(
  #       origindir = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/BAYES/Mixture Proof Tests/loci89",
  #       targetdir = paste("C:/Users/krshedd/BAYES/Bayes ", LETTERS[i], sep = ''),
  #       rep = i,
  #       rg = RGvec[RepGroup],
  #       chains = 5)
  #   }
  #    
  # Created by Kyle Shedd August 12, 2014 14:38:45 to be used in conjunction with RepeatedProofTest.GCL
  # Updated by Kyle Shedd September 29, 2014 13:06:54 in order to accomodate the different file structure using ProofTest.GCL in for loops with repeats named as "RG#" (i.e. Afognak1, Afognak2, etc.)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  x=NULL
  for(chain in 1:chains){x[chain]=paste(origindir,"/BAYES.control/",rg,rep,"Chain",chain,".ctl",sep='')}
  
  for(bayes in list("BAYES.baseline","BAYES.control","BAYES.mixture")){

    invisible(ifelse(bayes=="BAYES.baseline",
    assign("filestocopy",paste(origindir,"/",bayes,"/",rg,rep,".bse",sep='')),
    ifelse(bayes=="BAYES.control",assign("filestocopy",x),
           assign("filestocopy",paste(origindir,"/",bayes,"/",rg,rep,".mix",sep='')))
    ))

    invisible(file.copy(from=filestocopy,to=targetdir,overwrite=TRUE,recursive=FALSE,copy.mode=TRUE))
    
  }  
}