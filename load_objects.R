load_objects <- function(path) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function loads R objects saved with `dput`, it is a wrapper for `dget`.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   path - character vector of where the objects you wish to load reside
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns a character vector of objects loaded
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # setwd("V:/Analysis/1_SEAK/Sockeye/Mixture/Lynn Canal Inseason/2018/")
  # load_objects(path = "Objects")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  files_to_load <- list.files(path = path, pattern = ".txt", full.names = FALSE, recursive = FALSE)
  objects <- invisible(sapply(files_to_load, function(file) {
    obj <- unlist(strsplit(x = file, split = ".txt"))
    assign(x = obj, value = dget(file = paste(path, file, sep = "/")), pos = 1)
    obj
  }, USE.NAMES = FALSE))
  print(objects)
}