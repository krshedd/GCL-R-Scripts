load_objects <- function(path, pattern = NULL) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function loads R objects saved with `dput`, it is a wrapper for `dget`.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   path - character vector of where the objects you wish to load reside
  #   pattern - optinal argument so you can manually specify a pattern (i.e., specific object)
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns a character vector of objects loaded
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # setwd("V:/Analysis/1_SEAK/Sockeye/Mixture/Lynn Canal Inseason/2018/")
  # load_objects(path = "Objects", pattern = "^loci") - just loads loci objects from "Objects" dir
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  {
    if (is.null(pattern)) {
      pattern <- ".txt"
    } # defaults to .txt, if no pattern is specified.
    files_to_load <- list.files(path = path, pattern = pattern, full.names = FALSE, recursive = FALSE)
  objects <- invisible(sapply(files_to_load, function(file) {
    obj <- unlist(strsplit(x = file, split = ".txt"))
    assign(x = obj, value = dget(file = paste(path, file, sep = "/")), pos = 1)
    obj
  }, USE.NAMES = FALSE))
  print(objects)
  }
}