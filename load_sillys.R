load_sillys <- function(path, sillyvec = NULL) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function loads .gcl objects saved with `dput`, it is a wrapper for `dget`.
  # It assumes all .txt files in `path` are .gcls and will load all of them,
  # unless you specify `sillyvec`, in which case it will just load those.
  # Also assumes files do not have ".gcl" in filename.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   path - character vector of where the .gcl objects you wish to load reside
  #   sillyvec - (optional) character vector of sillys in `path`
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns a character vector of sillys loaded
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # setwd("V:/Analysis/1_SEAK/Sockeye/Mixture/Lynn Canal Inseason/2018/")
  # load_sillys(path = "Baseline genotpyes", sillyvec = c("SCKAT07E", "SCKAT07L"))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(sillyvec)) {
    files_to_load <- list.files(path = path, pattern = ".txt", full.names = FALSE, recursive = FALSE)
  } else {
    files_to_load <- paste0(sillyvec, ".txt")
    if(!all(files_to_load %in% list.files(path = path, pattern = ".txt", full.names = FALSE, recursive = FALSE))) {
      stop("Not all sillys in `sillyvec` are in `path`")
    }
  }
  objects <- invisible(sapply(files_to_load, function(file) {
    obj <- unlist(strsplit(x = file, split = ".txt"))
    silly <- paste0(obj, ".gcl")
    assign(x = silly, value = dget(file = paste(path, file, sep = "/")), pos = 1)
    obj
  }, USE.NAMES = FALSE))
  print(objects)
}