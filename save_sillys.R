save_sillys <- function(sillyvec, path) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function saves sillys as R objects with `dput`, it is a wrapper for `dput`.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   sillyvec - character vector of sillys you wish to save
  #   path - character vector of where you want to save sillys
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   All sillyvec are saved as "sillyvec.txt", the .gcl is removed
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # iris <- iris  # bringing into pos = 1
  # save_sillys(sillyvec = "SGILL15D15", path = "Raw genotypes")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sillyvec.gcl <- paste0(sillyvec, ".gcl")
  if(!all(sillyvec.gcl %in% ls(, pos = 1))) {
    stop(paste0("These sillys:\n", paste(setdiff(sillyvec, ls(, pos = 1)), collapse = "\n"), "\nare not in your workspace, hoser!!!"))
  }
  files_to_save <- paste0(path, "/", sillyvec, ".txt")
  sillyvec.mat = cbind(sillyvec.gcl, files_to_save)
  invisible(apply(sillyvec.mat, 1, function(obj) {
    dput(x = get(obj[1]), file = obj[2])
  } ))
}