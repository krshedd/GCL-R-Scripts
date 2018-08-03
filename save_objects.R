save_objects <- function(objects, path) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function saves R objects with `dput`, it is a wrapper for `dput`.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   objects - character vector of objects you wish to save
  #   path - character vector of where you want to save objects
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   All objects are saved as "objects.txt"
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # iris <- iris  # bringing into pos = 1
  # save_objects(objects = "iris", path = "Objects")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(!all(objects %in% ls(, pos = 1))) {
    stop(paste0("These objects:\n", paste(setdiff(objects, ls(, pos = 1)), collapse = "\n"), "\nare not in your workspace, hoser!!!"))
  }
  files_to_save <- paste0(path, "/", objects, ".txt")
  objects.mat = cbind(objects, files_to_save)
  invisible(apply(objects.mat, 1, function(obj) {
    dput(x = get(obj[1]), file = obj[2])
  } ))
}