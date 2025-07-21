# To output what R packages are loaded and used for a particular session
# Note that this function was taken and modified from the Data Quality Check Rmd
# The (other) R packages are concatenated and saved to the clipboard so they can be pasted into a report
otherPkgs <- names(sessionInfo()[["otherPkgs"]])
pkgs_vers <- vector("character", length(otherPkgs))
for (i in seq_along(otherPkgs)) {
  
  pkgs_vers[[i]] <- sessionInfo()[["otherPkgs"]][[otherPkgs[i]]][["Version"]]
  
}

print_pkgs <- vector("character", length(otherPkgs))
for (i in seq_along(otherPkgs)) {
  
  print_pkgs[[i]] <- paste0(otherPkgs[[i]], " (v", pkgs_vers[[i]], ")")
  
}

clipr::write_clip(paste0(print_pkgs, collapse = ", "))
