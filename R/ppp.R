.First.lib <- function(lib,pkg) {
  library.dynam("prodlim",pkg,lib)
}
.Last.lib <- function(lib)
  library.dynam.unload("prodlim",lib)
