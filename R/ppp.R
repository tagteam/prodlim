.First.lib <- function(lib,pkg) {
  require(survival)
  library.dynam("prodlim",pkg,lib)
}
.Last.lib <- function(lib)
  library.dynam.unload("prodlim",lib)
