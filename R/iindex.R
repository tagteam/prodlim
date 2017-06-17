iindex <- function (L,R,grid) {
  stopifnot((length(grid)>0)
            & (length(L)>0)
            & (length(R)>0))
  stopifnot(is.numeric(c(L,R,grid)))
  N <- length(L)
  NS <- length(grid)
  ind <- .C("iindexSRC",
            index = integer(N*NS),
            strata = integer(NS-1),
            as.double(L),
            as.double(R),
            as.double(grid),
            as.integer(N),
            as.integer(NS),
            PACKAGE="prodlim")
  strata <- ind$strata
  index <- ind$index[1:max(strata)]
  list(iindex=index,imax=strata)
}
