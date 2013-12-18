library(prodlim)
set.seed(17)
d <- SimSurv(10)
stopifnot(sum(d$status)==5)

library(prodlim)
SimSurv(10,cova=list("X1"=list("rnorm",1,1),"Z"=list("rnorm",-30,3)))
