library(prodlim)
library(pseudo)

# ------------------------------survival------------------------------

set.seed(177)
N <- 8888
ddd <- SimSurv(N)
ttt <- c(12,15,33)
fff <- prodlim(Hist(time,status)~1,data=ddd)

Rprof()
jack2 <- jackknife(fff,times=ttt)
Rprof(NULL)
summaryRprof()


system.time(jack <- pseudosurv(time=ddd$time,event=ddd$status,tmax=ttt))
system.time({jack2 <- jackknife.survival(fff,times=ttt)})
jack2 <- jackknife.survival(fff,times=ttt)

all(round(jack[,-c(1,2)],13)==round(jack2,13))

# ------------------------------comprisk------------------------------

set.seed(199)
N <- 2000
ddd <- prodlim:::SimCompRisk(N)
ttt <- c(12,15,33)
fff <- prodlim(Hist(time,cause)~1,data=ddd)
Rprof()
jack2 <- jackknife(fff,times=ttt)
Rprof(NULL)
summaryRprof()
system.time(jack <- pseudoci(time=ddd$time,event=ddd$cause,tmax=ttt))
system.time(jack2 <- jackknife(fff,times=ttt))
system.time({jack2.cause2 <- jackknife(fff,times=ttt,cause=2)})
all(round(jack[,c(3,5,7)],7)==round(jack2,7))
all(round(jack[,c(4,6,8)],7)==round(jack2.cause2,7))


# -------------------------further tests--------------------------------

library(prodlim)
library(pseudo)
set.seed(17)
N <- 8
ddd <- data.frame(time=sample(1:N),cause=rbinom(N,2,.5),X=rbinom(N,1,.5))
ttt <- c(5)
# ttt <- ddd$time
fff <- prodlim(Hist(time,cause)~1,data=ddd)
system.time(jack <- with(ddd,pseudoci(time,cause,ttt)))
system.time({jack2 <- jackknife.competing.risks(fff,times=ttt,cause=2)})
system.time({jack2a <- jackknife.competing.risks(fff,times=ttt,useC=FALSE)})
cbind(round(jack[,3],2),round(jack2,2))
cbind(round(jack[,45],2),round(jack2,2),round(jack2a,2))


## check individual 1
all(round(jack2[,1],9)==round(jack[,1],9))
     
## check all individuals
all(sapply(1:N,function(x){
  a <- round(jack[x,],8)
  b <- round(jack2[x,],8)
  # all(a[!is.na(a)]==b[!is.na(b)])
  all(a[!is.na(a)]==b[!is.na(a)])
}))
