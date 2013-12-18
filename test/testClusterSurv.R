library(prodlim)
if (!is.function("cluster")) cluster <- function(x)x
data(clusterTestData)
## clusterTestData <- data.frame(midtimeX=1:4,eventX=c(0,"pn","pn",0),patientid=c(1,1,2,2),AnyCrownFracture=c(1,1,2,2))
## clusterTestData <- data.frame(midtimeX=1:8,eventX=c(0,"pn","pn",0,0,"pn","pn",0),patientid=c(1,1,2,2,3,3,4,4),AnyCrownFracture=c(1,1,1,1,2,2,2,2))
## clusterTestData <- data.frame(midtimeX=1:8,eventX=c(0,"pn","pn",0,0,0,0,0),patientid=c(1,1,2,2,3,3,4,4),AnyCrownFracture=c(1,1,1,1,2,2,2,2))
a <- prodlim(Hist(midtimeX,eventX=="pn")~cluster(patientid)+AnyCrownFracture,data=clusterTestData)
b <- prodlim(Hist(midtimeX,eventX=="pn")~cluster(patientid),data=clusterTestData[clusterTestData$AnyCrownFracture=="Yes",])
c <- prodlim(Hist(midtimeX,eventX=="pn")~cluster(patientid),data=clusterTestData[clusterTestData$AnyCrownFracture=="No",])
d <- prodlim(Hist(midtimeX,eventX=="pn")~1,data=clusterTestData[clusterTestData$AnyCrownFracture=="2",])
summary(a)
summary(b)
summary(c)
summary(d)


