## 1.1. survival uncensored
h <- Hist(1:10)
stopifnot(attr(h,"cens.type")=="uncensored")
stopifnot(attr(h,"entry.type")==NULL)
stopifnot(attr(h,"model")=="survival")
stopifnot(colnames(h)==c("time","status"))

## 1.2. survival right censored
h <- Hist(time=1:10,event=c(0,1,0,0,0,0,0,0,0,1))
stopifnot(attr(h,"cens.type")=="uncensored")
stopifnot(attr(h,"entry.type")==NULL)
stopifnot(attr(h,"model")=="survival")
stopifnot(colnames(h)==c("time","status"))

## 1.3. survival right censored and left-truncated
h <- Hist(entry = 0:9, time=1:10,event=c(0,1,0,0,0,0,0,0,0,1))
stopifnot(attr(h,"cens.type")=="rightCensored")
stopifnot(attr(h,"entry.type")=="leftTruncated")
stopifnot(attr(h,"model")=="survival")
stopifnot(colnames(h)==c("entry", "time","status"))

## 1.4. survival: right censored and left-truncated
h <- Hist(entry = 0:9, time=1:10,event=c(0,1,0,0,0,0,0,0,0,1))
stopifnot(attr(h,"cens.type")=="rightCensored")
stopifnot(attr(h,"entry.type")=="leftTruncated")
stopifnot(attr(h,"model")=="survival")
stopifnot(colnames(h)==c("entry", "time","status"))

## 1.5. survival: interval censored and right censored and left-truncated
h <- Hist(entry = 0:9, time=list(1:10,b=2:11),event=c(0,1,0,1,0,0,0,0,0,1))
stopifnot(attr(h,"cens.type")=="intervalCensored")
stopifnot(attr(h,"entry.type")=="leftTruncated")
stopifnot(attr(h,"model")=="survival")
stopifnot(colnames(h)==c("entry", "L","R","status"))

