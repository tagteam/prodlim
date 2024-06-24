library(testthat)
library(prodlim)
library(data.table)
context("Prodlim")

test_that("Competing risks: asking for case for event that does not exist in data",{
    ##
    set.seed(10)
    d <- SimSurv(10)
    setDT(d)
    d[,event:=factor(event,levels=c(0,1),labels=c("0","2"))]
    f <- prodlim(Hist(time,event)~X1,data=d)
    predict(f,cause="2",times=8,newdata=data.frame(X1=1))
    expect_error(predict(f,cause="1",times=8,newdata=data.frame(X1=1)))
    set.seed(10)
    dd <- SimCompRisk(20)
    F <- prodlim(Hist(time,event)~X1,data=dd)
    predict(F,cause="1",times=4,newdata=data.frame(X1=0:1))
    expect_equal(lapply(predict(F,cause=2,times=4,newdata=data.frame(X1=0:1)),round,4),list(`X1=0`=0.0714,`X1=1`=0))
    expect_error(predict(F,cause=3,times=4,newdata=data.frame(X1=0:1)))
    expect_error(summary(F,cause=3))
    expect_error(plot(F,cause=3))
})

test_that("strata",{
    ## bug in version 1.5.1
    d <- data.frame(time=1:3,status=c(1,0,1),a=c(1,9,9),b=factor(c(0,1,0)))
    expect_output(print(prodlim(Hist(time,status)~b+factor(a),data=d)))
})


test_that("prodlim: print and summary  ",{
    library(lava)
    library(riskRegression)
    library(survival)
    # lava some data
    m <- crModel()
    addvar(m) <- ~X1+X2+X3+X4+X5+X6
    distribution(m,"X3") <- binomial.lvm()
    distribution(m,"X4") <- normal.lvm(mean=50,sd=10)
    distribution(m,"eventtime1") <- coxWeibull.lvm(scale=1/200)
    distribution(m,"censtime") <- coxWeibull.lvm(scale=1/1000)
    m <- categorical(m,K=4,eventtime1~X5,beta=c(1,0,0,0),p=c(0.1,0.2,0.3))
    m <- categorical(m,K=3,eventtime1~X1,beta=c(2,1,0),p=c(0.3,0.2))
    regression(m,to="eventtime1",from=c("X2","X4")) <- c(0.3,0)
    regression(m,to="eventtime2",from=c("X2","X4")) <- c(0.6,-0.07)
    set.seed(17)
    d <- sim(m,200)
    d$X1 <- factor(d$X1,levels=c(0,1,2),labels=c("low survival","medium survival","high survival"))
    d$X5 <- factor(d$X5,levels=c("0","1","2","3"),labels=c("one","two","three","four"))
    d$Event <- factor(d$event,levels=c("0","1","2"),labels=c("0","cause-1","cause-2"))
    d$status <- 1*(d$event!=0)
    # marginal Kaplan-Meier
    s0 <- prodlim(Hist(time,status)~1,data=d)
    expect_output(print(s0))
    a <- summary(s0,intervals=TRUE)
    b <- summary(s0,intervals=TRUE,format="list")$table
    attributes(a) <- attributes(b)
    expect_equal(a,b)
    stats::predict(s0,times=1:10)
    # 
    # stratified Kaplan-Meier vs subsetting
    #
    su <- prodlim(Hist(time,status)~1,data=d,subset=d$X1=="medium survival")
    s1 <- prodlim(Hist(time,status)~X1,data=d)
    S1 <- summary(s1,intervals=0,newdata=data.frame(X1=c("medium survival")))
    S2 <- summary(s1,intervals=0,newdata=data.frame(X1=c("medium survival","high survival","low survival")))
    S1a <- S2[S2$X1=="medium survival",]
    rownames(S1a) <- 1:NROW(S1a)
    expect_equal(S1,S1a)
    stats::predict(s1,times=0:10,newdata=data.frame(X1=c("medium survival","low survival","high survival")))
    suppressWarnings(A <- summary(su,intervals=0,times=S1$time))
    S1$X1 <- NULL
    attributes(A) <- attributes(S1)
    all.equal(A,S1)
    #
    # Beran estimator
    #
    s2 <- prodlim(Hist(time,status)~X2,data=d)
    expect_output(print(s2))
    suppressWarnings(summary(s2,intervals=TRUE))
    stats::predict(s2,times=0:10,newdata=data.frame(X2=quantile(d$X2)))
    #
    # Kaplan-Meier: two strata variable
    #
    s1a <- prodlim(Hist(time,status)~X1+X3,data=d)
    expect_output(print(s1a))
    stats::predict(s1a,times=0:10,newdata=expand.grid(X1=levels(d$X1),X3=unique(d$X3)))
    summary(s1a,intervals=TRUE,newdata=s1a$X[c(2,3),])
    #
    # Stratified Beran estimator
    #
    s3 <- prodlim(Hist(time,status)~X1+X2,data=d)
    expect_output(print(s3))
    stats::predict(s3,times=0:10,newdata=expand.grid(X1=levels(d$X1),X2=c(quantile(d$X2,0.05),median(d$X2))))
    suppressWarnings(summary(s3,intervals=TRUE))
    #
    # marginal Aalen-Johansen competing risks
    #
    f0 <- prodlim(Hist(time,event)~1,data=d)
    expect_output(print(f0))
    summary(f0,intervals=TRUE)
    stats::predict(f0,times=1:10)
    #
    # stratified Aalen-Johansen competing risks
    #
    f1 <- prodlim(Hist(time,event)~X1,data=d)
    expect_output(print(f1))
    summary(f1,intervals=TRUE,newdata=data.frame(X1=c("medium survival","high survival","low survival")))
    summary(f1,intervals=TRUE,cause=2,newdata=data.frame(X1=c("medium survival","low survival")))
    stats::predict(f1,times=0:10,newdata=data.frame(X1=c("medium survival","low survival","high survival")))
    #
    # Beran-Aalen-Johansen competing risks
    #
    f2 <- prodlim(Hist(time,event)~X2,data=d)
    expect_output(print(f2))
    suppressWarnings(summary(f2,intervals=TRUE))
    stats::predict(f2,times=0:10,newdata=data.frame(X2=quantile(d$X2)))
    #
    # 2-strata Aalen-Johansen competing risks
    #    
    f1a <- prodlim(Hist(time,event)~X1+X3,data=d)
    expect_output(print(f1a))
    summary(f1a,intervals=TRUE)
    stats::predict(f1a,times=0:10,newdata=expand.grid(X1=levels(d$X1),X3=unique(d$X3)))
    #
    # stratified Beran-Aalen-Johansen competing risks
    #    
    f3 <- prodlim(Hist(time,event)~X1+X2,data=d)
    expect_output(print(f3))
    suppressWarnings(summary(f3,intervals=TRUE))
    stats::predict(f3,times=0:10,newdata=expand.grid(X1=levels(d$X1),X2=c(quantile(d$X2,0.05),median(d$X2))))
})

test_that("prodlim vs survfit",{
    data(pbc)
    prodlim.0 <- prodlim(Hist(time,status!=0)~1,data=pbc)
    survfit.0 <- survfit(Surv(time,status!=0)~1,data=pbc)
    ttt <- sort(unique(pbc$time)[pbc$status!=0])
    ttt <- ttt[-length(ttt)]
    sum0.s <- summary(survfit.0,times=ttt)
    ## There is arounding issue with survfit:
    #---------------------------------------------------------------------------------------
    ## testdata <- data.frame(time=c(16.107812,3.657545,1.523978),event=c(0,1,1))
    ## sum0 <- summary(survfit(Surv(time,event)~1,data=testdata),times=sort(testdata$time))
    ## testdata$timeR <- round(testdata$time,1)
    ## sum1 <- summary(survfit(Surv(timeR,event)~1,data=testdata),times=sort(testdata$time))
    ## sum0
    ## sum1
    ## sum0 != sum1
    ## summary(survfit.0,times=c(0,0.1,0.2,0.3))
    #---------------------------------------------------------------------------------------
    result.survfit <- data.frame(time=sum0.s$time,n.risk=sum0.s$n.risk,n.event=sum0.s$n.event,surv=sum0.s$surv,std.err=sum0.s$std.err,lower=sum0.s$lower,upper=sum0.s$upper)
    result.prodlim <- data.frame(summary(prodlim.0,times=ttt)[,c("time","n.risk","n.event","n.lost","surv","se.surv","lower","upper")])
    ## cbind(result.survfit[,c("time","n.risk","n.event","surv")],result.prodlim[,c("time","n.risk","n.event","surv")])
    a <- round(result.survfit$surv,8)
    b <- round(result.prodlim$surv[!is.na(result.prodlim$se.surv)],8)
    expect_equal(a,b)
    expect_equal(round(result.survfit$std.err,8),round(result.prodlim$se.surv[!is.na(result.prodlim$se.surv)],8))
    pbc <- pbc[order(pbc$time,-pbc$status),]
    set.seed(17)
    boot <- sample(1:NROW(pbc),size=NROW(pbc),replace=TRUE)
    boot.weights <- table(factor(boot,levels=1:NROW(pbc)))
    s1 <- prodlim(Hist(time,status>0)~1,data=pbc,caseweights=boot.weights)
    ## plot(s1,col=1,confint=FALSE,lwd=8)
    s2 <- prodlim(Hist(time,status>0)~1,data=pbc[sort(boot),])
    ## plot(s2,add=TRUE,col=2,confint=FALSE,lwd=3)
    expect_equal(summary(s1,intervals=1,times=seq(0,3500,500)),summary(s2,intervals=1,times=seq(0,3500,500)))
})


test_that("weigths, subset and smoothing",{
    d <- SimSurv(100)
    f1 <- prodlim(Hist(time,status)~X2,data=d)
    f2 <- prodlim(Hist(time,status)~X2,data=d,caseweights=rep(1,100))
    expect_equal(f1$surv,f2$surv)
    d <- SimSurv(100)
    d <- data.frame(d, group = c(rep(1, 70), rep(0,30)))
    f1a <- prodlim(Hist(time,status)~X2,data=d, caseweights = rep(1, 100), subset = d$group==1,bandwidth=0.1)
    f1b <- prodlim(Hist(time,status)~X2,data=d[d$group==1, ], caseweights = rep(1, 100)[d$group==1], bandwidth=0.1)
    f1a$call <- f1b$call
    expect_equal(f1a,f1b)
    f1 <- prodlim(Hist(time,status)~X1,data=d, subset = d$group==1)
    f2 <- prodlim(Hist(time,status)~X1,data=d,caseweights=d$group)
    expect_equal(unique(f1$surv),unique(f2$surv))
    expect_equal(predict(f1,newdata = d[1, ], times = 5),
                 predict(f2, newdata = d[1, ], times = 5))
})

test_that("weights and delay",{
    library(survival)
    library(survey)
    ## library(SmoothHazard)
    library(etm)
    pbc <- pbc[order(pbc$time,-pbc$status),]
    ## pbc$randprob<-fitted(biasmodel)
    ## pbc$randprob <- as.numeric(pbc$sex=="m")+0.1
    set.seed(17)
    pbc$randprob <- abs(rnorm(NROW(pbc)))
    dpbc <- svydesign(id=~id, weights=~randprob, strata=NULL, data=pbc)
    survey.1<-svykm(Surv(time,status>0)~1, design=dpbc)
    ## plot(survey.1,lwd=8)
    prodlim.1 <- prodlim(Hist(time,status>0)~1,data=pbc,caseweights=pbc$randprob)
    ## plot(prodlim.1,add=TRUE,col=2,confint=FALSE)
    pbc$entry <- round(pbc$time/5)
    survfit.delay <- survfit(Surv(entry,time,status!=0)~1,data=pbc)
    prodlim.delay <- prodlim(Hist(time,status!=0,entry=entry)~1,data=pbc)
    a <- summary(survfit.delay)
    b <- summary(prodlim.delay,times=a$time)
    expect_equal(a$surv,b$surv)
    ## plot(survfit.delay,lwd=8)
    ## plot(prodlim.delay,lwd=4,col=2,add=TRUE,confint=FALSE)
    pbc0 <- pbc
    pbc0$entry <- round(pbc0$time/5)
    survfit.delay.edema <- survfit(Surv(entry,time,status!=0)~edema,data=pbc0)
    ## survfit.delay.edema.0.5 <- survfit(Surv(entry,time,status!=0)~1,data=pbc0[pbc0$edema==0.5,])
    prodlim.delay.edema <- prodlim(Hist(time,status!=0,entry=entry)~edema,data=pbc0)
    ## prodlim.delay.edema.0.5 <- prodlim(Hist(time,status!=0,entry=entry)~1,data=pbc0[pbc0$edema==0.5,])
    ## plot(survfit.delay.edema,conf.int=FALSE,col=1:3,lwd=8)
    ## plot(prodlim.delay.edema,add=TRUE,confint=FALSE,col=c("gray88","orange",5),lwd=4)
    ## a <- summary(survfit.delay.edema)
    ## b <- summary(prodlim.delay.edema,times=a$time)
    ## expect_equal(a$surv,b$surv)    
    data(abortion)
    cif.ab.etm <- etmCIF(Surv(entry, exit, cause != 0) ~ 1,abortion,etype = cause,failcode = 3)
    cif.ab.prodlim <- prodlim(Hist(time=exit, event=cause,entry=entry) ~ 1,data=abortion)
    ## plot(cif.ab.etm,lwd=8,col=3)
    ## plot(cif.ab.prodlim,add=TRUE,lwd=4,col=5,cause=3)
    x <- prodlim(Hist(time=exit, event=cause,entry=entry) ~ 1,data=abortion)
    x0 <- etmCIF(Surv(entry, exit, cause != 0) ~ 1,data=abortion,etype = cause)
    ## graphics::par(mfrow=c(2,2))
    cif.ab.etm <- etmCIF(Surv(entry, exit, cause != 0) ~ 1,abortion,etype = cause,failcode = 3)
    cif.ab.prodlim <- prodlim(Hist(time=exit, event=cause,entry=entry) ~ 1,data=abortion)
                                        # cause 3
    ## plot(cif.ab.etm, ci.type = "bars", pos.ci = 24, col = c(1, 2), lty = 1,which.cif=3,lwd=8)
    ## plot(cif.ab.prodlim,add=TRUE,cause=3,confint=TRUE,col=2)
                                        # cause 2
    ## plot(cif.ab.etm, ci.type = "bars", pos.ci = 24, col = c(1, 2), lty = 1,which.cif=2,lwd=8)
    ## plot(cif.ab.prodlim,add=TRUE,cause=2,confint=TRUE,col=2)
                                        # cause 1
    ## plot(cif.ab.etm, ci.type = "bars", pos.ci = 24, col = c(1, 2), lty = 1,which.cif=1,lwd=8)
    ## plot(cif.ab.prodlim,add=TRUE,cause=1,confint=TRUE,col=2)
    cif.ab.etm <- etmCIF(Surv(entry, exit, cause != 0) ~ group,abortion,etype = cause,failcode = 3)
    names(cif.ab.etm[[1]])
    head(cbind(cif.ab.etm[[1]]$time,cif.ab.etm[[1]]$n.risk))
    cif.ab.prodlim <- prodlim(Hist(time=exit, event=cause,entry=entry) ~ group,data=abortion)
    ## plot(cif.ab.etm, ci.type = "bars", pos.ci = 24, col = c(1, 2), lty = 1, curvlab = c("Control", "Exposed"),lwd=8)
    ## plot(cif.ab.prodlim,add=TRUE,cause=3,confint=FALSE,col="yellow")
    testdata <- data.frame(entry=c(1,5,2,8,5),exit=c(10,6,4,12,33),event=c(0,1,0,1,0))
    cif.test.etm <- etmCIF(Surv(entry, exit, event) ~ 1,data=testdata,etype = event,failcode = 1)
    cif.test.survival <- survfit(Surv(entry, exit, event) ~ 1,data=testdata)
    cif.test.prodlim <- prodlim(Hist(exit,event,entry=entry)~1,data=testdata)
    ## plot(cif.test.etm, ci.type = "bars", pos.ci = 24, lwd=5)
    ## plot(cif.test.etm, ci.type = "bars", pos.ci = 24, lwd=5)
    ## plot(cif.test.prodlim,add=TRUE,cause=2,col=2,confint=TRUE,type="cuminc")
    ## simulate data from an illness-death model
    ## mod <- idmModel(K=10,schedule=0,punctuality=1)
    ## regression(mod,from="X",to="lifetime") <- log(2)
    ## regression(mod,from="X",to="waittime") <- log(2)
    ## regression(mod,from="X",to="illtime") <- log(2)
    ## set.seed(137)
    ## we round the event times to have some ties
    ## testdata <- round(sim(mod,250),1)
    ## the data enter with delay into the intermediate state (ill)
    ## thus, to estimate the absolute risk cumulative incidence of
    ## the absorbing state (death) after illness we 
    ## have left-truncated data
    ## illdata <- testdata[testdata$illstatus==1,]
    ## illdata <- illdata[order(illdata$lifetime,-illdata$seen.exit),]
    ## sindex(jump.times=illdata$illtime,eval.times=illdata$lifetime)
    ## F <- prodlim(Hist(lifetime,status,entry=illtime)~1,data=illdata[1:5,])
    ## f <- survfit(Surv(illtime,lifetime,status)~1,data=illdata[1:5,],type="kaplan-meier")
    ## survfit.delayed.ill <- survfit(Surv(illtime,lifetime,seen.exit)~1,data=illdata)
    ## prodlim.delayed.ill <- prodlim(Hist(lifetime,seen.exit,entry=illtime)~1,data=illdata)
    ## plot(survfit.delayed.ill,lwd=5)
    ## plot(prodlim.delayed.ill,lwd=2,col=2,add=TRUE)
})

test_that("left truncation: survival",{
    library(prodlim)
    library(data.table)
    library(survival)
    dd <- data.table(entry=c(1,1,56,1,1,225,277,1647,1,1),
                     time=c(380,46,217,107,223,277,1638,2164,45,40),
                     status=c(1,0,1,1,0,0,0,1,0,1))
    ## --------------------------------------------------------------
    ## by convention in case of ties 
    ## entry happens after events and after censoring
    ## --------------------------------------------------------------
    prodlim.delayed <- prodlim(Hist(time,status,entry=entry)~1,data=dd)
    data.table(time=prodlim.delayed$time,n.risk=prodlim.delayed$n.risk,n.event=prodlim.delayed$n.event,n.lost=prodlim.delayed$n.lost)
    summary(prodlim.delayed,times=c(0,10,56,267,277,1000,2000))    
    survfit.delayed <- survfit(Surv(entry,time,status)~1,data=dd)
    summary(prodlim.delayed,times=c(0,10,40),intervals=TRUE)
    summary(survfit.delayed,times=c(0,1,10,40,50))
    summary.survfit.delayed <- summary(survfit.delayed,times=c(0,10,56,267,277,1000,2000))
    summary.prodlim.delayed <- summary(prodlim.delayed,times=c(0,10,56,267,277,1000,2000),intervals=1)
    expect_equal(as.numeric(summary.survfit.delayed$surv),
                 as.numeric(unlist(summary.prodlim.delayed[,"surv"])))
    ## FIXME: lifetab does not handle delayed entry
    ##        and shows wrong numbers at risk before the
    ##        first event time
    ## expect_equal(as.numeric(summary.survfit.delayed$n.risk),
    ## as.numeric(summary.prodlim.delayed[,"n.risk"]))
    
})


context("Clustered survival data")
test_that("clustersurv",{
    library(prodlim)
    ## if (!is.function("cluster")) cluster <- function(x)x
    clusterTestData <- data.frame(midtimeX=1:8,eventX=c(0,"pn","pn",0,0,0,0,0),patientid=c(1,1,2,2,3,3,4,4),AnyCrownFracture=c(1,1,1,1,2,2,2,2))
    a <- prodlim(Hist(midtimeX,eventX=="pn")~cluster(patientid)+AnyCrownFracture,data=clusterTestData)
    b <- prodlim(Hist(midtimeX,eventX=="pn")~cluster(patientid),data=clusterTestData[clusterTestData$AnyCrownFracture==1,])
    c <- prodlim(Hist(midtimeX,eventX=="pn")~cluster(patientid),data=clusterTestData,subset=clusterTestData$AnyCrownFracture==1)
    d <- prodlim(Hist(midtimeX,eventX=="pn")~1,data=clusterTestData[clusterTestData$AnyCrownFracture==2,])
    expect_equal(round(as.numeric(summary(a)[AnyCrownFracture==1][["se.surv"]]),5),c(0,0.20951,0.10476,0.10476,NA,NA,NA,NA))
    expect_equal(summary(b), summary(c))
})

context("Construction of pseudovalues")
test_that("pseudo",{
    library(prodlim)
    library(pseudo)
    # comparison to pseudoci
    # make sure we get the same
    # results with both packages
    set.seed(17)
    N <- 80
    ddd <- SimCompRisk(80)
    ttt <- c(3,5,10)
    # ttt <- ddd$time
    fff <- prodlim(Hist(time,event)~1,data=ddd)
    system.time(jack <- with(ddd,pseudoci(time,event,ttt)))
    system.time({jack2 <- jackknife.competing.risks(fff,times=ttt,cause=1)})
    ## check individual 2
    expect_true(all(round(jack2[,2],9)==round(jack[[3]]$cause1[,2],9)))
    ## check all individuals
    expect_true(all(sapply(1:N,function(x){
        a <- round(jack[[3]]$cause1[x,],8)
        b <- round(jack2[x,],8)
        # all(a[!is.na(a)]==b[!is.na(b)])
        all(a[!is.na(a)]==b[!is.na(a)])
    })))
    ## the pseudoci function seems only slightly slower
    ## for small sample sizes (up to ca. 200) but
    ## much slower for large sample sizes:
    set.seed(17)
    N <- 80
    ddd <- SimCompRisk(80)
    ttt <- c(3,5,10)
    # ttt <- ddd$time
    fff <- prodlim(Hist(time,event)~1,data=ddd)
    system.time(jack <- with(ddd,pseudoci(time,event,ttt)))
    system.time({jack2 <- jackknife.competing.risks(fff,times=ttt,cause=1)})
    expect_true(all(round(jack2[,1],9)==round(jack$pseudo$cause1[,1],9)))
    ## set.seed(17)
    ## N <- 2000
    ## ddd <- SimCompRisk(2000)
    ## ttt <- c(3,5,10)
    ## fff <- prodlim(Hist(time,event)~1,data=ddd)
    ## a <- system.time(jack <- with(ddd,pseudoci(time,event,ttt)))
    ## b <- system.time({jack2 <- jackknife.competing.risks(fff,times=ttt,cause=1)})
})
