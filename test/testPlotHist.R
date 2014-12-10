# {{{  plot.Hist

## pdf("/tmp/testPlot.pdf")
pdf("~/tmp/testPlot.pdf")
library(prodlim)

## or simply create some data
SurvFrame <- data.frame(time=1:10,status=sample(0:1,10,TRUE))
SurvHist <- with(SurvFrame,Hist(time,status))
plot(SurvHist)

## two competing risks
comprisk.model <- data.frame(time=1:3,status=1:3)
CRHist <- with(comprisk.model,Hist(time,status,cens.code=2))
plot(CRHist)
plot(CRHist,box2.label="This\nis\nstate 2",arrow1.label=paste(expression(gamma[1](t))))
plot(CRHist,box3.label="Any\nLabel",arrow2.label="any\nlabel")

## change the layout
plot(CRHist,
     box1.label="Alive",box2.label="Dead\n cause 1",box3.label="Dead\n cause 2",
     arrow1.label=paste(expression(gamma[1](t)),arrow2.label=paste(expression(eta[2](t)))),
     box1.col=2,
     box2.col=3,
     box3.col=4,
     nrow=2,ncol=3,box1.row=1,box1.column=2,box2.row=2,box2.column=1,box3.row=2,box3.column=3)

## more competing risks
comprisk.model2 <- data.frame(time=1:4,status=1:4)
CRHist2 <- with(comprisk.model2,Hist(time,status,cens.code=2))
plot(CRHist2)

## illness-death models
illness.death.frame <- data.frame(time=1:4,
				  from=c("Disease\nfree","Disease\nfree","Diseased","Disease\nfree"),
				  to=c("0","Diseased","Dead","Dead"))
IDHist <- with(illness.death.frame,Hist(time,event=list(from,to)))
plot(IDHist)

## illness-death with recovery
illness.death.frame2 <- data.frame(time=1:5,from=c("Disease\nfree","Disease\nfree","Diseased","Diseased","Disease\nfree"),to=c("0","Diseased","Disease\nfree","Dead","Dead"))
IDHist2 <- with(illness.death.frame2,Hist(time,event=list(from,to)))
plot(IDHist2)

## 4 state model
x=data.frame(from=c(1,2,1,3,4),to=c(2,1,3,4,1),time=1:5)
y=with(x,Hist(time=time,event=list(from=from,to=to)))
plot(y)

# in a competing risk model multiple types of failure occur
# --------------------------------------------------------------------

compRiskHistory1 <- Hist(time=1:5,event=c("cause 1","cause 2","cause 3", 0,"cause 4"))
plot(compRiskHistory1,layout=list(nrow=3,ncol=3,box.pos=list(c(2,2),c(1,1),c(1,3),c(3,1),c(3,3))),arrow.lab.offset=5)

compRiskHistory2 <- Hist(time=1:5,event=c("Cancer","Heart attack","Car\naccident", 0,"Murder"))
plot(compRiskHistory2,layout=list(nrow=3,ncol=3,box.pos=list(c(2,2),c(1,1),c(1,3),c(3,1),c(3,3))),arrow.lab.offset=5,state.cex=1.8,xbox.rule=.1,ybox.rule=.4,tagBoxes=TRUE)

# in a multi state model events occur in a certain order
# --------------------------------------------------------------------

## illness-death model without recovery
illness.death.frame <- data.frame(time=1:4,
		    from=c("Disease-free","Disease-free",
		      "Diseased","Disease-free"),
		    to=c("0","Diseased","Dead","Dead"))

IDHist <- with(illness.death.frame,Hist(time,event=list(from,to)))
plot(IDHist,ybox.rule=4,xbox.rule=.3,state.cex=1.3,enum=TRUE,arrow.lab.side=c(-1,-1,1))


## illness-death with recovery
illness.death.frame2 <- data.frame(time=1:5,
				   from=c("Disease\nfree","Disease\nfree",
                                     "Diseased","Diseased","Disease\nfree"),
				   to=c("0","Diseased","Disease\nfree",
                                                           "Dead","Dead"))
IDHist2 <- with(illness.death.frame2,Hist(time,event=list(from,to)))
plot(IDHist2,ybox.rule=1.3,xbox.rule=.3,state.cex=2,arrow.lab.offset=c(13,13,8,10),enum=TRUE,verbose=FALSE)

## change the layout of the graphic

plot(IDHist2,ybox.rule=1.3,xbox.rule=.3,state.cex=2,enum=TRUE,verbose=FALSE,layout=list(ncol=3,nrow=2,box.pos=list(c(1,1),c(2,2),c(1,3))),arrow.lab.side=c(-1,1,-1,1),arrow.lab.offset=c(15,15,10,10))


## Bordeaux: 4 state model
x=data.frame(from=c("1","1","2","2","3"),to=c("2","4","4","3","4"),time=1:5)
y=with(x,Hist(time=time,event=list(from=from,to=to)))
plot(y,stateLabels=c("Healthy","Pre-diagnosed","Demented","Dead"),cex=1.3,ncol=3,nrow=5,box1.column=1,box2.column=2,box3.column=3,box4.column=2,box1.row=2,box2.row=1,box3.row=2,box4.row=5)

# }}}
dev.off()
## system(paste("evince ","/tmp/testPlot.pdf&"),intern=FALSE)
system(paste("evince ","~/tmp/testPlot.pdf&"),intern=FALSE)
