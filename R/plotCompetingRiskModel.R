plotCompetingRiskModel <- function(labels,horizontal=TRUE,...){
  if (missing(labels)) labels <- c("Disease\nfree","Cause1","Cause2")
  nTrans <- length(labels)-1
  crHist <- Hist(time=1:nTrans,event=list(from=rep("1",nTrans),to=labels[-1]))
  if (horizontal==TRUE){
    columns <- c(1,rep(3,nTrans))
    names(columns) <- paste("box",1:length(labels),".column",sep="")
    ## do not specify the row of box 1
    ## this box will be centered in column 1
    rows <- 1:nTrans
    names(rows) <- paste("box",2:length(labels),".row",sep="")
    ncol <- 3
    nrow <- nTrans
  }
  else{
    nrow <- 3
    if (nTrans/2==round(nTrans/2)){
      ncol <- nTrans+1
      midCol <- ceiling(ncol/2)
      columns <- c(midCol,(1:ncol)[-midCol])
      names(columns) <- paste("box",1:length(labels),".column",sep="")
      rows <- c(1,rep(3,nTrans))
      names(rows) <- paste("box",1:length(labels),".row",sep="")
    }
    else{
      ncol <- nTrans
      columns <- c(nTrans+1/2,1:nTrans)
      names(columns) <- paste("box",1:length(labels),".column",sep="")
      rows <- c(1,rep(3,nTrans))
      names(rows) <- paste("box",2:length(labels),".row",sep="")
    }
  }
  do.call("plot.Hist",c(list(x=crHist,labels=labels,nrow=nrow,ncol=ncol,...),columns,rows))
}
