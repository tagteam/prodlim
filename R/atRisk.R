atRisk <- function(x,
                   newdata,
                   times,
                   line,
                   col,
                   interspace,
                   cex,
                   labels,
                   pos,
                   adj,
                   dist,
                   adjust.labels=TRUE,
                   ...){
  if (missing(times)) times <- seq(0,x$maxtime,x$maxtime/10)
  if (x$model=="competing.risks")
    px <- lifeTab(object=x,times=times,newdata=newdata,stats=NULL)[[1]]
  else if (x$model=="survival"){
    px <- lifeTab(object=x,times=times,newdata=newdata,stats=NULL)
  }
  if (is.matrix(px) || is.data.frame(px))
    sumx <- lapply(data.frame(px)[,grep("n.risk",colnames(px)),drop=FALSE],function(x)x)
  else
    sumx <- lapply(px,function(v){
      u <- v[,grep("n.risk",colnames(v)),drop=FALSE]
      if (NCOL(u)>1){
        ulist <- lapply(1:NCOL(u),function(i)u[,i])
        names(ulist) <- colnames(u)
        ulist
      }
      else
        u
    })
  if (is.list(sumx[[1]]))
    sumx <- unlist(sumx,recursive=FALSE)
  if (all(sapply(sumx,NCOL))==1)
    nlines <- length(sumx)
  if (missing(line))
    line <- par()$mgp[2] + dist +
      (0:(2*nlines-1)) *interspace
  if (missing(cex)) cex <- 1
  if (missing(pos)) pos <- min(times)
  if (missing(adj)) adj <- 1.5
  
  if (missing(labels))
    if (length(names(sumx)==nlines))
      labels <- paste("[",names(sumx),"]",sep="")
    else
      labels <- c("No.   \nat-risk",rep("",nlines-1))
  # labeling the no. at-risk below plot
  # --------------------------------------------------------------------
  if (is.null(adjust.labels) || adjust.labels==TRUE){
    labels <- format(labels,justify="left")}
  if (length(col)==nlines/2) ## 1 cluster level
    col <- c(col,col)
  lapply(1:nlines,function(y){
    mtext(text=as.character(sumx[[y]]),
          side=1,
          at=times,
          line=rep(line[y],length(times)),
          col=rep(col[y],length(times)),
          cex=cex,
          outer=FALSE,
          xpd=NA,
          ...)
    mtext(text=labels[y],
          side=1,
          at=pos,
          col=col[y],
          line=line[y],
          adj=adj,
          cex=cex,
          outer=FALSE,
          xpd=NA,
          ...)
  })
}
