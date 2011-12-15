# {{{ Header

plot.prodlim <- function(x,
                         type,
                         cause=1,
                         newdata,
                         add = FALSE,
                         col,
                         lty,
                         lwd,
                         ylim,
                         xlim,
                         xlab="Time",
                         ylab,
                         legend=TRUE,
                         marktime=FALSE,
                         confint=TRUE,
                         automar,
                         atrisk=ifelse(add,FALSE,TRUE),
                         timeOrigin=0,
                         axes=TRUE,
                         background=TRUE,
                         percent=TRUE,
                         minAtrisk=0,
                         ...){

  # }}}
  # {{{  backward compatibility
  ##   args=match.call(expand=TRUE)
  ##   args[[1]]=list
  allArgs <- match.call()
  if (missing(type)){
    type=allArgs[[match("what",names(allArgs))]]
  }
  # }}}
  # {{{  extracting a list of lines to draw

  cens.type <- x$cens.type    # uncensored, right or interval censored
  if (cens.type=="intervalCensored") confint <- FALSE
  if (cens.type=="intervalCensored") atrisk <- FALSE
  cotype <- x$covariate.type       # no, discrete, continuous or both
  model <- x$model                 # survival, competing risks or multi-state
  clusterp <- !is.null(x$clustervar)
  if (missing(type)||is.null(type)){
    type <- switch(model,"survival"="surv","competing.risks"="cuminc","multi.states"="hazard")
    if (!is.null(x$reverse) && x$reverse==TRUE && model=="survival") type <- "cuminc"
  }
  else
    type <- match.arg(type,c("surv","cuminc","hazard"))
  if (model=="competing.risks" && type=="surv") stop("To plot the event-free survival curve, please fit a suitable model: prodlim(Hist(time,status!=0)~....")
  
  if (cens.type=="intervalCensored")
    plot.times <- sort(unique(x$time[2,]))
  else{
    plot.times <- sort(unique(x$time))
    if (plot.times[1]>timeOrigin) plot.times <- c(timeOrigin,plot.times)
  }
  if (length(x$clustervar)>0)
    nRisk <- x$n.risk[,1]
  else
    nRisk <- x$n.risk
  if (minAtrisk>0 && any(nRisk<=minAtrisk)){

    if (all(nRisk<=minAtrisk)){
      return(plot(0,0,type="n",xlim=c(0, max(plot.times)),ylim=c(0, 1),axes=FALSE))
    }
    criticalTime <- min(x$time[nRisk<=minAtrisk])
    plot.times <- plot.times[plot.times<criticalTime]
  }
  if (missing(newdata)) newdata <- x$X
  if (NROW(newdata)>10) newdata <- newdata[c(1,round(median(1:NROW(newdata))),NROW(newdata)),,drop=FALSE]
  
  if (length(cause)!=1){
    warning("Currently only the cumulative incidence of a single cause can be plotted in one go. Use argument add=TRUE to add the lines of the other causes. For now I use the first cause")
    cause <- cause[1]
  }
  ## Y <- predict(x,times=plot.times,newdata=newdata,level.chaos=1,type=type,cause=cause,mode="list")
  startValue=ifelse(type=="surv",1,0)
  stats=list(c(type,startValue))
  if (model=="survival" && type=="cuminc") {
    startValue=1
    stats=list(c("surv",startValue))
  }
  if (confint==TRUE)
    stats=c(stats,list(c("lower",startValue),c("upper",startValue)))
  if (x$cens.type=="intervalCensored")
    stop("FIXME")
  sumX <- lifeTab(x,
                  times=plot.times,
                  newdata=newdata,
                  stats=stats,
                  percent=FALSE)

  if (model=="competing.risks") sumX=sumX[[cause]]
  if (cotype==1) sumX=list(sumX)
  if (model=="survival" && type=="cuminc"){
    Y <- lapply(sumX,function(x)1-x[,"surv"])
    nlines <- length(Y)}
  else{
    if (NROW(newdata)==1) sumX <- list(sumX)
    Y <- lapply(sumX,function(x)x[,type])
    nlines <- length(Y)
  }

  # }}}
  # {{{  getting default arguments for plot, atrisk, axes, legend, confint, marktime

  if (missing(ylab)) ylab <- switch(type,"surv"=ifelse(x$reverse==TRUE,"Censoring probability","Survival probability"),"cuminc"="Cumulative incidence","hazard"="Cumulative hazard")
  if (missing(xlab)) xlab <- "Time"
  if (missing(xlim)) xlim <- c(0, max(plot.times))
  if (missing(ylim)) ylim <- c(0, 1)
  if (missing(lwd)) lwd <- rep(3,nlines)
  if (missing(col)) col <- 1:nlines
  if (missing(lty)) lty <- rep(1, nlines)
  if (length(lwd) < nlines) lwd <- rep(lwd, nlines)
  if (length(lty) < nlines) lty <- rep(lty, nlines)
  if (length(col) < nlines) col <- rep(col, nlines)
  
  background.DefaultArgs <- list(xlim=xlim,ylim=ylim,horizontal=seq(0,1,.25),vertical=NULL,bg="white",fg="gray88")
  axis1.DefaultArgs <- list()
  axis2.DefaultArgs <- list(at=seq(0,1,.25),side=2)
  lines.DefaultArgs <- list(type="s")
  plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,xlab = xlab,ylab = ylab)
  marktime.DefaultArgs <- list(x=Y,nlost=lapply(sumX,function(x)x[,"n.lost"]),times=plot.times,pch="I",col=col)
  atrisk.DefaultArgs <- list(x=x,newdata=newdata,interspace=1,dist=.3,col=col,times=seq(0,min(x$maxtime,xlim[2]),min(x$maxtime,xlim[2])/10))
  legend.DefaultArgs <- list(legend=names(Y),
                             lwd=lwd,
                             col=col,
                             lty=lty,
                             cex=1.5,
                             bty="n",
                             y.intersp=1.3,
                             trimnames=TRUE,
                             x="topright")
  confint.DefaultArgs <- list(x=x,newdata=newdata,type=type,citype="shadow",times=plot.times,cause=cause,density=55,col=col[1:nlines],lwd=rep(2,nlines),lty=rep(3,nlines))

  # }}}
  # {{{  backward compatibility

  if (match("legend.args",names(args),nomatch=FALSE)){
    legend.DefaultArgs <- c(args[[match("legend.args",names(args),nomatch=FALSE)]],legend.DefaultArgs)
    legend.DefaultArgs <- legend.DefaultArgs[!duplicated(names(legend.DefaultArgs))]
  }
  if (match("confint.args",names(args),nomatch=FALSE)){
    confint.DefaultArgs <- c(args[[match("confint.args",names(args),nomatch=FALSE)]],confint.DefaultArgs)
    confint.DefaultArgs <- confint.DefaultArgs[!duplicated(names(confint.DefaultArgs))]
  }
  if (match("atrisk.args",names(args),nomatch=FALSE)){
    atrisk.DefaultArgs <- c(args[[match("atrisk.args",names(args),nomatch=FALSE)]],atrisk.DefaultArgs)
    atrisk.DefaultArgs <- atrisk.DefaultArgs[!duplicated(names(atrisk.DefaultArgs))]
  }
  smartA <- SmartControl(call=  list(...),
                         keys=c("plot","lines","atrisk","legend","confint","background","marktime","axis1","axis2"),
                         ignore=c("x","type","cause","newdata","add","col","lty","lwd","ylim","xlim","xlab","ylab","legend","marktime","confint","automar","atrisk","timeOrigin","percent","axes","atrisk.args","confint.args","legend.args"),
                         defaults=list("plot"=plot.DefaultArgs,"atrisk"=atrisk.DefaultArgs,"lines"=lines.DefaultArgs,"legend"=legend.DefaultArgs,"confint"=confint.DefaultArgs,"marktime"=marktime.DefaultArgs,"background"=background.DefaultArgs,"axis1"=axis1.DefaultArgs,"axis2"=axis2.DefaultArgs),
                         forced=list("plot"=list(axes=FALSE),"axis1"=list(side=1)),
                         ignore.case=TRUE,
                         replaceDefaults=FALSE,
                         verbose=TRUE)

# }}}
  # {{{  setting margin parameters 

  if (atrisk==TRUE){
    oldmar <- par()$mar
    if (missing(automar) || automar==T){
      ##        bottomMargin =  margin line (in 'mex' units) for xlab
      ##                        + distance of xlab from xaxis
      ##                        + distance of atrisk numbers from xlab
      ##                        + number of atrisk lines
      ##                        + one extra line below the bottom number atrisk line
      ##      leftSideMargin =  margin line + atrisk.lab
      bottomMargin <- par()$mgp[2] + smartA$atrisk$dist+ ifelse(clusterp,2,1)*nlines + 1
      newmar <- par()$mar + c(bottomMargin,0,0,0)
      par(mar=newmar)
    }
  }

  # }}}
  # {{{  plot and backGround
  if (!add) {
    do.call("plot",smartA$plot)
    ##     if (background==TRUE && match("bg",names(smartA$background),nomatch=FALSE)){
    ## par(bg=smartA$background$bg)
    ##     }
    if (background==TRUE){
      do.call("backGround",smartA$background)
    }
  }
  # }}}
  # {{{  axes

  if (!add) {
    if (axes){
      do.call("axis",smartA$axis1)
      if (percent & is.null(smartA$axis2$labels))
        smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
      do.call("axis",smartA$axis2)
    }
  }
  if (atrisk==TRUE) par(mar=oldmar) ## reset

  # }}}
  # {{{  pointwise confidence intervals
  if (confint==TRUE) {
    ## if (verbose==TRUE){print(smartA$confint)}
    do.call("confInt",smartA$confint)
  }

  # }}}
  # {{{  adding the lines 
  lines.type <- smartA$lines$type
  nix <- lapply(1:nlines, function(s) {
    lines(x = plot.times,
          y = Y[[s]],
          type = lines.type,
          col = col[s],
          lty = lty[s],
          lwd = lwd[s])
  })
# }}}
  # {{{  marks at the censored times

  if (marktime==TRUE){
    if (model %in% c("survival","competing.risks")){
      do.call("markTime",smartA$marktime)
    }
    else{
      message("Marking the curves at censored times is not yet available for multi-state models.")
    }
  }

# }}}
  # {{{  adding the no. of individuals at risk

  if (atrisk==TRUE && !add){
    if (hit <- match("at",names(smartA$atrisk),nomatch=FALSE)){
      if (match("atrisk.times",names(list(...)),nomatch=FALSE)){
        warning("Atrisk argument clash: remove either 'atrisk.at' or 'atrisk.times'.")
      }
      else{
        names(smartA$atrisk)[hit] <- "times"
        smartA$atrisk <- smartA$atrisk[!duplicated(names(smartA$atrisk))]
      }
    }
    do.call("atRisk",smartA$atrisk)
  }
  # }}}
  # {{{  legend

  if(legend==TRUE && !add && !is.null(names(Y))){

    if (smartA$legend$trimnames==TRUE){
      smartA$legend$legend <- sapply(strsplit(names(Y),"="),function(x)x[[2]])
      smartA$legend$title <- unique(sapply(strsplit(names(Y),"="),function(x)x[[1]]))
    }
    smartA$legend <- smartA$legend[-match("trimnames",names(smartA$legend))]
    save.xpd <- par()$xpd
    par(xpd=TRUE)
    do.call("legend",smartA$legend)
    par(xpd=save.xpd)
  }

# }}}
  invisible(x)

}
