"predict.prodlim" <- function(object,
                              times,
                              newdata,
                              level.chaos=1,
                              type=c("surv","cuminc","list"),
                              mode="list",
                              bytime=FALSE,
                              cause=1,
                              ...){
  if (length(times)==0) stop("Argument 'times' has length 0")
  if (missing(type))
    type <- switch(object$model,"survival"="surv","competing.risks"="cuminc","list")
  else
    type <- switch(type,"survival"="surv","surv"="surv","incidence"="cuminc","cuminc"="cuminc","list")
  
  if (type=="surv"){
    predictSurv(object=object,
                times=times,
                newdata=newdata,
                level.chaos=level.chaos,
                mode=mode,
                bytime=bytime)
  }
  else{
    if (type=="cuminc"){
      predictCuminc(object=object,
                    times=times,
                    newdata=newdata,
                    level.chaos=level.chaos,
                    mode=mode,
                    cause=cause)
    }
    else{
      predictList(object=object,
                  times=times,
                  newdata=newdata,
                  level.chaos=level.chaos)
    }
  }
}

"predictList" <- function(object,times,newdata,level.chaos=1){
  if (missing(times)) stop("Argument times is missing.")
  NT <- length(times)
  order.times <- order(times)
  unsorted.times <- times
  times <- times[order.times]
  if (object$cens.type=="intervalCensored")
    jTimes <- object$time[2,]
  else
    jTimes <- object$time

  # no factors
  # --------------------------------------------------------------------
  if (object$covariate.type==1){
    tindex <- sindex(jump.times=jTimes,eval.times=times)
    tindex[times>object$maxtime] <- NA
    if (level.chaos==2)
      indices <- list(time=tindex[order(order.times)],strata=1)
    else
      indices <- list(time=tindex,strata=1)
    dimensions <- list(time=NT,strata=1)
    predictors <- NULL
    names.strata <- NULL
  }
  else {
    # conditional on factors
    # --------------------------------------------------------------------
    if (missing(newdata)) stop("Argument newdata is missing.")
    NX <- NROW(object$X)
    fit.X <- object$X
    ## strata.vars <- sapply(strsplit(grep("strata",names(fit.X),val=TRUE),"strata."),function(x)x[2])
    ## NN.vars <- sapply(strsplit(grep("NN",names(object$X),val=TRUE),"NN."),function(x)x[2])
    strata.vars <- object$discrete.predictors
    NN.vars <- object$continuous.predictors
    X.formula <- update(formula(object$formula),NULL~.)
    ## delete.response(terms(formula(object$formula)))
    iid <- is.null(object$clustervar)
    if (!iid){
      find.clu <- match(object$clustervar,all.vars(X.formula))
      X.formula <- drop.terms(terms(X.formula),find.clu)
    }
    if (!all(match(all.vars(X.formula),names(newdata),nomatch=FALSE)))
      stop("Arg newdata does not contain all the covariates used for fitting. \n\nrequested: ", as.character(X.formula))
    requested.X <- newdata[,all.vars(X.formula),drop=FALSE]
    NR <- NROW(requested.X)
    requested.names <- extract.name.from.special(names(requested.X))
    names(requested.X) <- requested.names
    check.vars <- match(c(strata.vars,NN.vars),requested.names,nomatch=FALSE)
    if (length(strata.vars)==0){
      requested.strata <- rep(1,NR)
      fit.strata <- rep(1,NX)
      freq.strata <- NX
    }
    else{
      # strata
      # --------------------------------------------------------------------
      requested.strata <- do.call("paste",c(requested.X[,strata.vars,drop=FALSE],sep="\r"))
      ## fit.strata <- factor(do.call("paste",c(fit.X[,paste("strata",strata.vars,sep="."),drop=FALSE],sep="\r")))
      fit.strata <- factor(do.call("paste",c(fit.X[,strata.vars,drop=FALSE],sep="\r")))
      ## changed Tue Sep 16 10:45:01 CEST 2008
      ## fit.levels <- unique(fit.strata)
      fit.levels <- as.character(unique(fit.strata))
      ## changed Tue Sep 16 10:45:01 CEST 2008
      if (!all(unique(requested.strata) %in% (fit.levels)))
        stop(paste("Not all values of newdata strata variables occur in fit:\nrequested:",
                   paste(unique(requested.strata),collapse=","),
                   "\nfitted:",
                   paste(fit.levels,collapse=",")))
      NS <- length(fit.levels)
      fit.strata <- factor(fit.strata,levels=unique(fit.strata),labels=1:NS)    
      requested.strata <- factor(requested.strata,levels=fit.levels,labels=1:NS)
      freq.strata <- cumsum(tabulate(fit.strata))
    }
    # neighborhoods
    # --------------------------------------------------------------------
    switch(length(NN.vars)+1,
           {requested.NN <- NULL
            fit.NN <- NULL
            new.order <- order(requested.strata)},
           {requested.NN <- requested.X[,NN.vars,drop=TRUE]
            fit.NN <- fit.X[,NN.vars,drop=TRUE]
            new.order <- order(requested.strata,requested.NN)
          },
           stop("Currently only one continuous covariate allowed."),
           stop("Currently only one continuous covariate allowed."))
    # findex identifies the individual strata neighborhood combination 
    # --------------------------------------------------------------------
    findex <- .C("findex",
                 index=integer(NR),
                 as.integer(as.integer(length(NN.vars)>0)),
                 as.integer(requested.strata[new.order]),
                 as.integer(freq.strata),
                 as.double(requested.NN[new.order]),
                 as.double(fit.NN),
                 as.integer(NR),
                 as.integer(NT),
                 NAOK=FALSE,
                 PACKAGE="prodlim")$index
    if (level.chaos==2) stop("Need sorted times with strata.")
    if (level.chaos==1){# do NOT sort by factors
      predictors <- requested.X
      findex <- findex[order(new.order)]
    }
    else{
      predictors <- requested.X[new.order,,drop=FALSE]
    }
    # pindex identifies the predicted probabilities
    # --------------------------------------------------------------------
    pindex <- .C("pred_index",
                 index=integer(NT*NR),
                 as.double(times),
                 as.double(jTimes),
                 as.integer(object$first.strata[findex]),
                 as.integer(object$size.strata[findex]),
                 as.integer(NR),
                 as.integer(NT),
                 NAOK=FALSE,
                 PACKAGE="prodlim")$index
    pindex[pindex==-1] <- NA
    indices <- list(time=pindex,strata=findex)
    dimensions <- list(time=NT,strata=NR)
    names.strata <- apply(do.call("cbind",lapply(names(requested.X),function(n){
      paste(n,format(requested.X[,n],digits=2),sep="=")})),1,paste,collapse=", ")
    ##     print(names.strata)
    predictors <- predictors
  }
  if (level.chaos==2) times <- unsorted.times
  else times <- times
  out <- list(times=times,
              predictors=predictors,
              indices=indices,
              dimensions=dimensions,
              names.strata=names.strata)
  out
}

predictSurv <- function(object,
                        times,
                        newdata,
                        level.chaos=1,
                        mode="list",
                        bytime=FALSE){
  p <- predict(object,
               newdata=newdata,
               level.chaos=level.chaos,
               times=times,type="list")
  NT <- p$dimensions$time
  NR <- p$dimensions$strata
  pindex <- p$indices$time
  if (object$covariate.type==1){
    psurv <- c(1,object$surv)[pindex+1]
  }
  else{
    if (bytime==FALSE){
      psurv <- split(c(1,object$surv)[pindex+1],
                     rep(1:NR,rep(NT,NR)))
      names(psurv) <- p$names.strata
    }
    else{
      psurv <- split(c(1,object$surv)[pindex+1],rep(1:NT,NR))
      names(psurv) <- paste("t",times,sep="=")
    }
  }
  if (mode=="matrix" && NR>1) {
    psurv <- do.call("rbind",psurv)
  }
  psurv
}

"predictCuminc" <- function(object,
                            times,
                            newdata,
                            level.chaos=1,
                            mode="list",
                            cause,
                            ...){
  #  if (object$model!="competing.risks") stop("This object is not a competing.risks model.")
  p <- predict(object,newdata=newdata,level.chaos=level.chaos,times=times,type="list")
  NT <- p$dimensions$time
  NR <- p$dimensions$strata
  pindex <- p$indices$time
  if (object$model=="survival")
    object$cuminc <- list("1"=1-object$surv)
  if (missing(cause))
    cause <- 1:NCOL(object$cuminc)
  if (is.character(cause))
    cause <- match(cause,names(object$cuminc))
  stopifnot(is.numeric(cause) || any(cause>NROW(object$cuminc)))
  out <- lapply(cause,function(thisCause){
    if (NR == 1){
      pcuminc <- c(0,object$cuminc[[thisCause]])[pindex+1]
      if (mode=="matrix")
        pcuminc <- matrix(pcuminc,nrow=1)
    }
    else{
      pcuminc <- split(c(0,object$cuminc[[thisCause]])[pindex+1],
                       rep(1:NR,rep(NT,NR)))
      names(pcuminc) <- p$names.strata
      if (mode=="matrix" && NR>1) {
        pcuminc <- do.call("rbind",pcuminc)
      }
    }
    pcuminc})
  if (length(cause)==1){
    out[[1]]}
  else{
    names(out) <- names(object$cuminc)[cause]
    out}
}
