#' Predicting event probabilities from product limit estimates
#' 
#' Evaluation of estimated survival or event probabilities at given times and
#' covariate constellations.
#' 
#' Predicted (survival) probabilities are returned that can be plotted,
#' summarized and used for inverse of probability of censoring weighting.
#' 
#' @aliases predict.prodlim predictSurv predictAbsrisk predictCuminc
#' @param object A fitted object of class "prodlim".
#' @param times Vector of times at which to return the estimated probabilities (survival or absolute event risks).
#' @param newdata A data frame with the same variable names as those that
#' appear on the right hand side of the 'prodlim' formula.  If there are
#' covariates this argument is required.
#' @param level.chaos Integer specifying the sorting of the output: `0' sort by
#' time and newdata; `1' only by time; `2' no sorting at all
#' @param type Choice between "surv","risk","cuminc","list":
#' 
#' "surv": predict survival probabilities only survival models
#' 
#' "risk"/"cuminc": predict absolute risk, i.e., cumulative incidence function.
#' 
#' "list": find the indices corresponding to times and newdata. See value.
#' 
#' Defaults to "surv" for two-state models and to "risk" for competing risk
#' models.
#' @param mode Only for \code{type=="surv"} and \code{type=="risk"}. Can
#' either be "list" or "matrix". For "matrix" the predicted probabilities will
#' be returned in matrix form.
#' @param bytime Logical. If TRUE and \code{mode=="matrix"} the matrix with
#' predicted probabilities will have a column for each time and a row for each
#' newdata. Only when \code{object$covariate.type>1} and more than one time is
#' given.
#' @param cause Character (other classes are converted with \code{as.character}).
#' The cause for predicting the absolute risk of an event, i.e., the cause-specific cumulative
#' incidence function, in competing risk models. At any time after time zero this is the absolute risk of 
#' an event of type \code{cause} to occur between time zero and \code{times} . 
#' @param \dots Only for compatibility reasons.
#' @return \code{type=="surv"} A list or a matrix with survival probabilities
#' for all times and all newdata.
#' 
#' \code{type=="risk"} or \code{type=="cuminc"} A list or a matrix with cumulative incidences for all
#' times and all newdata.
#' 
#' \code{type=="list"} A list with the following components:
#' 
#' \item{times}{The argument \code{times} carried forward}
#' 
#' \item{predictors}{The relevant part of the argument \code{newdata}.}
#' \item{indices}{ A list with the following components
#' 
#' \code{time}: Where to find values corresponding to the requested times
#' \code{strata}: Where to find values corresponding to the values of the
#' variables in newdata.  Together time and strata show where to find the
#' predicted probabilities.  } \item{dimensions}{ a list with the following
#' components: \code{time} : The length of \code{times} \code{strata} : The
#' number of rows in \code{newdata} \code{names.strata} : Labels for the
#' covariate values.  }
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{predictSurvIndividual}}
#' @keywords survival
#' @examples
#' 
#' 
#' dat <- SimSurv(400)
#' fit <- prodlim(Hist(time,status)~1,data=dat)
#' 
#' ## predict the survival probs at selected times 
#' predict(fit,times=c(3,5,10))
#' 
#' ## NA is returned when the time point is beyond the
#' ## range of definition of the Kaplan-Meier estimator:
#' predict(fit,times=c(-1,0,10,100,1000,10000))
#' 
#' ## when there are strata, newdata is required 
#' ## or neighborhoods (i.e. overlapping strata)
#' mfit <- prodlim(Hist(time,status)~X1+X2,data=dat)
#' predict(mfit,times=c(-1,0,10,100,1000,10000),newdata=dat[18:21,])
#' 
#' ## this can be requested in matrix form
#' predict(mfit,times=c(-1,0,10,100,1000,10000),newdata=dat[18:21,],mode="matrix")
#' 
#' ## and even transposed
#' predict(mfit,times=c(-1,0,10,100,1000,10000),newdata=dat[18:21,],mode="matrix",bytime=TRUE)
#' 
#' @export 
"predict.prodlim" <- function(object,
                              times,
                              newdata,
                              level.chaos=1,
                              type=c("surv","risk","cuminc","list"),
                              mode="list",
                              bytime=FALSE,
                              cause,
                              ...){
    if (!missing(newdata)&&!is.null(newdata)&&inherits(newdata,"data.table"))
        data.table::setDF(newdata)
    if (missing(cause)){
        cause <- attr(object$model.response,"states")
    } else {
        cause <- checkCauses(cause,object)
    }
    if (length(times)==0) stop("Argument 'times' has length 0")
    if (missing(type))
        type <- switch(object$model,"survival"="surv","competing.risks"="risk","list")
    else
        type <- switch(type,"survival"="surv","surv"="surv","risk"="risk","cuminc"="risk","list")
  
  if (type=="surv"){
    predictSurv(object=object,
                times=times,
                newdata=newdata,
                level.chaos=level.chaos,
                mode=mode,
                bytime=bytime)
  }
  else{
    if (type=="risk"){
      predictAbsrisk(object=object,
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
      requested.X <- NULL
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
          stop("Arg newdata does not contain all the covariates used for fitting. \n\nfitted variables: ", paste(all.vars(X.formula),collapse=", "),"\nnewdata contains:",ifelse(length(names(newdata))==0," nothing",names(newdata)))
      requested.X <- newdata[,all.vars(X.formula),drop=FALSE]
      NR <- NROW(requested.X)
      requested.names <- extract.name.from.special(names(requested.X))
      names(requested.X) <- requested.names
      check.vars <- match(c(strata.vars,NN.vars),requested.names,nomatch=FALSE)
      if (length(strata.vars)==0){
          requested.strata <- rep(1,NR)
          fit.strata <- rep(1,NX)
          freq.strata <- NX
      } else{
          # strata
          # --------------------------------------------------------------------
          ## changed 09 Dec 2014 (16:44) -->
          ## requested.strata <- do.call("paste",c(requested.X[,strata.vars,drop=FALSE],sep="\r"))
          fit.strata <- interaction(fit.X[,strata.vars,drop=FALSE],sep=":",drop=TRUE)
          requested.strata <- interaction(requested.X[,strata.vars,drop=FALSE],sep=":",drop=TRUE)
          fit.levels <- as.character(unique(fit.strata))
          ## <-- changed 09 Dec 2014 (16:44)          
          ## before version 1.5.1
          ## fit.strata <- factor(do.call("paste",c(fit.X[,strata.vars,drop=FALSE],sep="\r")))
          ## fit.levels <- unique(fit.strata)
          if (!all(unique(requested.strata) %in% (fit.levels))){
              stop(paste("Not all values of newdata strata variables occur in fit:\nrequested:",
                         paste(unique(requested.strata),collapse=","),
                         "\nfitted:",
                         paste(fit.levels,collapse=",")))
          }
          NS <- length(fit.levels)
          ## fit.strata <- factor(fit.strata,levels=unique(fit.strata),labels=1:NS)
          fit.strata <- factor(fit.strata,levels=levels(fit.strata),labels=1:NS)
          requested.strata <- factor(requested.strata,levels=fit.levels,labels=1:NS)
          freq.strata <- cumsum(tabulate(fit.strata))
      }
      # neighborhoods
      # --------------------------------------------------------------------
      switch(length(NN.vars)+1,
      {requested.NN <- NULL
          fit.NN <- NULL
          new.order <- order(requested.strata)},
      {
          requested.NN <- requested.X[,NN.vars,drop=TRUE]
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
      if (level.chaos==2) stop("Need to sort the times if there are strata.")
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
      ## bug fix (10 Oct 2013 (10:08)):
      ##        order of names needs to
      ##        obey level.chaos
      names.strata <- apply(do.call("cbind",lapply(names(requested.X),function(n){
          if(is.numeric(requested.X[,n]))
              paste(n,format(requested.X[,n],digits=2),sep="=")
          else 
              paste(n,requested.X[,n],sep="=")})),1,paste,collapse=", ")
      if (level.chaos==0) {names.strata <- names.strata[new.order]}
  }
  if (level.chaos==2) times <- unsorted.times
  else times <- times
  out <- list(times=times,
              predictors=predictors,
              indices=indices,
              dimensions=dimensions,
              strata=requested.X,
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
        ## psurv <- cbind(p$strata,do.call("rbind",psurv))
        psurv <- do.call("rbind",psurv)
    } 
    psurv
}

"predictAbsrisk" <- function(object,
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
    if (object$model=="survival"){
        object$cuminc <- list("1"=1-object$surv)
        cause <- 1
    }
    if (object$model=="competing.risks"){
        if (missing(cause))
            cause <- attributes(object$model.response)$states
        else
            cause <- checkCauses(cause,object)
    }
    out <- lapply(cause,function(thisCause){
        if (NR == 1){
            prisk <- c(0,object$cuminc[[thisCause]])[pindex+1]
            if (mode=="matrix")
                prisk <- matrix(prisk,nrow=1)
        }
        else{
            prisk <- split(c(0,object$cuminc[[thisCause]])[pindex+1],
                             rep(1:NR,rep(NT,NR)))
            names(prisk) <- p$names.strata
            if (mode=="matrix" && NR>1) {
                prisk <- do.call("rbind",prisk)
            }
        }
        prisk})
    if (length(cause)==1){
        out[[1]]
    } else{
        names(out) <- cause
        out
    }
}
