SimCompRisk <- function(N,
                        NC=2,
                        cens,
                        cova,
                        verbose=1,
                        ...){
  # {{{ process arguments  
  cumincNames <- paste("cuminc",1:NC,sep="")
  default.cuminc.args <- lapply(cumincNames,function(n){
    list(surv.model="Cox-Weibull",shape=1,baseline=1,link="exp",coef=c(1,-1),transform=NULL)
  })
  names(default.cuminc.args) <- cumincNames
  default.cens.args <- list(model="Cox-exponential",baseline=1/100,link="exp",max=NULL,type="right",coef=NULL,transform=NULL)
  default.cova.list <- list(X1=list("rnorm",mean=0,sd=2),X2=list("rbinom",size=1,prob=.5))
  smartA <- SmartControl(call=match.call(),keys=c(cumincNames,"cens","cova"),ignore=c("N",cumincNames,"cens","cova","verbose"),defaults=c(default.cuminc.args,list("cens"=default.cens.args,"cova"=default.cova.list)),ignore.case=FALSE,replaceDefaults=c(sapply(cumincNames,function(x)FALSE),c("cens"=FALSE,"cova"=length(grep("cova\\.*",names(match.call())))>0)),verbose=TRUE)
  cuminc1 <- smartA[[grep("^cuminc1",names(smartA))]]
  cuminc2 <- smartA[[grep("^cuminc2",names(smartA))]]
  cens <- smartA$cens
  # }}}
  # {{{ resolving covariates
  X.matrix <- resolveX(object=smartA$cova,N=N)
  NP <- NCOL(X.matrix)
  cova <- smartA$cova
  # }}}
  # {{{ sampling the cause and latent times
  lp1 <- resolveLinPred(X=X.matrix,coef=cuminc1$coef)
  lp2 <- resolveLinPred(X=X.matrix,coef=cuminc2$coef)
  elp <- (exp(lp1)/(1+exp(lp1)))
  cause <- do.call("rbinom",list(n=N,size=1,prob=elp))+1
  T1 <- SimSurvInternal(N=sum(cause==1),model=cuminc1$surv.model,link=cuminc1$link,baseline=cuminc1$baseline,linpred=lp1[cause==1],shape=cuminc1$shape)
  T2 <- SimSurvInternal(N=N-sum(cause==1),model=cuminc2$surv.model,link=cuminc2$link,baseline=cuminc2$baseline,linpred=lp2[cause==2],shape=cuminc2$shape)
  # }}}
  # {{{ censoring 
  if (length(cens$type)>0)
    censType <- match(cens$type,c("interval","right"),nomatch=0)
  censnotwanted <- ((!missing(cens) && is.logical(cens) && cens==FALSE) || censType==0)
  if (censnotwanted==TRUE)
    cens.time <- rep(Inf,N)
  else{ # special links between right censoring and covariates
    if (length(cens$transform)>0){
      censSpecials <- TRUE
      cens.X <- transformX(X=X.matrix,transform=cens$transform,transName="f")
    }
    else{
      censSpecials <- FALSE
      cens.X <- X.matrix
    }
    linpred.cens <- resolveLinPred(X=cens.X,coef=cens$coef,verbose=verbose)
    cens.time <- SimSurvInternal(N=N,model=cens$model,link=cens$link,baseline=cens$baseline,linpred=linpred.cens,shape=cens$shape)
  }
  if (is.numeric(cens$max)) cens.time <- pmin(cens.time, cens$max)
  # }}}
  # {{{ event time and censoring status

  surv.time <- numeric(N)
  surv.time[cause==1] <- T1
  surv.time[cause==2] <- T2
  status <- as.numeric(surv.time <= cens.time)

  # }}}
  # {{{ return data frame
  out <- data.frame(cbind(time = pmin(surv.time, cens.time),status = status,cause=cause*(status>0),ucens.cause=cause,ucens.time=surv.time,cens.time=cens.time))
  if (!is.null(X.matrix)) out <- cbind(out,X.matrix)
  out
  # }}}
}

## lp <- c(X %*% beta)
# lp <- c(X %*% beta)
## cause1 <- rbinom(n,1,p1x)
## time <- (-1)*log(1-((1-(1-u*p1x)^(exp(-lRR)))/p1))
## if (trans=="logistic") # what model is this
##   lRR <- c(X %*% beta); RR <- exp(lRR);
##   cause1 <- rbinom(n,1,p1x)
##   time <- (-1)*log(1-exp(log(u*p1x/(1-u*p1x))-lRR)*(1-p1)/p1)
 
