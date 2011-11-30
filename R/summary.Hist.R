summary.Hist <- function(object, verbose=TRUE,...){
  D <- object[,"status",drop=TRUE]
  states <- attr(object,"states")
  cens.code <- attr(object,"cens.code")
  # {{{ resolving events and model states 
  model <- attr(object,"model")
  model.string <- paste("response of a", model,"model")
 if (model=="multi.states"){
    from <- object[,"from"]
    to <- object[,"to"]
    code.from <- getEvent(object,mode="factor",column="from")
    code.to <- getEvent(object,mode="factor",column="to")
    state.types <- factor(as.numeric(match(states,unique(code.from),nomatch=0)!=0) + 2*as.numeric(match(states,unique(code.to),nomatch=0)!=0),levels=c(1,2,3))
    names(state.types) <- states
    levels(state.types) <- c("initial","absorbing","transient")
    state.types <- table(state.types)
  }
  else{
    from <- rep("initial",NROW(object))
    code.to <- getEvent(object,mode="factor",column=ifelse(model=="survival","status","event"))
    code.from <- factor(from)
    state.types <- c(1,length(states))
    names(state.types) <- c("initial","absorbing")
  }
  # }}}
  # {{{ transition frame
  ##   trans.frame <- unique(data.frame(from=code.from,to=code.to),MARGIN=1)
  trans.frame <- data.frame(from=code.from,to=code.to)
  Transitions <- apply(cbind(as.character(code.from),as.character(code.to)),1,paste,collapse=" -> ")
  obnoxious.factor.levels <- unique(Transitions)
  Transitions <- factor(Transitions,obnoxious.factor.levels)
  transitions <- table(Transitions)
  summary.out <- list(states=state.types,transitions=transitions,trans.frame=trans.frame)
  if (verbose==TRUE){
    state.table <- as.matrix(transitions)
    colnames(state.table) <- c("Freq")
  }
  # }}}
  # {{{ resolving the censoring mechanism
  if (verbose==TRUE){
    ## event time
    cens.type <- attr(object,"cens.type")
    ## cens.string <- capitalize(cens.type)
    cens.string <- cens.type
    Observations <- switch(cens.type,"intervalCensored"=factor(D,levels=c(1,2,0),labels=c("exact.time","intervalCensored","rightCensored")),"rightCensored"=factor(D,levels=c(1,0),labels=c("event","right.censored")),"uncensored"=factor(D,labels=c("event")))
    Freq <- table(Observations)
    ## entry time
    entry.type <- attr(object,"entry.type")
    if (entry.type!="")
      entry.string <- paste("with ",entry.type," entry time")
    else
      entry.string <- ""
    cat(paste("\n",cens.string," ",model.string," ",entry.string,"\n",sep=""))
    cat("\nNo.Observations:",NROW(object),"\n\nPattern:\n")
    switch(model,"survival"={
      prmatrix(cbind(names(Freq),Freq),
               quote=FALSE,
               rowlab=rep("",NROW(Freq)))},
           "competing.risks"={
             events <- getEvent(object)
             prout <- table("Cause"=events,as.character(Observations))
             print(prout)
           },
           "multi.states"={
             x=table(Transitions,Observations)
             aaa=sapply(strsplit(rownames(x)," -> "),function(x)x[1])
             bbb=sapply(strsplit(rownames(x)," -> "),function(x)x[1])
             print(x[order(aaa,bbb),,drop=FALSE])
           })
  }
  # }}}
  invisible(summary.out)
}

## capitalize <- function(x) {
  ## s <- strsplit(x, " ")[[1]]
  ## paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
## }
