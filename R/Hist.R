"Hist" <- function(time,
                   event,
                   entry=NULL,
                   id=NULL,
                   cens.code="0",
                   addInitialState=FALSE) {
  
  # {{{ resolving the `time' argument

  if (is.matrix(time)) time <- data.frame(time)
  if (class(time)=="list"){
    if (length(time) !=2 || length(time[[1]])!=length(time[[2]]))
      stop("Argument time has a wrong format")
    time <- data.frame(time)
  }
  if (is.data.frame(time)){
    cens.type <- "intervalCensored"
    L <- time[[1]]
    R <- time[[2]]
    N <- length(L)
    stopifnot(is.numeric(L))
    stopifnot(is.numeric(R))
    stopifnot(all(L<=R || is.na(R)))
    status <- rep(2,N)
    status[L==R] <- 1
    status[is.infinite(R) | is.na(R) | (L!=R & as.character(R)==cens.code)] <- 0
    ## the last part of the condition achieves to things:
    ##     1. for multi-state models allow transitions to a censored state
    ##     2. to ignore this, if an event occured exactly at time 0 and 0 is the cens.code
    R[status==0] <- Inf
  }
  else{
    stopifnot(is.numeric(time))
    cens.type <- "rightCensored"
    N <- length(time)
    status <- rep(1,N) ## temporary dummy
  }
  # }}}
  # {{{ resolving the `entry' argument
  if (is.null(entry))
    entry.type <- ""
  else{
    if (is.matrix(entry)) entry <- data.frame(entry)
    if (class(entry)=="list"){
      if (length(entry) !=2 || length(entry[[1]])!=length(entry[[2]]))
        stop("Argument entry has a wrong format")
      entry <- data.frame(entry)
    }
    if (is.data.frame(entry)){
      entry.type <-"intervalCensored"
      U <- entry[[1]]
      V <- entry[[2]]
      stopifnot(is.numeric(U))
      stopifnot(is.numeric(V))
      stopifnot(all(!is.na(U))|all(!is.na(V)))
    }
    else{
      stopifnot(is.numeric(entry))
      if (is.null(id))
        entry.type <- "leftTruncated"
      else
        entry.type <- "exact"
    }}
  ## check if entry < exit
  if (cens.type=="intervalCensored")
    if (entry.type=="intervalCensored")
      stopifnot(all(V<L))
    else{
      stopifnot(is.null(entry) | all(entry<=L))
    }
  else
    if (entry.type=="intervalCensored")
      stopifnot(all(V<=time))
    else
      stopifnot(is.null(entry)| all(entry<=time))
  
  # }}}
  # {{{ resolving the argument `event' 
  if (missing(event)){
    model <- "onejump"
    event <-  rep(1,N)
    warning("Argument event is missing:\nassume observations of a survival model\nand only one event per subject")
  }
  else{
    if (is.matrix(event)) event <- data.frame(event)
    if ((is.vector(event) & class(event)!="list")|| is.factor(event))
      stopifnot(length(event)==N)
    if (class(event)=="list"){
      if (length(event) !=2 || length(event[[1]])!=length(event[[2]]))
        stop("Argument event has a wrong format")
      event <- data.frame(event)
    }
    if (!is.data.frame(event)){
      if (is.null(id)){
        model <- "onejump"
        if (is.logical(event)) event <- as.numeric(event)
        status[is.na(event) | is.infinite(event) | as.character(event)==cens.code] <- 0
      }
      else{
        ## inFormat <- "longitudinal"
        stopifnot(is.numeric(id) || is.factor(id))
        model <- "multi.states"
        if (cens.type=="intervalCensored"){
          stop("Dont know the order of transitions for interval censored observations.")
        }
        else{
          if (addInitialState==TRUE){
            time <- c(rep(0,length(unique(id))),time)
            if (is.factor(event)){
              event <- factor(c(rep("initial",length(unique(id))),as.character(event)),levels=c("initial",levels(event)))
            }
            else{
              stopifnot(match("initial",unique(event),nomatch=0)==0)
              event <- c(rep("initial",length(unique(id))),event)
            }
            id <- c(unique(id),id)
            ## status <- c(rep(cens.code,length(unique(id))),status)
          }
          # 1. sort the observations by id and time
          sorted <- order(id,time)
          time <- time[sorted]
          ## status <- status[sorted] consists only of 1's
          id <- id[sorted]
          event <- event[sorted]
          # time <- time[duplicated(id)] ## remove the resp. first time
          # status <- status[duplicated(id)]
          if (length(unique(id))!=sum(time==0))
            stop("There are ",length(unique(id))," different individuals (id's), but the state at time 0 is available for ",sum(time==0)," id's.")
          initialState <- event[time==0]
          last.id <- c(diff(id) != 0, 1)
          first.id <- c(1, diff(id) != 0)
          from <- factor(event[last.id!=1])
          to <- factor(event[first.id!=1])
          id <- id[time!=0]
          time <- time[time!=0]
          # 2. get back to the original order
          ### cannot easily get back since
          ### length(time) < sorted
          ## time <- time[sorted]
          ## id <- id[sorted]
          ## event <- event[sorted]
          status <- rep(1,length(to))
          status[is.na(to) | is.infinite(to) | as.character(to)==cens.code] <- 0
        }
      }
    }
    else{
      ## inFormat <- "from2to"
      model <- "multi.states"
      from <- event[[1]]
      to <- event[[2]]
      status[is.na(to) | is.infinite(to) | as.character(to)==cens.code] <- 0
      if (length(unique(from))==1){
        model <- "onejump"
        event <- to
        if (is.logical(to)) to <- as.numeric(to)
        status[is.na(to) | is.infinite(to) | as.character(event)==cens.code] <- 0
      }
    }
  }
  ## if (all(status==0)) warning("All observations are censored")
  if (all(status==1)) cens.type <- "uncensored"
  
  if(model=="onejump"){

    # }}}
  # {{{  2-state and competing.risks models
    if (is.factor(event)){
      event <- factor(event) # drop unused levels
      states <- levels(event)
      ## states <- states[match(state.order,states)]
    }
    else{
      states <- sort(as.character(unique(event)))
    }
  
    states <- as.character(states[states!=cens.code])
  
  if (length(states)>1)
    model <- "competing.risks"
  else
      model <- "survival"
    
    if (cens.type=="intervalCensored"){
      if (model=="survival"){
        if (entry.type=="intervalCensored")
          history <- cbind(U=U,V=V,L=L,R=R,status=status)
        else
          history <- cbind(entry = entry,L=L,R=R,status=status)
      }
      else{
        if (entry.type=="intervalCensored")
          history <- cbind(U=U,V=V,L=L,R=R,status=status,event=as.integer(factor(event,levels=c(states,cens.code))))
        else
          history <- cbind(entry = entry,L=L,R=R,status=status,event=as.integer(factor(event,levels=c(states,cens.code))))
      }
    }
    else{
      if (model=="survival"){
        if (entry.type=="intervalCensored")
          history <- cbind(U=U,V=V,time=time,status=status)
        else
          history <- cbind(entry = entry,time=time,status=status)
      }
      else{
        if (entry.type=="intervalCensored")
          history <- cbind(U=U,V=V,time=time,status=status,event=as.integer(factor(event,levels=c(states,cens.code))))
        else
          history <- cbind(entry = entry,time=time,status=status,event=as.integer(factor(event,levels=c(states,cens.code))))
      }
    }
  } else{
    # }}}
    # {{{  multi.state models

    if (any(as.character(from)==as.character(to))) stop("Data contain transitions from state x to state x")

    eventISfactor <- as.numeric(is.factor(from)) + as.numeric(is.factor(to))

    if (eventISfactor==1) stop("Components of event have different classes")
    
    if (eventISfactor==2)
      states <- unique(c(levels(from),levels(to)))
    else
      states <- as.character(unique(c(from,to)))
    
    states <- as.character(states[states!=cens.code])
    ## states <- states[match(state.order,states)]
    if (cens.code %in% levels(from)) stop(paste("The Cens.code",cens.code," identifies censored data, but is found amoung the `from' state of some transitions")

    if (cens.type=="intervalCensored"){
      if (entry.type=="intervalCensored")
        history <- cbind(U=U,V=V,L=L,R=R,status=status,from=as.integer(factor(from,levels=c(states,cens.code))),to=as.integer(factor(to,levels=c(states,cens.code))))
      else{
        history <- cbind(entry = entry,L=L,R=R,status=status,from=as.integer(factor(from,levels=c(states,cens.code))),to=as.integer(factor(to,levels=c(states,cens.code))))
      }
    }
    else{
      if (entry.type=="intervalCensored")
        history <- cbind(U=U,V=V,time=time,status=status,from=as.integer(factor(from,levels=c(states,cens.code))),to=as.integer(factor(to,levels=c(states,cens.code))))
      else{
        history <- cbind(entry = entry,time=time,status=status,from=as.integer(factor(from,levels=c(states,cens.code))),to=as.integer(factor(to,levels=c(states,cens.code))))
      }
    }
  }

  # }}}
  # {{{ add id

  if (!is.null(id)) history <- cbind(history,id)
  
  # }}}
  # {{{ class and attributes
  rownames(history) <- NULL
  class(history) <- c("Hist")
  attr(history,"states") <- states
  attr(history,"cens.type") <- cens.type
  attr(history,"cens.code") <- as.character(cens.code)
  attr(history,"model") <- model
  ## print(entry.type)
  attr(history,"entry.type") <- entry.type
  history
  # }}}
}

subset.Hist <- function(x,subset,select,drop){
  if (missing(select)){
    xx <- x
    class(xx) <- "matrix"
    xx <- subset(xx,subset=subset,drop=drop)
    attr(xx,"class") <- attr(x,"class")
    attr(xx,"states") <- attr(x,"states")
    attr(xx,"model") <- attr(x,"model")
    attr(xx,"cens.type") <- attr(x,"cens.type")
    attr(xx,"cens.code") <- attr(x,"cens.code")
    attr(xx,"entry.type") <- attr(x,"entry.type")
    xx
  }
  else{
    class(x) <- "matrix"
    NextMethod("subset")
  }
}

"[.Hist" <- function(x,i,j,drop=FALSE){
  if (missing(j)){
    xx <- x
    class(xx) <- "matrix"
    xx <- xx[i,,drop=drop]
    class(xx) <- "Hist"    
    attr(xx,"class") <- attr(x,"class")
    attr(xx,"states") <- attr(x,"states")
    attr(xx,"model") <- attr(x,"model")
    attr(xx,"cens.type") <- attr(x,"cens.type")
    attr(xx,"cens.code") <- attr(x,"cens.code")
    attr(xx,"entry.type") <- attr(x,"entry.type")
    xx
  }
  else{
    class(x) <- "matrix"
    ## x[i,j,drop=drop]
    NextMethod("[")
  }
}

# does not work
# as.data.frame.Hist <- function(x,...){
#  class(x) <- "matrix"
#  as.data.frame(x)
# }
  

is.na.Hist <- function(x) {
  as.vector( (1* is.na(unclass(x)))%*% rep(1, ncol(x)) >0)
}

str.Hist <- function(x){
  class(x) <- "matrix"
  str(x)
}

