getEvent <- function(object,mode="factor",column="event"){
  model <- attr(object,"model")
  if (model=="multi.state")
    stop("Dont know how to extract events from a multi.state model")
  cens.code <- attr(object,"cens.code")
  states <- attr(object,"states")
  E <- factor(object[,column],
              levels=1:(length(states)+1),
              labels=c(as.character(states),"unknown"))
  switch(mode,"character"=as.character(E),"numeric"=as.numeric(E),E)
}
