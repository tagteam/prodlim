transformX <- function(X,transform,transName="f"){
  specials <- length(transform)>0
  if (specials){
    found.link <- match(names(transform),colnames(X),nomatch=FALSE)
    which.link <- match(colnames(X),names(transform),nomatch=FALSE)>0
    if (any(!(found.link)))
      stop("Argument transform must be a named list matching the covariate names")
    XX <- do.call("cbind",lapply(colnames(X),function(p){
      if (match(p,names(transform),nomatch=FALSE)){
        sapply(X[,p],transform[[p]])
      }
      else
        X[,p]
    }))
    colnames(XX)[which.link] <- paste(transName,colnames(XX)[which.link],sep=".")
    XX
  }
  else
    NULL
}
