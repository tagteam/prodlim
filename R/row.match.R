"row.match" <-
  function(x, table, nomatch=NA){
    if (class(table)=="matrix") table <- as.data.frame(table)
    if (is.null(dim(x))) x <- as.data.frame(matrix(x,nrow=1))
    cx <- do.call("paste",c(x[,,drop=FALSE],sep="\r"))
    ct <- do.call("paste",c(table[,,drop=FALSE],sep="\r"))
    match(cx,ct,nomatch=nomatch)
  }
