markTime <- function(x,times,nlost,pch,col,...){
  mtimeList=lapply(1:length(x),function(i){
    who=nlost[[i]]>0 & !is.na(nlost[[i]])
    mark.x=times[who]
    mark.y=x[[i]][who]
    if (length(col)<length(x)) mcol=col[1] else mcol=col[i]
    if (length(pch)<length(x)) mpch=pch[1] else mpch=pch[i]
    points(x=mark.x,y=mark.y,col=mcol,pch=mpch,...)
    ##       cbind(mark.x,mark.y)
    invisible(NULL)
  })
}


