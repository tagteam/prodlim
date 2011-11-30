backGround <- function(xlim,
                       ylim,
                       bg="white",
                       fg="gray77",
                       horizontal=NULL,
                       vertical=NULL,
                       border="black"){
  U <- par("usr")
  if (missing(xlim))
    xlim <- c(U[1],U[2])
  if (missing(ylim))
    ylim <- c(U[3],U[4])
  # background
  if (!is.null(bg)){
    rect(U[1],U[3],U[2],U[4],col=bg[1], border=border)
    if (length(bg)>1){
      ybot <- sort(unique(c(ylim[1],horizontal,ylim[2])))
      NR <- length(ybot)
      bcol <- rep(bg,length.out=NR)
      nix <- sapply(1:(NR-1),function(r){
        ## for (r in 1:(NR-1)){
        ## rect(xleft=xlim[1],xright=xlim[2],ybottom=ybot[r],ytop=ybot[r+1],col=bcol[r],border=FALSE)
        ## polygon(x=c(xlim[1],xlim[1],xlim[2],xlim[2],xlim[1]),
        polygon(x=c(U[1],U[1],U[2],U[2],U[1]),
                y=c(ybot[r],ybot[r+1],ybot[r+1],ybot[r],ybot[r]),
                col=bcol[r],
                border=FALSE)
        ## do NOT specify: density=100 as this slows this down!
      })
    }
  }
  # grid 
  if (length(fg)>0){
    if (length(vertical)>0)
      abline(v=vertical,col=fg)
    if (length(horizontal)>0)
      abline(h=horizontal,col=fg)
  }
}
