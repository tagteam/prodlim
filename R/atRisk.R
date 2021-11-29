#' Drawing numbers of subjects at-risk of experiencing an event below
#' Kaplan-Meier and Aalen-Johansen plots.
#' 
#' This function is invoked and controlled by \code{plot.prodlim}.
#' 
#' This function should not be called directly. The arguments can be specified
#' as \code{atRisk.arg} in the call to \code{plot.prodlim}.
#' 
#' @param x an object of class `prodlim' as returned by the
#' \code{prodlim} function.
#' @param newdata see \code{plot.prodlim}
#' @param times Where to compute the atrisk numbers.
#' @param line Distance of the atrisk numbers from the inner plot.
#' @param col The color of the text.
#' @param labelcol The color for the labels. Defaults to col.
#' @param interspace Distance between rows of atrisk numbers.
#' @param cex Passed on to \code{mtext} for both atrisk numbers and
#' labels.
#' @param labels Labels for the at-risk rows.
#' @param title Title for the at-risk labels
#' @param titlecol The color for the title. Defaults to 1 (black).
#' @param pos The value is passed on to the \code{mtext} argument
#' \code{at} for the labels (not the atrisk numbers).
#' @param adj Passed on to \code{mtext} for the labels (not the atriks
#' numbers).
#' @param dist If \code{line} is missing, the distance of the upper
#' most atrisk row from the inner plotting region: par()$mgp[2].
#' @param xdist Distance in x-axis direction to define the distance between the labels
#' and the numbers at-risk. Deftaults to \code{strwidth("MM",cex=cex)}.
#' @param adjust.labels If \code{TRUE} the labels are left adjusted.
#' @param show.censored If \code{TRUE} the cumulative number of subjects
#'                      lost to follow up is shown in parentheses.
#' @param unit The graphical coordinate systems unit to convert from when line2user is calling \code{grconvertX} and \code{grconvertY}.
#'             Default is \code{'npc'}
#' @param ... Further arguments that are passed to the function
#' \code{mtext}.
#' @return Nil
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{plot.prodlim}}, \code{\link{confInt}},
#' \code{\link{markTime}}
#' @keywords survival
#' @export
# {{{ header 

atRisk <- function(x,
                   newdata,
                   times,
                   line,
                   col,
                   labelcol=NULL,
                   interspace,
                   cex,
                   labels,
                   title="",
                   titlecol=NULL,
                   pos,
                   adj,
                   dist,
                   xdist,
                   adjust.labels=TRUE,
                   show.censored=FALSE,
                   unit='npc',
                   ...){

    # }}} 
    # {{{ start values
    if (missing(times)) times <- seq(0,x$maxtime,x$maxtime/10) else times <- sort(unique(times))
    if (missing(cex)) cex <- 1
    if (missing(xdist)) xdist=strwidth("MM",cex=cex)
    if (missing(pos)) pos <- NULL
    if (missing(adj)) adj <- 0
    clusterp <- length(x$clustervar[[1]])>0
    # }}}
    # {{{ find numbers at risk at given times
    if (x$model=="competing.risks"){
        px <- lifeTab(object=x,times=times,cause=getStates(x)[1],newdata=newdata,stats=NULL,intervals=FALSE,format="dt")
    }
    else if (x$model=="survival"){
        px <- lifeTab(object=x,times=times,newdata=newdata,stats=NULL,intervals=FALSE,format="dt")
    }
    if (clusterp)
        nrisk.vars <- grep("n.risk",colnames(px))
    else
        nrisk.vars <- "n.risk"
    xvars  <- data.table::key(px)
    xvars <- xvars[xvars!="cause"]
    if (length(xvars)>0){
        xdata <- px[,xvars,with=FALSE]
        xdegen <- sapply(xdata,function(x)length(unique(x))==1)
        xvars <- xvars[!xdegen]
    }
    if (length(xvars)>0){
        xdata <- xdata[,xvars,with=FALSE]
        xstrata <- apply(do.call("cbind",lapply(xvars,function(n){xdata[[n]]})),1,paste,collapse=", ")
        if (clusterp){
            number.atrisk <- c(split(px[[nrisk.vars[[1]]]],xstrata),
                               split(px[[nrisk.vars[[2]]]],xstrata))
        } else{
            number.atrisk <- split(px[[nrisk.vars[[1]]]],xstrata)
        }
    } else{
        if (clusterp){
            xdata <- data.table::data.table(c("Subjects: ","Clusters: "))
            colnames(xdata) <- ""
            number.atrisk <- c(list(px[[nrisk.vars[[1]]]]),list(px[[nrisk.vars[[2]]]]))
        }else{
            xdata <- data.table::data.table("Subjects: ")
            colnames(xdata) <- ""
            number.atrisk <- list(px[[nrisk.vars[[1]]]])
        }
    }
    nlines <- length(number.atrisk)*NCOL(number.atrisk[1])
    # }}}
    # {{{ labels 
    if (!missing(labels) && labels[[1]]!="fixme"){
        if (is.numeric(labels) || is.character(labels) || !(data.table::is.data.table(labels)))
            labels <- data.table::data.table(labels)
        if (!is.data.table(labels))
            xdata <- data.table(labels)
        else
            xdata <- labels
        n.columns <- NCOL(xdata)
    }else{
        xdata <- unique(xdata)
        if (clusterp) {
            xdata <- cbind("V1"=rep(c("Subjects: ","Clusters: "),rep(NROW(xdata),2)),data.table::rbindlist(list(xdata,xdata)))
            data.table::setnames(xdata,"V1","")
        }
        n.columns <- NCOL(xdata)
        isnum <- sapply(xdata,is.numeric)
        if (any(isnum)){
            for (j in (1:n.columns)[isnum]){
                data.table::set(xdata,j=j,value=format(xdata[[j]],digits=2))
            }
        }
    }
    # title of labels column(s)
    if (title[[1]]!=FALSE && title[[1]] !="" && length(title)>0) {
        if (length(title)==length(names(xdata)))
            names(xdata) <- title
        else{
            warning(paste0("Argument 'atRisk.title' has the wrong length.\nIt should have length ",length(names(xdata)),", i.e., one label for each column of the 'atrisk.labels'."))
        }
    }
    # distance from plot 
    if (missing(line)){
        line <- par()$mgp[1]+dist+(1:nlines)*c(1,rep(interspace,nlines-1))
    }
    # }}}
    # {{{ the actual numbers below the plot

    ## color of labels for clustered data
    if (clusterp && (length(col)==nlines/2))
        col <- rep(col,rep(2,length(col)))
    atrisk.figures.lines <- lapply(1:nlines,function(y){
        if (show.censored==FALSE){
            atrisk.figures <- as.character(number.atrisk[[y]])
        }else{
            if (x$model=="competing.risks"){
                qx <- lifeTab(object=x,times=times,cause=getStates(x)[1],newdata=newdata,stats=NULL,intervals=TRUE,format="dt")
            }
            else if (x$model=="survival"){
                qx <- lifeTab(object=x,times=times,newdata=newdata,stats=NULL,intervals=TRUE,format="dt")
            }
            if (length(xvars)>0){
                ncens <- split(qx[["n.lost"]],xstrata)
            }else{
                ncens <- list(qx[["n.lost"]])
            }
            if (show.censored=="cumulative"){
                atrisk.figures <- paste0(as.character(number.atrisk[[y]])," (",cumsum(ncens[[y]]),")")
            }else{# in interval
                atrisk.figures <- paste0(as.character(number.atrisk[[y]])," (",ncens[[y]],")")
            }
        }
    })
    lapply(1:nlines,function(y){
        atrisk.figures <- atrisk.figures.lines[[y]]
        text(labels=atrisk.figures,
             adj=c(0.5,0),
             x=times,
             y=rep(line2user(line[y],side=1,unit=unit),length(times)),
             col=rep(col[y],length(times)),
             cex=cex,
             xpd=NA,
             ...)
        if (is.null(labelcol)){
            lcol <- col[y]
        } else {
            if (is.na(labelcol[y]))
                lcol <- labelcol[1]
            else
                lcol <- labelcol[y]
        }
        if (length(pos)==0) pos <- min(times)-xdist
        if(title!=FALSE && title !="" && length(title)>0){
            column.widths <- cumsum(c(0,xdist+rev(sapply(1:n.columns,function(j){
                max(strwidth(c(names(xdata)[[j]],xdata[[j]])))}))))
        } else{
            column.widths <- cumsum(c(0,xdist+rev(sapply(1:n.columns,function(j){
                max(strwidth(xdata[[j]]))}))))
        }
        # reverse column order
        xdata <- xdata[,n.columns:1,with=FALSE]
        for (j in 1:n.columns){
            if (y==1 && title!=FALSE && title !="" && length(title)>0){
                text(labels=c(names(xdata)[[j]],as.character(xdata[y,j,with=FALSE][[1]])),
                     x=rep(pos-column.widths[[j]],2),
                     col=c(ifelse(length(titlecol)==0 || is.na(titlecol[1]),1,titlecol[1]),lcol),
                     y=c(line2user(line[1]-1,side=1,unit=unit),line2user(line[1],side=1,unit=unit)),
                     adj=c(1,0),
                     cex=cex,
                     xpd=NA,
                     ...)
            }else{
                text(labels=as.character(xdata[y,j,with=FALSE][[1]]),
                     x=pos-column.widths[[j]],
                     col=lcol,
                     y=line2user(line[y],side=1,unit=unit),
                     adj=c(1,0),
                     cex=cex,
                     xpd=NA,
                     ...)
            }
        }
    })

    # }}}
}
