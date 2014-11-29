##' Extract design matrix and data specials from a model.frame
##'
##' The function \code{untangle.specials} of the survival function does a similar job.
##' @title Extract a design matrix and specials from a model.frame 
##' @param data A model.frame
##' @param xlev Passed to model.matrix
##' @param dropIntercept If TRUE drop intercept term from the design
##' matrix
##' @param maxOrder An error is produced if special variables are
##' involved in interaction terms of order higher than max.order.
##' @param unspecialsDesign A logical value: if \code{TRUE} apply \code{\link{model.matrix}}
##' to unspecial covariates. If \code{FALSE} extract unspecial covariates from data.
##' @param specialsFactor A character vector containing special
##' variables which should be coerced into a single factor. If
##' \code{TRUE} all specials are treated in this way.
##' @param specialsDesign A character vector containing special
##' variables which should be transformed into a design matrix via
##' \code{\link{model.matrix}}.  If \code{TRUE} all specials are
##' treated in this way.
##' @param stripSpecialNames If TRUE strip the special from variable
##' name, i.e., use X instead of strata(X).
##' @return A list which contains
##'   - the design matrix with the levels of the variables stored in attribute 'levels' 
##'   - separate data.frames which contain the values of the special variables.
##' @seealso \code{\link{EventHistory.frame}} model.frame terms model.matrix .getXlevels  
##' @examples
##'
##' f <- formula(y~x+ID(z))
##' set.seed(8)
##' d <- data.frame(y=rnorm(5),x=factor(c("a","b","b","a","c")),z=c(2,2,7,7,7))
##' ID <- function(x)x
##' t <- terms(f,special="ID",data=d)
##' m <- model.frame(t,d)
##' md <- model.design(m,specialsFactor=TRUE)
##' md
##' md <- model.design(m,specialsFactor=TRUE,unspecialsDesign=FALSE)
##' md
##'
##' # special function with argument
##' treat <- function(x,...) x
##' f2 <- formula(y~x+treat(z,arg=2)+treat(u,arg=-1))
##' set.seed(8)
##' d <- data.frame(y=rnorm(5),u=1:5,x=factor(c("a","b","b","a","c")),z=c(2,2,7,7,7))
##' t2 <- terms(f2,special="treat",data=d)
##' m2 <- model.frame(t2,d)
##' md2 <- model.design(m2,specialsFactor=TRUE)
##' 
##' library(survival)
##' data(pbc)
##' tt <- terms(Surv(time,status!=0)~factor(edema)*age+strata(I(log(bili)>1))+strata(sex)+cluster(id),
##' special=c("strata","cluster"),data=pbc[1:10,])
##' dd <- model.frame(tt,data=pbc[1:10,])
##' model.design(dd)
##' model.design(dd,dropIntercept=TRUE)
##' 
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
model.design <- function(data,
                         xlev,
                         dropIntercept=FALSE,
                         maxOrder=1,
                         unspecialsDesign=TRUE,
                         specialsFactor=TRUE,
                         specialsDesign=FALSE,
                         stripSpecialNames=TRUE){
    # {{{ analyse the terms
    terms <- attr(data,"terms")
    if (dropIntercept) attr(terms, "intercept") <- 1
    design <- attr(terms,"factor")
    varnames <- rownames(design)
    termsOrder <- attr(terms,"order")
    wantedSpecials <- attr(terms,"specials")
    existingSpecials <- wantedSpecials[sapply(wantedSpecials,length)!=0]
    specials <- names(existingSpecials)
    names(specials) <- specials
    if (is.logical(specialsDesign) && (specialsDesign==TRUE)){
        specialsDesign <- specials
    }
    if (is.logical(specialsFactor) && (specialsFactor==TRUE)){
        specialsFactor <- specials
    }
    # }}}
    if (length(specials)>0){
        # {{{ extract information about specials
        specialInfo <- lapply(specials,function(spc){
            pos <- existingSpecials[[spc]]
            ff <- apply(design[pos,,drop=FALSE],2,sum)
            terms <- seq(ff)[ff>0]
            if (any(termsOrder[terms]>maxOrder))
                stop(paste(spc,
                           " can not be used in an interaction of order higher than ",
                           maxOrder,
                           sep=""),call.=FALSE)
            ## extract additional arguments from term.lables 
            spc.vnames <- varnames[pos]
            list(vars=varnames[pos],terms=as.vector(terms))
        })
        specialTerms <- unlist(lapply(specialInfo,function(x)x$terms))
        termLabels <- attr(terms,"term.labels")
        ## only specials
        if (length(termLabels) == length(specialTerms))
            unspecialTerms <- NULL
        else
            unspecialTerms <- drop.terms(terms,specialTerms)
        # }}}
        # {{{ loop over specials
        specialFrames <- lapply(specials,function(sp){
            Info <- specialInfo[[sp]]
            spTerms <- terms[Info$terms]
            spLevels <- .getXlevels(spTerms,data)
            if (stripSpecialNames==TRUE)
                names(spLevels) <- sub("[^(]*\\((.*)\\)","\\1",names(spLevels))
            if (sp %in% specialsDesign){
                spMatrix <- model.matrix(spTerms,data=data,xlev=spLevels)[,-1,drop=FALSE]
                if (stripSpecialNames==TRUE)
                    colnames(spMatrix) <- sub("[^(]*\\((.*)\\)","\\1",colnames(spMatrix))
                attr(spMatrix,"levels") <- spLevels
                spMatrix
            }else{
                spData <- data[,Info$vars,drop=FALSE]
                cnames <- colnames(spData)
                if (stripSpecialNames==TRUE)
                    cnames <- sub("[^(]*\\((.*)\\)","\\1",cnames)
                colnames(spData) <- cnames
                if (sp %in% specialsFactor){
                    ## force into a single factor
                    if (NCOL(spData)>1) {
                        spData <- data.frame(apply(spData,1,paste,collapse=", "))
                        names(spData) <- paste(cnames,collapse=", ")
                    }
                }
                attr(spData,"levels") <- spLevels
                spData
            }
        })
        # }}}
        # {{{ unspecials
        if (!is.null(unspecialTerms)){
            if(missing(xlev))
                xlev <- .getXlevels(unspecialTerms,data)
            if (unspecialsDesign==TRUE){
                X <- model.matrix(unspecialTerms,data,xlev=xlev)
                if (dropIntercept) X <- X[,-1,drop=FALSE]
            }
            else{
                X <- model.frame(unspecialTerms,data)
            }
        } else {
            X <- NULL
            xlev <- NULL
        }
        attr(X,"levels") <- xlev
        c(list(design=X),specialFrames)
        # }}}
    }else{
        # {{{ no specials
        if (unspecialsDesign==TRUE){
            X <- model.matrix(terms,data,xlev=xlev)
            if (dropIntercept) X <- X[,-1,drop=FALSE]
        }else {
            X <- model.frame(delete.response(terms),data)
        }
        if (missing(xlev))
            xlev <- .getXlevels(terms,data)
        attr(X,"levels") <- xlev
        list(design=X)
        # }}}
    }
}
