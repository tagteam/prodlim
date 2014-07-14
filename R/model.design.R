##' Extract design matrix and data specials from a model.frame
##'
##' The function \code{untangle.specials} of the survival function does a similar job.
##' @title Extract a design matrix and specials from a model.frame 
##' @param data A model.frame
##' @param xlev Passed to model.matrix
##' @param drop.intercept If TRUE drop intercept term from the design
##' matrix
##' @param max.order An error is produced if special variables are
##' involved in interaction terms of order higher than max.order.
##' @param specials.to.factor A character vector containing special variables which should be coerced into a single factor.
##' @return A list which contains
##'   - the design matrix X
##'   - the xlevels of X
##'   - separate data.frames with the values
##'     of the special variables.
##' @seealso model.frame terms model.matrix .getXlevels
##' @examples
##'
##' f <- formula(y~x+ID(z))
##' d <- data.frame(y=rnorm(5),x=factor(rbinom(5,1,0.3)+3),z=c(2,2,7,7,7))
##' ID <- function(x)x
##' t <- terms(f,special="ID",data=d)
##' m <- model.frame(t,d)
##' md <- model.design(m,specials.to.factor=TRUE)
##' md
##' 
##' library(survival)
##' data(pbc)
##' tt <- terms(Surv(time,status!=0)~factor(edema)*age+I(log(bili)>1)+strata(sex)+cluster(id),special=c("strata","cluster"),data=pbc[1:10,])
##' dd <- model.frame(tt,data=pbc[1:10,])
##' model.design(dd)
##' model.design(dd,drop.intercept=TRUE)
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
model.design <- function(data,
                         xlev,
                         drop.intercept=FALSE,
                         max.order=1,
                         specials.to.factor=TRUE){
    terms <- attr(data,"terms")
    if (drop.intercept) attr(terms, "intercept") <- 1
    design <- attr(terms,"factor")
    varnames <- rownames(design)
    terms.order <- attr(terms,"order")
    wanted.specials <- attr(terms,"specials")
    existing.specials <- wanted.specials[!sapply(wanted.specials,is.null)]
    specials <- names(existing.specials)
    names(specials) <- specials
    if (is.logical(specials.to.factor) && (specials.to.factor==TRUE)){
        specials.to.factor <-specials
    }
    if (length(specials)>0){
        special.info <- lapply(specials,function(spc){
            pos <- existing.specials[[spc]]
            ff <- apply(design[pos,,drop=FALSE],2,sum)
            terms <- seq(ff)[ff>0]
            if (any(terms.order[terms]>max.order))
                stop(paste(spc," can not be used in an interaction of order higher than ",max.order,sep=""),call.=FALSE)
            list(vars=varnames[pos],terms=terms)
        })
        special.terms <- sapply(special.info,function(x)x$terms)
        ## FIXME: there must be a better way to check if there are no
        ##        terms left after removing the specials
        tryx <- try(xterms <- drop.terms(terms,special.terms),silent=TRUE)
        if (class(tryx)[1]=="try-error")
            xterms <- NULL
        special.frames <- lapply(specials,function(sp){
            x <- special.info[[sp]]
            sp.data <- data[,x$vars,drop=FALSE]
            if (sp %in% specials.to.factor){
                if (NCOL(sp.data)>1) sp.data <- apply(sp.data,1,paste,collapse=".")
                ## force into factor
                cname <- colnames(sp.data)
                sp.data[[1]] <- data.frame(factor(sp.data[[1]]))
                colnames(sp.data[[1]]) <- cname
                sp.data[[1]]
            }else{
                sp.data
            }
        })
        if (!is.null(xterms)){
            if(missing(xlev))
                xlev <- .getXlevels(xterms,data)
            X <- model.matrix(xterms,data,xlev=xlev)
            if (drop.intercept) X <- X[,-1,drop=FALSE]
        } else {
            X <- NULL
            xlev <- NULL
        }
        c(list(X=X,xlevels=xlev),special.frames)

    }else{
        X <- model.matrix(terms,data,xlev=xlev)
        if (missing(xlev))
            xlev <- .getXlevels(terms,data)
        if (drop.intercept) X <- X[,-1,drop=FALSE]
        list(X=X,xlevels=xlev)
    }
}
