##' Extract from a vector of character strings the names of special functions and auxiliary arguments
##' 
##' The names of the special functions may not be nested, i.e., c(treat, treatment) will not work.
##' @title Parse special terms 
##' @param x Vector of character strings, usually the column names of
##' the design matrix obtained with \code{\link{model.design}}.
##' @param specials A vector with character strings providing the
##' names of the special arguments.
##' @param specialArgumentNames A named list with one element for each
##' special of argument
##' @return A named list of parsed arguments. The names of the list are the variable names.
##' @seealso model.design
##' @examples
##' parseSpecialNames("treat(Z)",specials="treat")
##' parseSpecialNames("treat(Z,u=2)",specials="treat",specialArgumentNames=list("treat"="u"))
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
##' @export
parseSpecialNames <- function(x,specials,specialArgumentNames){
    if (missing(specialArgumentNames)) {
        specialArgumentNames <- lapply(specials,function(x)NULL)
    }
    else {
        if (!(all(match(names(specialArgumentNames),specials,nomatch=0))))
            stop("Mispecified argument specialArgumentNames")
    }
    split.regexp <- paste("(",paste(specials,collapse="|"),")\\(|)$",sep="")
    found.specials <- grep(specials,x,value=FALSE)
    spc.order <- specials[found.specials]
    listTerms <- strsplit(x,split.regexp)
    ## if length is 1 then term is unspecial
    ## isSpecial <- sapply(listTerms,length)
    # check for further arguments
    listTermsWithArguments <- unlist(lapply(listTerms,function(x){
        if (length(x)<2) NULL
        else
            strsplit(x[[2]],"[ ]*,[ ]*")}), recursive=FALSE)
    specialVarnames <- sapply(listTermsWithArguments,function(x){x[[1]]})
    if (length(problem <- grep("=",specialVarnames,value=TRUE))>0)
        stop(paste("Problematic variable name '",problem,"'. Variable names used in specials may not contain '='.",sep=""))
    specialArguments <- lapply(listTermsWithArguments,function(x){
        if (length(x)==1) NULL else x[2:length(x)]
    })
    names(specialArguments) <- specialVarnames
    if (length(specialArguments)>0){
        specialArgumentList <- lapply(1:length(specialArguments),function(i){
            args <- specialArguments[[i]]
            spc <- spc.order[[i]]
            if (!is.null(args)){
                fullvalue <- strsplit(args,"=")
                fullvalue <- lapply(fullvalue,function(x){ ## remove whitespace
                    gsub(" ","",x)
                })
                givennames <- sapply(fullvalue,function(x){
                    if (length(x)==1)
                        ""
                    else
                        x[[1]]
                })
                values <- lapply(fullvalue,function(x){
                    if (length(x)==1)
                        x[[1]]
                    else
                        x[[2]]
                })
                if (is.null(specialArgumentNames[[spc]]))
                    wantednames <- paste("Arg",1:length(args),sep=".")
                else{
                    wantednames <- specialArgumentNames[[spc]]
                    if(length(wantednames)<length(args)) stop("Too many arguments for special function ",spc)
                }
                realnames <- givennames[givennames!=""]
                thismatch <- match(realnames,wantednames,nomatch=0)
                if (length(realnames)>0)
                    if (!all(thismatch))
                        stop("Argument(s) ",
                             realnames,
                             " is not an argument of  ",
                             spc,
                             ". Valid argument(s): ",
                             paste(wantednames,collapse=", "))
                names(values) <- givennames
                nadd <- length(wantednames)-length(values)
                if (nadd>0){
                    values <- c(values,rep(NA,nadd))
                }
                thatmatch <- match(wantednames,names(values),nomatch=0)
                names(values)[names(values)==""] <- wantednames[thatmatch==0]
                values
            }
            else NULL
        })
        names(specialArgumentList) <- names(specialArguments)
        specialArgumentList
    }else{NULL}
}
