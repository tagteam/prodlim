\name{SmartControl}
\alias{SmartControl}
\title{
  Function to facilitate the control of arguments passed to subroutines.
}
\description{
  Many R functions need to pass several arguments to 
  several different subroutines. Such arguments can are given as part of
  the three magic dots "...". The function SmartControl reads the dots
  together with a list of default values and returns for each subroutine
  a list of arguments. 
}
\usage{
SmartControl(call,
	     keys,
	     ignore,
	     defaults,
	     forced,
	     split,
	     ignore.case=TRUE,
	     replaceDefaults,
	     verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{call}{
    A list of named arguments, as for example can be obtained via \code{list(...)}. 
  }
  \item{keys}{
    A vector of names of subroutines.
}
  \item{ignore}{
    A list of names which are removed from the argument \code{call}
    before processing.
}
  \item{defaults}{
    A named list of default argument lists for the subroutines.
  }
  \item{forced}{
    A named list of forced arguments for the subroutines. 
  }
  \item{split}{
    Regular expression used for splitting keys from arguments. Default
    is \code{"\\."}.
  }
  \item{ignore.case}{
    If \code{TRUE} then all matching and splitting is not case sensitive. 
  }
  \item{replaceDefaults}{
    If \code{TRUE} default arguments are replaced by given
    arguments. Can also be a named list with entries for each subroutine.
  }
  \item{verbose}{
    If \code{TRUE} warning messages are given for arguments in
    \code{call} that are not ignored via argument \code{ignore}
    and that do not match any \code{key}.
  }
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
% \value{
% A list with arguments for all subroutines.
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
% }
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{Thomas Alexander Gerds <tag@biostat.ku.dk>}
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot.prodlim}}}
\examples{

dumPlot = function(...){
 ## set defaults
 plot.DefaultArgs=list(x=0,y=0,type="n")
 lines.DefaultArgs=list(x=1:10,lwd=3)
 ## apply smartcontrol
 x=SmartControl(call=list(...),
           defaults=list("plot"=plot.DefaultArgs, "lines"=lines.DefaultArgs),
        ignore.case=TRUE,keys=c("plot","axis2","lines"),
             forced=list("plot"=list(axes=FALSE),"axis2"=list(side=2)))
 ## call subroutines
 do.call("plot",x$plot)
 do.call("lines",x$lines)
 do.call("axis",x$axis2)
}
dumPlot(plot.ylim=c(0,5),plot.xlim=c(0,20),lines.lty=3,axis2.At=c(0,3,4))
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Graphics}
