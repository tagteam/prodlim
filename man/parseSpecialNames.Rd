\name{parseSpecialNames}
\alias{parseSpecialNames}
\title{Parse special terms}
\usage{
parseSpecialNames(x, special, arguments)
}
\arguments{
  \item{x}{Vector of character strings.}

  \item{special}{A character string: the name of the
  special argument.}

  \item{arguments}{A vector which contains the arguments of
  the special function}
}
\value{
A named list of parsed arguments. The names of the list are
the special variable names, the elements are lists of
arguments.
}
\description{
Extract from a vector of character strings the names of
special functions and auxiliary arguments
}
\details{
Signals an error if an element has more arguments than
specified by argument arguments.
}
\examples{
## ignore arguments
parseSpecialNames("treat(Z)",special="treat")
## set default to 0
parseSpecialNames(c("log(Z)","a","log(B)"),special="log",arguments=list("base"=0))
## set default to 0
parseSpecialNames(c("log(Z,3)","a","log(B,base=1)"),special="log",arguments=list("base"=0))
## different combinations of order and names
parseSpecialNames(c("log(Z,3)","a","log(B,1)"),
                  special="log",
                  arguments=list("base"=0))
parseSpecialNames(c("log(Z,1,3)","a","log(B,u=3)"),
                  special="log",
                  arguments=list("base"=0,"u"=1))
parseSpecialNames(c("log(Z,u=1,base=3)","a","log(B,u=3)"),
                  special="log",
                  arguments=list("base"=0,"u"=1))
parseSpecialNames(c("log(Z,u=1,base=3)","a","log(B,base=8,u=3)"),
                  special="log",
                  arguments=list("base"=0,"u"=1))
parseSpecialNames("treat(Z,u=2)",
                  special="treat",
                  arguments=list("u"=1,"k"=1))
parseSpecialNames(c("treat(Z,1,u=2)","treat(B,u=2,k=3)"),
                  special="treat",
                  arguments=list("u"=NA,"k"=NULL))
## does not work to set default to NULL:
parseSpecialNames(c("treat(Z,1,u=2)","treat(B,u=2)"),
                  special="treat",
                  arguments=list("u"=NA,"k"=NULL))
}
\author{
Thomas A. Gerds <tag@biostat.ku.dk>
}
\seealso{
model.design
}

