\name{meanNeighbors}
\alias{meanNeighbors}
\title{
Helper function to obtain running means for prodlim objects.
}
\description{
Compute average values of a variable according to neighborhoods.
}
\usage{
meanNeighbors(x, y, ...)
}
\arguments{
  \item{x}{
    Object of class \code{"neighborhood"}.
  }
  \item{y}{
    Vector of numeric values.
}
  \item{\dots}{
    Not used.
  }
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
% \value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
% }
% \references{
%% ~put references to the literature/web site here ~
% }
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{neighborhood}}
}
\examples{
meanNeighbors(x=1:10,y=c(1,10,100,1000,1001,1001,1001,1002,1002,1002))
}
\keyword{survival }
