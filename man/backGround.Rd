\name{backGround}
\alias{backGround}
\title{
  Background and grid color control.
}
\description{
  Some users like background colors, and it may be helpful to have grid
  lines to read off e.g. probabilities from a Kaplan-Meier graph. Both
  things can be controlled with this function. However, it mainly serves \code{\link{plot.prodlim}}.
}
\usage{
backGround(xlim, ylim, bg="white", fg, horizontal = NULL, vertical = NULL, border = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xlim}{
    Limits for the xaxis, defaults to par("usr")[1:2].
  }
  \item{ylim}{
    Limits for the yaxis, defaults to par("usr")[3:4].
  }
  \item{bg}{
    Background color. Can be multiple colors which are then switched
    at each horizontal line.
  }
  \item{fg}{
    Grid line color.
  }
  \item{horizontal}{
    Numerical values at which horizontal grid lines are plotted.
  }
  \item{vertical}{
    Numerical values at which vertical grid lines are plotted.
  }
  \item{border}{
    The color of the border around the background.
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
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

plot(0,0)
backGround(bg="beige",fg="red",vertical=0,horizontal=0)

plot(0,0)
backGround(bg=c("yellow","green"),fg="red",xlim=c(-1,1),ylim=c(-1,1),horizontal=seq(0,1,.1))
backGround(bg=c("yellow","green"),fg="red",horizontal=seq(0,1,.1))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
