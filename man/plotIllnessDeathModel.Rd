\name{plotIllnessDeathModel}
\alias{plotIllnessDeathModel}
\title{
Plotting an illness-death-model.
}
\description{
Plotting an illness-death-model using \code{plot.Hist}.
}
\usage{
plotIllnessDeathModel(stateLabels, style = 1, recovery = FALSE, ...)
}
\arguments{
  \item{stateLabels}{
    Labels for the three boxes.
  }
  \item{style}{
    Either \code{1} or anything else, switches 
    the orientation of the graph. Hard to explain in words, see examples.
  }
  \item{recovery}{
    Logical. If \code{TRUE} 
    there will be an arrow from the illness state to the initial state.
  }
  \item{\dots}{
    Arguments passed to plot.Hist.
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
% \value{
% %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
% }
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

 \seealso{
\code{\link{plotCompetingRiskModel}}, \code{\link{plot.Hist}}
 }
\examples{
plotIllnessDeathModel()
plotIllnessDeathModel(style=2)
plotIllnessDeathModel(style=2,stateLabels=c("a","b\nc","d"),box1.col="yellow",box2.col="green",box3.col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival }
