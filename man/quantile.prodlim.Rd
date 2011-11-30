\name{quantile.prodlim}
\alias{quantile.prodlim}
\title{
Quantiles for Kaplan-Meier and Aalen-Johannsen estimates.
}
\description{
  Quantiles for Kaplan-Meier and Aalen-Johannsen estimates.
}
\usage{
quantile.prodlim(x, q, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Object of class \code{"prodlim"}.
  }
  \item{q}{
    Quantiles. Vector of values between 0 and 1.
  }
  \item{\dots}{
    not used
  }
}
% \details{
% }
% \value{
% %%  ~Describe the value returned
% %%  If it is a LIST, use
% %%  \item{comp1 }{Description of 'comp1'}
% %%  \item{comp2 }{Description of 'comp2'}
% %% ...
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

% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{

\dontrun{
d=SimSurv(100)
f=prodlim(Hist(time,status)~X2,data=d)
# median
quantile(f,.5)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}