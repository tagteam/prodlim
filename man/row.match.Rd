\name{row.match}
\alias{row.match}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identifying rows in a matrix or data.frame}
\description{
  Function for finding matching rows between two matrices or
  data.frames. First the matrices or data.frames are vectorized by row
  wise pasting together the elements. Then it uses the function match.
  Thus the function returns a vector with the row numbers of (first) matches of
  its first argument in its second.
}
\usage{
row.match(x, table, nomatch = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector or matrix whose rows are to be matched}
  \item{table}{Matrix or data.frame that contain the rows to be matched against.}
  \item{nomatch}{the value to be returned in the case when no match is found. 
          Note that it is coerced to 'integer'.}
}
\value{
  A vector of the same length as 'x'.
}
\author{Thomas A. Gerds}
\seealso{\code{match}}
\examples{
tab <- data.frame(num=1:26,abc=letters)
x <- c(3,"c")
row.match(x,tab)
x <- data.frame(n=c(3,8),z=c("c","h"))
row.match(x,tab)
}
\keyword{misc}