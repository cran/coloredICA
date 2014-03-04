\name{rerow}
\alias{rerow}
\title{
Standardization of a matrix reordering the rows
}
\description{
This function reorders and standardizes the rows of a matrix.
}
\usage{
rerow(w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{
Matrix to standardize.
}
}
\details{
The standardization is done in such a way that every row has length 1, the largest absolute value of the row has a positive sign and the rows are ordered decreasingly according to their largest value.
}
\value{
It returns the standardized matrix.
}
%\references{
%references
%}
\author{
Lee, S., Shen, H., Truong, Y. and Zanini, P.
}
\note{
It is auxiliary for \code{\link{cICA}} and \code{\link{scICA}} functions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{amari_distance}}, \code{\link{cICA}}, \code{\link{scICA}}
}

%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
