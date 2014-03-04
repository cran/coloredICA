\name{grad}
\alias{grad}
\title{
Gradient
}
\description{
This function evaluates the gradient of the objective function for the spectral density local maximum likelihood estimator.
}
\usage{
grad(x, omega, l_period, n, freq, h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Current estimate
}
  \item{omega}{
Frequency at which the spectral density estimate is evaluated.
}
  \item{l_period}{
Vector of length \code{n} with the log-periodogram evaluations at the \code{n} Fourier frequencies.
}
  \item{n}{
Number of points in the analyzed lattice.
}
  \item{freq}{
\code{n} \eqn{\times} 2 matrix with the \code{n} Fourier frequencies.
}
  \item{h}{
Kernel bandwidth.
}
}
\details{
In \code{\link{scICA}} function the maximization for the spectral density local maximum likelihood estimator is obtained through the Newton-Raphson algorithm. This function returns the gradient needed in the optimization method. See \code{\link{locmulti}} for further details.
}
\value{
It returns a gradient vector of length 3.
}
%\references{
%references
%}
\author{
Lee, S., Shen, H., Truong, Y. and Zanini, P.
}
\note{
It is auxiliary for \code{\link{scICA}} function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{scICA}}, \code{\link{locmulti}}, \code{\link{kern}}.
}

%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
