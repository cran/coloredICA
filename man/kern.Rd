\name{kern}
\alias{kern}
\title{
Kernel evaluation
}
\description{
This function returns a vector of length \eqn{n} with the evaluation of an exponential Kernel at a point \code{x0} \eqn{\in R^2} for the \eqn{n} Fourier frequencies.
}
\usage{
kern(x0, h, freq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{
Point \eqn{\in R^2} at which the Kernel is evaluated.
}
  \item{h}{
Kernel bandwidth.
}
  \item{freq}{
\code{n} \eqn{\times} 2 matrix with the \code{n} Fourier frequencies, where \code{n} is the number of points in the analyzed lattice.
}
}
\details{
This function returns a vector of length \eqn{n} with the evaluation of an exponential Kernel at a point \code{x0} \eqn{\in R^2} for the \eqn{n} Fourier frequencies. In particular the \eqn{k}-th vector element is
\deqn{e^{-\frac{(\boldsymbol{\omega}_k - x_0)^T(\boldsymbol{\omega}_k - x_0)}{h^2}}.}
}
\value{
It returns a list containing the following component:
\item{v}{vector containing the \eqn{n} exponential Kernel evaluations.}
}
%\references{
%references
%}
\author{
Lee, S., Shen, H., Truong, Y. and Zanini, P.
}
\note{
It is auxiliary for \code{\link{locmulti}}, \code{\link{grad}} and \code{\link{hess}} functions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{scICA}}, \code{\link{locmulti}}, \code{\link{grad}}, \code{\link{hess}}.
}

%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
