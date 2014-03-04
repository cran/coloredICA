\name{locmulti}
\alias{locmulti}
\title{
Local polynomial spectral density estimation
}
\description{
This function implements a local polinomial estimation for the log spectral density at a point \code{x0} \eqn{\in R^2}.
}
\usage{
locmulti(x0, l_period, n, freq, h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{
Point \eqn{\in R^2} at which the spectral density estimate is evaluated.
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
\code{\link{locmulti}} function is auxiliary for the nonparametric estimation of the sources spectral density step of the \code{\link{scICA}} function. \code{\link{locmulti}} function implements the initial estimates for the local maximum likelihood estimator of the log spectral density \eqn{m(\code{x0})} at a point \code{x0} \eqn{\in R^2}. To obtain an estimate of \eqn{m(\code{x0})} the local likelihood function
\deqn{\sum_{k}\left(Y_k - a - \textbf{b}'(\boldsymbol{\omega}_k - x0) - e^{Y_k - a - \textbf{b}'(\boldsymbol{\omega}_k - x0)} \right)K_H(\boldsymbol{\omega}_k - x0)}
is constructed, where \eqn{Y_k} denotes the log-periodogram value at the Fourier frequency \eqn{\boldsymbol{\omega}_k}, \eqn{K_H} a surface kernel and \eqn{H=(h,h)}. The local maximum estimator \eqn{\widehat{m}_{LK}(\code{x0})} is \eqn{\widehat{a}} in the maximizer \eqn{(\widehat{a},\widehat{\textbf{b}})}. The estimate is implemented directly in the \code{\link{scICA}} function through a Newton-Rapshon algorithm. The initialization for the Newton-Rapshon algorithm is derived through a local polynomial approximation implemented in this \code{\link{locmulti}} function. In particular the following function is minimized to find a local polynomial approximation for \eqn{m(\code{x0})}
\deqn{\sum_{k}\left(Y_k - a - \textbf{b}'(\boldsymbol{\omega}_k - x0) \right)^{2}K_H(\boldsymbol{\omega}_k - x0)}
and the minimizer \eqn{\widehat{a}} is used as an initial value in order to obtain the local maximum likelihood estimator \eqn{\widehat{m}_{LK}(\code{x0})}.
}
\value{
It returns a list containing the following component:
\item{ahat}{local polynomial estimate of the log spectral density at \code{x0}.}
}
\references{
Shen, H., Truong, Y., Zanini, P. (2014). Independent Component Analysis for Spatial Processes on a Lattice. \emph{MOX report 38/2014}, Department of Mathematics, Politecnico di Milano.

Fan, J., Kreutzberger, E. (1998). Automatic Local Smoothing for Spectral Density Estimation. \emph{Scandinavian Journal of Statistics}, \bold{25}, 359--369.
}
\author{
Lee, S., Shen, H., Truong, Y. and Zanini, P.
}
\note{
It is auxiliary for \code{\link{scICA}} function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{scICA}}, \code{\link{kern}}
}

%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
