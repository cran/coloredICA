\name{amari_distance}
\alias{amari_distance}
\title{
Amari error
}
\description{
This function measures the Amari error between two matrices.
}
\usage{
amari_distance(Q1, Q2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q1}{
First matrix.
}
  \item{Q2}{
Second matrix.
}
}
\details{
The Amari error \eqn{D(Q1|Q2)} between two M x M matrices \eqn{Q1} and \eqn{Q2} is evaluated through

\deqn{D(Q1|Q2)=\frac{1}{2M(M-1)}\sum_{j=1}^M\left(\frac{\sum_{i}|a_{ij}|}{max_{i}|a_{ij}|}-1\right)+\frac{1}{2M(M-1)}\sum_{i=1}^M\left(\frac{\sum_{j}|a_{ij}|}{max_{j}|a_{ij}|}-1\right),}

where \eqn{Q2} is invertible and \eqn{a_{ij}} is the \eqn{ij}th element of \eqn{Q1Q2^{-1}}.

}
\value{
It returns the Amari error between two matrices \eqn{Q1} and \eqn{Q2}.
}
\references{
Amari, S., Cichocki, A., Yang, H. et al. (1996). A New Learning Algorithm for Blind Signal Separation. \emph{Advances in Neural Information Processing Systems}, \bold{8}, 757--763

Bach, F., Jordan, M. (2003). Kernel Independent Component Analysis. \emph{Journal of Machine Learning Research}, \bold{3}, 1--48
}
\author{
Lee, S., Shen, H., Truong, Y. and Zanini, P.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cICA}}, \code{\link{scICA}}
}

\examples{

M <- 4
A <- matrix(rnorm(M*M),M,M)
B <- matrix(rnorm(M*M),M,M)

amari_distance(A,B)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
