\name{scICA}
\alias{scICA}
\title{
spatial colored Independent Component Analysis
}
\description{
This function implements the spatial colored Independent Component Analysis (scICA) algorithm, where sources are treated as spatial stochastic processes on a lattice.
}
\usage{
scICA(Xin, M = dim(Xin)[1], Win = diag(M), tol = 1e-04, maxit = 20, nmaxit = 1, 
unmixing.estimate = "eigenvector", n1, n2, nx01 = n1, nx02 = n2, h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xin}{
Data matrix with \code{p} rows (representing variables) and \code{n} columns (representing observations).
}
  \item{M}{
Number of components to be extracted.
}
  \item{Win}{
Initial guess for the unmixing matrix \code{W}. Dimensions need to be \code{M x M}.
}
  \item{tol}{
Tolerance used to establish the convergence of the algorithm.
}
  \item{maxit}{
Maximum number of iterations.
}
  \item{nmaxit}{
If the algorithm does not converge, it is run again with a new initial guess for the unmixing matrix \code{W}. This operation is done \code{nmaxit} times.
}
  \item{unmixing.estimate}{
The method used in the unmixing matrix estimation step. The two allowed choices are \code{"eigenvector"} and \code{"newton"} (see Details).
}
  \item{n1}{
Number of rows of the lattice.
}
  \item{n2}{
Number of columns of the lattice.
}
  \item{nx01}{
Number of rows of the lattice where the spectral density is evaluate. Default value is \code{n1}.
}
  \item{nx02}{
Number of columns of the lattice where the spectral density is evaluate. Default value is \code{n2}.
}
  \item{h}{
Kernel bandwidth used for the nonparametric estimation of the sources spectral densities.
}
}
\details{
In the Independent Component Analysis approach, the data matrix \eqn{X} is considered to be a linear combination of independent components, i.e. \eqn{X = AS}, where rows of \eqn{S} contain the unobserved realizations of the independent components and \eqn{A} is a linear mixing matrix. According to classical ICA procedures data matrix \eqn{X} is centered and, then, whitened by projecting the data onto its principal component directions, i.e. \eqn{X \rightarrow KX = \widetilde{X}} where \eqn{K} is a \code{M x p} pre-whitening matrix. The scICA algorithm then estimates the unmixing matrix \eqn{W}, with \eqn{W\widetilde{X} = S}, according to the procedure described below. Then, defining \eqn{\widetilde{W}=WK}, the mixing matrix \eqn{A} is recovered through \eqn{A=\widetilde{W}^T(\widetilde{W}\widetilde{W}^T)^{-1}}.

Spatial colored Independent Component Analysis assumes that the independent sources are spatial stochastic processes on a lattice. To perform ICA, the Whittle log-likelihood is exploited. In particular the log-likelihood is written in function of the unmixing matrix \eqn{W} and the spectral densities \eqn{f_{S_j}} of the spatial autocorrelated sources as follows:
\deqn{l(W,\boldsymbol{f}_{\boldsymbol{S}};\widetilde{X})=\sum_{j=1}^p\sum_{k=1}^n\left(\frac{\boldsymbol{e}_j^T W \widetilde{\boldsymbol{f}}(r_k,\widetilde{X})W^T \boldsymbol{e}_j}{f_{S_j}(r_k)}+\ln f_{S_j}(r_k)\right) + n\ln|\det(W)|.}
Due to whitening, \eqn{W} is orthogonal and the last term of the objective function can be dropped. The orthogonality of the unmixing matrix \eqn{W} can be imposed in two different ways, setting the argument \code{unmixing.estimate}. In this way the estimate of the unmixing matrix \eqn{W} can be found according two different procedures:
\itemize{
\item as described in Shen et al. (2014). A penalty term is added to the objective function. In particular \eqn{\tau\bold{w}'_{j}C_{j}\bold{w}_{j}}, where \eqn{\bold{w}'_{j}} is the \eqn{j}th column of \eqn{W}, \eqn{C_j=\sum_{k\neq j}\bold{w}_k\bold{w}'_k} and \eqn{\tau} is a tuning parameter. The matrix \eqn{C_j} provides an orthogonality constraint in the sense that \eqn{\bold{w}'_{j}C_{j}\bold{w}_{j}=\sum_{k\neq j}}. In this way the objective function assumes a symmetric and positive-definite form and the argmin correspond to the lower eigenvalue. This choice is obtained setting \code{unmixing.matrix = "eigenvector"}.
\item as described in Lee et al. (2011). The orthogonality constraint is considered performing the minimization of the objective function according a Newton-Raphson method with Lagrange multiplier. This choice is obtained setting \code{unmixing.matrix = "newton"}.
}
Independently from the choice of the technique to minimize the objective function, the scICA algorithm is based on an iterative procedure. While the Amari error is greater than \code{tol} and the number of iteration is less or equal than \code{maxit}, the two following steps are repeated:
\itemize{
\item nonparametric estimation of the sources spectral density through a multidimensional local linear kernel estimator \eqn{\widehat{m}_{LK}} (see Shen et al. (2014) for further details).
\item estimate the unmixing matrix \eqn{W} according the method selected in \code{unmixing.estimate}.
}
}
\value{
A list containing the following components:
\item{W}{Estimate of the \code{M x M} unmixing matrix in the whitened space.}
\item{K}{pre-whitening matrix that projects data onto the first \code{M} principal components. Dimensions are \code{M x p}.}
\item{A}{Estimate of the \code{p x M} mixing matrix.}
\item{S}{Estimate of the \code{M x n} source matrix.}
\item{X}{Original \code{p x n} data matrix.}
\item{iter}{number of iterations.}
\item{NInv}{number of times the algorithm is rerun after it does not achieve convergence.}
\item{den}{Estimate of the spectral density of the sources. Dimensions are \code{M x n}.}
}
\references{
Shen, H., Truong, Y., Zanini, P. (2014). Independent Component Analysis for Spatial Processes on a Lattice. \emph{MOX report 38/2014}, Department of Mathematics, Politecnico di Milano.

Lee, S., Shen, H., Truong, Y., Lewis, M., Huang, X. (2011). Independent Component Analysis Involving Autocorrelated Sources With an Application to Funcional Magnetic Resonance Imaging. \emph{Journal of the American Statistical Association}, \bold{106}, 1009--1024.
}
\author{
Lee, S., Shen, H., Truong, Y. and Zanini, P.
}
\note{
Note that source matrix \code{S} and spectral density matrix \code{den} are \code{n x M} matrices. Every row, that should be in a \code{n1 x n2} grid, has been vectorized in a \code{n} vector by column, with \code{n = n1 x n2}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cICA}}
}

\examples{
## Not run:

require(fastICA)

n1=20
n2=20
M=2

# Fist source

sigma1=2
S1=matrix(0,n1,n2)
for (i in 1:n1){
S1[i,]=rnorm(n2,i*2,0.2)
}
for (j in 1:n2){
S1[,j]=S1[,j]+rnorm(n1,j*2,0.2)
}
S1=S1+matrix(rnorm(n1*n2,0,sigma1),n1,n2)

image(1:n2,1:n1,t(S1[n1:1,]),xlab="",ylab="",main="Source 1")
contour(1:n2,1:n1,t(S1[n1:1,]),add=TRUE)

# Second source

val1=1
val2=1.2
val3=1.5
val4=2
sigma2=0.1 

S2=matrix(0,n1,n2)
S2[2:5,4:10]=val1
S2[3:4,5:9]=val3
S2[13:18,16:19]=val2
S2[14:17,17:18]=val4
S2=S2+matrix(rnorm(n1*n2,0,sigma2),n1,n2)

image(1:n2,1:n1,t(S2[n1:1,]),xlab="",ylab="",main="Source 2")
contour(1:n2,1:n1,t(S2[n1:1,]),add=TRUE)

# Generating data matrix X

A = rerow(matrix(runif(M^2)-0.5,M,M))
W = solve(A)
S=rbind(as.vector(S1),as.vector(S2))
X = A \%*\% S

# Solving Blind Source Separation problem with three different methods

cica = cICA(X,tol=0.001)
## scica = scICA(X,n1=n1,n2=n2,h=(2*pi/10),tol=0.001)
fica = fastICA(t(X),2)

amari_distance(t(A),t(cica$A))
## amari_distance(t(A),t(scica$A))
amari_distance(t(A),fica$A)

Shat1=cica$S
## Shat2=scica$S
Shat3=t(fica$S)

par(mfrow=c(2,2))
image(t(S1[n1:1,]),xlab="",ylab="")
contour(t(S1[n1:1,]),add=TRUE)
image(t(S2[n1:1,]),xlab="",ylab="")
contour(t(S2[n1:1,]),add=TRUE)
image(t(matrix(Shat1[1,],n1,n2)[n1:1,]),xlab="",ylab="")
contour(t(matrix(Shat1[1,],n1,n2)[n1:1,]),add=TRUE)
image(t(matrix(Shat1[2,],n1,n2)[n1:1,]),xlab="",ylab="")
contour(t(matrix(Shat1[2,],n1,n2)[n1:1,]),add=TRUE)

## par(mfrow=c(2,2))
## image(t(S1[n1:1,]),xlab="",ylab="")
## contour(t(S1[n1:1,]),add=TRUE)
## image(t(S2[n1:1,]),xlab="",ylab="")
## contour(t(S2[n1:1,]),add=TRUE)
## image(t(matrix(Shat2[1,],n1,n2)[n1:1,]),xlab="",ylab="")
## contour(t(matrix(Shat2[1,],n1,n2)[n1:1,]),add=TRUE)
## image(t(matrix(Shat2[2,],n1,n2)[n1:1,]),xlab="",ylab="")
## contour(t(matrix(Shat2[2,],n1,n2)[n1:1,]),add=TRUE)

par(mfrow=c(2,2))
image(t(S1[n1:1,]),xlab="",ylab="")
contour(t(S1[n1:1,]),add=TRUE)
image(t(S2[n1:1,]),xlab="",ylab="")
contour(t(S2[n1:1,]),add=TRUE)
image(t(matrix(Shat3[1,],n1,n2)[n1:1,]),xlab="",ylab="")
contour(t(matrix(Shat3[1,],n1,n2)[n1:1,]),add=TRUE)
image(t(matrix(Shat3[2,],n1,n2)[n1:1,]),xlab="",ylab="")
contour(t(matrix(Shat3[2,],n1,n2)[n1:1,]),add=TRUE)

## End (Not run)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
