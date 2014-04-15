\name{scio}
\alias{scio}
\title{Sparse Column-wise Inverse Operator}
\description{
Estimates a sparse inverse covariance matrix using Sparse Column-wise
Inverse Operator
}
\usage{
scio(S, lambda,  thr=1e-4, maxit=1e4, pen.diag=F, sym=T)
}
\arguments{
\item{S}{Input covariance matrix of size p by p (symmetric).}
\item{lambda}{(Non-negative) regularization parameter for the lasso
  penalty.  Can be a scalar or a matrix of size p by p.}
\item{thr}{Threshold for convergence. Iterations stop when the maximum
  change in two successive updates is less than \code{thr}. Default value is 1e-4.}
\item{maxit}{Maximum number of iterations for each column computation. Default 10,000.}
\item{pen.diag}{Whether the diagonal should be penalized.  Default
  False.}
\item{sym}{Whether the return values should be symmetrized.  Default True.}
}

\details{
  This is a fast, nonparametric approach to estimate sparse inverse covariance
  matrices, with possibly really large dimensions.  Details of this procedure are
  described in the reference.
}
\value{
  A list with components:
  \item{w}{Estimated inverse covariance matrix}
}

\references{
  Weidong Liu and Xi Luo (2012). Fast and Adaptive Sparse Precision
  Matrix Estimation in High Dimensions. arXiv:1203.3896.
}
\examples{
set.seed(100)
x<-matrix(rnorm(50*20),ncol=4)
s<- var(x)
a<-scio(s, lambda=.01)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

