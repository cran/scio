\name{scio.refit}
\alias{scio.refit}
\title{Refitted Sparse Column-wise Inverse Operator}
\description{
  Refitted SCIO Estimators using Penalized Likelihood
}
\usage{
scio.refit(S, Omega,thr=1e-4, pkg = c("QUIC","glasso"),...)
}
\arguments{
\item{S}{Input covariance matrix of size p by p (symmetric).}
\item{Omega}{Estimated inverse covariance matrices. Can be a matrix of
  size p by p from \code{scio} or a collection of matrices from
  \code{sciopath}.}
\item{thr}{Tolerance. Small entries in magnitude  (<\code{thr}) in \code{Omega}
 are treated as zeros when refitting the precision matrix with the same
 support as the \code{scio} or \code{sciopath} outputs.  Default 1e-4.}
\item{pkg}{R packge to be used for refitting. Default \code{QUIC}.}
\item{...}{Additional options passed on to \code{QUIC}, which is the
  only likelihood solver called in the current release. More solvers
  will be included in future releases.}
}

\details{
  This implements the refitting procedure discussed in Cai, Liu, and Luo
  (2011). The current version uses the \code{QUIC} solver for the
  penalized likelihood criterion. More solvers will be added.
}
\value{
  A list with one component:
  \item{w}{Estimated inverse covariance matrix when a single
    \code{Omega} matrix is supplied, or an array of matrices when a 3
    dimensional array 
    of \code{Omega} is supplied.}
}

\references{
  Weidong Liu and Xi Luo (2012). Fast and Adaptive Sparse Precision
  Matrix Estimation in High Dimensions. arXiv:1203.3896.

  Tony Cai, Weidong Liu, and Xi Luo (2011). A Constrained L1 Minimization Approach to Sparse Precision Matrix Estimation. Journal of the American Statistical Association, 106(494), 594-607.
}
 
\examples{
set.seed(100)
x<-matrix(rnorm(50*20),ncol=4)
s<- var(x)
a<-scio(s, lambda=.01)

require(glasso) 
b <- scio.refit(s, a$w)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

