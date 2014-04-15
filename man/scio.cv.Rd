\name{scio.cv}
\alias{scio.cv}
\title{Sparse Column-wise Inverse Operator}
\description{
Cross validated estimates of a sparse inverse covariance matrix using Sparse Column-wise
Inverse Operator
}
\usage{
scio.cv(X, lambda.max=1, alpha=0.95, cv.maxit=1e2, ...)
}
\arguments{
\item{X}{Input data of dimension samples (n)  x variables (p).}
\item{lambda.max}{Maximum lambda to start with in CV, which is decreased
  by mutliplying \code{alpha} in each iteration.}
\item{alpha}{Scaling factor to decrease \code{lambda} by
  multiplication.} 
\item{cv.maxit}{Maximum number of CV iterations. Default 1e2.}
\item{...}{Other option parameters in \code{scio}.}
}

\details{
  This is a fast, nonparametric approach to estimate sparse inverse covariance
  matrices, with possibly really large dimensions.  Details of this procedure are
  described in the reference.

  This function does a simple cross validation based on likelihood. 
}
\value{
  A list with components:
  \item{w}{Estimated inverse covariance matrix}
  \item{lambda.cv}{CV selected lambda}
}

\references{
  Weidong Liu and Xi Luo (2012). Fast and Adaptive Sparse Precision
  Matrix Estimation in High Dimensions. arXiv:1203.3896.
}
\examples{
set.seed(100)
x<-matrix(rnorm(50*20),ncol=4)
a<-scio.cv(x)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

