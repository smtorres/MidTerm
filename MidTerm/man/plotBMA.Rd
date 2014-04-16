\name{plotBMA}
\alias{plotBMA}
\title{plotBMA}
\usage{
  plotBMA(x, y, g, parallel = FALSE)
}
\arguments{
  \item{x}{A matrix of covariates to be included in the
  regressions}

  \item{y}{A vector of values with the outcome variable}

  \item{g}{A numeric value indicating the desired hyper-g
  prior}

  \item{parallel}{A logical vector specifying multi-core
  computations}
}
\value{
  A pdf with a graph per covariate showing a distribution
  of the coefficients computed in every model for each
  variable, the probability of the coefficient of being
  zero and the expected value of the coefficient.
}
\description{
  Plots regression statistics and Bayesian Model Averaging
  statistics
}
\examples{
cov.mat<-matrix(c(1,4,1,2,6,7,8,2,4,6,2,7,0,7,2,7,3,7,1,6), ncol=4)
out.vec<-cov.mat[,2]*3 + sample(1:5,5,replace=TRUE)
plotBMA(x=cov.mat, y=out.vec, g=3)
}
\author{
  Michelle Torres
}
\seealso{
  \code{\link{summBMA}} \code{\link{fitBMA}}
}

