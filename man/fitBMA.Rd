\name{fitBMA}
\alias{fitBMA}
\title{fitBMA}
\usage{
fitBMA(x, y, g = NULL, parallel = FALSE)
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
A list with the following elements \item{Stats}{A matrix
with the coefficients and R-squared obtained from all
regressions} \item{PostModOdds}{A matrix containing the
posterior model odds for each model} \item{PostExp}{A
vector with the posterior expected value of each
coefficient} \item{PostProb}{A vector with the posterior
probability that each coefficient is non-zero}
}
\description{
Computes regression statistics and Bayesian Model Averaging
statistics
}
\examples{
cov.mat<-matrix(c(1,4,1,2,6,7,8,2,4,6,2,7,0,7,2,7,3,7,1,6), ncol=4)
out.vec<-cov.mat[,2]*3 + sample(1:5,5,replace=TRUE)
fitBMA(x=cov.mat, y=out.vec, g=3)
}
\author{
Michelle Torres
}
\seealso{
summBMA plotBMA
}

