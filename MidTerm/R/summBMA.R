#' summBMA
#' 
#' Prints a table with Bayesian Model Averaging statistics
#' 
#' @param x A matrix of covariates to be included in the regressions
#' @param y A vector of values with the outcome variable 
#' @param g A numeric value indicating the desired hyper-g prior
#' @param parallel A logical vector specifying multi-core computations
#' 
#' @return A table with the following elements for the coefficient of each variable in the covariate matrix
#' \item{PostExp}{Number indicating the posterior expected value of each coefficient}
#' \item{PostProb}{Number indicating the posterior probability that each coefficient is non-zero}
#' 
#' @author Michelle Torres
#' 
#' @examples
#' 
#' cov.mat<-matrix(c(1,4,1,2,6,7,8,2,4,6,2,7,0,7,2,7,3,7,1,6), ncol=4)
#' out.vec<-cov.mat[,2]*3 + sample(1:5,5,replace=TRUE)
#' summBMA(x=cov.mat, y=out.vec, g=3)
#' 
#' @seealso \code{\link{fitBMA}} \code{\link{plotBMA}}
#' 
#' @export
summBMA<-function(x,y,g, parallel=FALSE){
  if(parallel==TRUE){output<-fitBMA(x,y,g, parallel=TRUE)}
  if(parallel==FALSE){output<-fitBMA(x,y,g)}
  Coefficient<-c("PostExp", "PostProb != 0")
  coeff.names<-paste0("x", 1:ncol(x))
  post.exp<-as.matrix(output[[3]])
  post.prob<-t(as.matrix(output[[4]]))
  table.results<-as.data.frame(cbind(post.exp, post.prob))
  colnames(table.results)<-Coefficient
  print(table.results)
}