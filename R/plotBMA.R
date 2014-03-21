#'plotBMA
#' 
#' Plots regression statistics and Bayesian Model Averaging statistics
#' 
#' @param x A matrix of covariates to be included in the regressions
#' @param y A vector of values with the outcome variable 
#' @param g A numeric value indicating the desired hyper-g prior
#' @param parallel A logical vector specifying multi-core computations
#' 
#' @return A pdf with a graph per covariate showing a distribution of the coefficients computed in every model for each variable, the probability of the coefficient of being zero and the expected value of the coefficient.
#' 
#' 
#' @author Michelle Torres
#' 
#' @examples
#' 
#' cov.mat<-matrix(c(1,4,1,2,6,7,8,2,4,6,2,7,0,7,2,7,3,7,1,6), ncol=4)
#' out.vec<-cov.mat[,2]*3 + sample(1:5,5,replace=TRUE)
#' plotBMA(x=cov.mat, y=out.vec, g=3)
#' 
#' @seealso \code{\link{summBMA}} \code{\link{fitBMA}}
#' 
#' @export
plotBMA<-function(x,y,g,parallel=FALSE){
    nvars<-ncol(x) 
    
   if(parallel==TRUE){output<-fitBMA(x,y,g, parallel=TRUE)}
   if(parallel==FALSE){output<-fitBMA(x,y,g)}
  #PLOT
  bm.plot<-function(x,z,w){
    nm<-colnames(x)
    x<-as.matrix(x)
    len<-nrow(x)
    prob<-z
    post.exp<-w
    x2<-x[1:len,]
    x2<-x[!is.na(x)==TRUE]
    x.lim1<-min(x2,w,0.5)-1 
    x.lim2<- max(x2+1)
    y.lim2<-(max(density(x2)$y, 1-prob))+0.1
    plot(density(x2), main=nm, xlim=c(x.lim1, x.lim2),ylim=c(0,y.lim2))
    segments(0,0,0,1-prob, lwd=2, col="red")
    abline(v=w, lty=2)
    legend(x.lim2-(.33*abs(x.lim2-x.lim1)), y.lim2, legend=c("Coefficients", "1-Post. prob", "Expected value"),
    lty=c(1,1,2), col=c("black", "red", "black"))
  }
  pdf(file='plotsBMA.pdf')
  for(i in 1:nvars){
    bm.plot(output[[1]][i],output[[4]][i],output[[3]][i])
  }
  #apply(mat.coeffs, 2, FUN=bm.plot)
  dev.off()
  print("PDF document plotsBMA.pdf created")
}