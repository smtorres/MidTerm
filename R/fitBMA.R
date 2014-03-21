#' fitBMA
#' 
#' Computes regression statistics and Bayesian Model Averaging statistics
#' 
#' @param x A matrix of covariates to be included in the regressions
#' @param y A vector of values with the outcome variable 
#' @param g A numeric value indicating the desired hyper-g prior
#' @param parallel A logical vector specifying multi-core computations
#' 
#' @return A list with the following elements
#' \item{Stats}{A matrix with the coefficients and R-squared obtained from all regressions}
#' \item{PostModOdds}{A matrix containing the posterior model odds for each model}
#' \item{PostExp}{A vector with the posterior expected value of each coefficient}
#' \item{PostProb}{A vector with the posterior probability that each coefficient is non-zero}
#' 
#' @author Michelle Torres
#' 
#' @examples
#' 
#' cov.mat<-matrix(c(1,4,1,2,6,7,8,2,4,6,2,7,0,7,2,7,3,7,1,6), ncol=4)
#' out.vec<-cov.mat[,2]*3 + sample(1:5,5,replace=TRUE)
#' fitBMA(x=cov.mat, y=out.vec, g=3)
#' 
#' @seealso \code{\link{summBMA}} \code{\link{plotBMA}}
#' 
#' @export
  fitBMA<-function(x,y,g=NULL, parallel=FALSE){
  #Standardize covariates
  cov<-apply(x,2,FUN=function(x)(x-mean(x))/sd(x))
  num.cov<-ncol(x)
  num.obs<-nrow(x)
  #Standardize outcome variable
  out<-(y-mean(y))/sd(y)
  #Packages required
  require(combinat)
  require(plyr)
  #Empty list: vars --> possible combinations
  vars<-list()
  #Add to "vars" all the possible combination of covariates
  for (k in 1:num.cov){
    vars[[k]]<-combn(num.cov,k)
  }
  vars[[num.cov]]<-as.matrix(vars[[num.cov]])
  #Empty list: covs.coeffs --> coefficients from regression
  covs.coeffs<-list()
  #FUNCTION: covs.store (matrix covariates, dependent variable, index coefficients to compute)
  covs.store<-function(ivs, dv, index){  
    #Create names of variables "x" and assign them to matrix
    cov.names<-paste0("x", 1:ncol(ivs))
    colnames(ivs)<-cov.names
    #Vector of names of covariates to include in the regression and new dependent variable
    xnam<-paste0("x", index)
    dep<-dv
    #Create vectors with x values per variable
    for (i in 1:ncol(ivs)){
      assign(paste("x",i,sep=""), ivs[,i])
    }
    #Create the formula to be computed inside the function
    (fmla <- as.formula(paste("dep ~ ", paste(xnam, collapse= "+"), "-1")))
    #Run regression and store coefficients in "coeffs"
    coeffs<-summary(lm(fmla))$coefficients[,1]
    #Create a list "output" with the coefficients obtained 
    output<-t(as.matrix(coeffs))
    colnames(output)<-xnam
    output<-as.data.frame(output)
    return(output)
  }
  #Apply function "covs.store" to the matrix "x" based on the combinations in "vars" 
  if(parallel==TRUE){
  for (i in 1:length(vars)){
    covs.coeffs[[i]]<-alply(vars[[i]], .margins=2, .fun=covs.store, ivs=x, dv=y, .parallel=TRUE)
  }
  }
  
  if(parallel==FALSE){
    for (i in 1:length(vars)){
      covs.coeffs[[i]]<-alply(vars[[i]], .margins=2, .fun=covs.store, ivs=x, dv=y)
    }
  }
  #Collapse elements of a list (Level 2)
  reg.coeffs<-lapply(covs.coeffs, rbind.fill)
  #Collapse elements of a list (Level 1)
  reg.coeffs<-rbind.fill(reg.coeffs)
  #FUNCTION: r.store (Run regressions and stores R.SQUARED)
  r.store<-function(ivs, dv, index){
    # Create names of covariates and assign to matrix
    cov.names<-paste0("x", 1:ncol(ivs))
    colnames(ivs)<-cov.names
    # Create the name of the selected variables to be included
    xnam<-paste0("x", index)
    dep<-dv
    #Create vectors of valuex x for each variable
    for (i in 1:ncol(ivs)){
      assign(paste("x",i,sep=""), ivs[,i])
    }
    #Formula to compute regression
    (fmla <- as.formula(paste("dep ~ ", paste(xnam, collapse= "+"),"-1")))
    #Store R-squared from regression
    r2<-summary(lm(fmla))$r.squared
    return(r2)
  }
  #Empty list: r.reg --> Vector with R-squared values for each regression
  r.reg<-list()
  for (i in 1:length(vars)){
    r.reg[[i]]<-apply(vars[[i]], 2, FUN=r.store, ivs=x, dv=y)
  }
  r.reg<-as.matrix(unlist(r.reg))
  colnames(r.reg)<-"R-squared"
  #Final list with matrix of coefficients and vector of R-squared values
  results.model<-cbind(reg.coeffs, r.reg)
  
  ##BMA--> z: vector with regression results + r^2, g: prior, n: number of observations
  BMA.1<-function(z,g,n){
    mod<-z[is.na(z)==FALSE]
    p<-length(mod)-1
    R2<- tail(as.numeric(z), n=1)
    #B[M_k : M_0]
    base.mod<-(1+g)^((n-p-1)/2) * (1+g*(1-R2))^(-(n-1)/2)
    return(base.mod)
  }
  
  #Apply BMA.1 to the matrix with all possible models
  BMk<-NULL
  BMk<-apply(results.model, 1, FUN=BMA.1, g=g, n=num.obs)
  #Calculation of p(M_k|Y)
  PostOdds<-BMk/sum(BMk) 
  #Matrix with Coefficients, R-squared and Post Odds
  results.model.post<-cbind(results.model, PostOdds)

  #BMA.2--> w: vector with regresion results + r^2, g: prior
  BMA.2<-function(w,g){
    #Create a vector with only the coefficients in the model
    mod<-w[is.na(w)==FALSE]
    length.w<-length(as.numeric(mod))
    coeffs.vars<-w[1:(length(w)-2)]
    E.Betak<-(g/(g+1))*coeffs.vars
      ##Posterior expected value
        post.odd<-as.numeric(w[length(w)])
        E.Beta.Y.pre<-post.odd*E.Betak
    return(E.Beta.Y.pre)
  }
  
  #Apply to all the models
  E.Beta.Y<-apply(results.model.post,1,FUN=BMA.2, g=g)
  E.Beta.Y<-apply(E.Beta.Y,1,FUN=function(x) sum(x, na.rm=TRUE))

  #Posterior Probability
  post.prob<-NULL
  for (i in 1:(ncol(results.model.post)-2)){
    post.prob[i]<-sum(PostOdds[is.na(results.model.post[,i])==FALSE])
  }
  #Organize final output
    name.models<-paste("Model",1:nrow(results.model))
    rownames(results.model)<-name.models
    PostOdds<-as.matrix(PostOdds)
    rownames(PostOdds)<-name.models
    colnames(PostOdds)<-"PostOdds"
    post.prob<-t(as.matrix(post.prob))
    colnames(post.prob)<-paste0("x", 1:ncol(x))
    final.results<-list(results.model, PostOdds, E.Beta.Y, post.prob)
  names(final.results)<-c("Stats", "PostModOdds", "PostExp", "PostProb")
  return(final.results)
}
