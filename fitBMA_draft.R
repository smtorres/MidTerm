fitBMA<-function(x,y,g=NULL){
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
  for (i in 1:length(vars)){
    covs.coeffs[[i]]<-apply(vars[[i]], 2, FUN=covs.store, ivs=x, dv=y)
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
return(PostOdds)
}