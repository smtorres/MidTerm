setwd("/Users/michelletorres/Desktop/MidTerm")

library(devtools)
library(roxygen2)
current.code <- as.package("MidTerm")
load_all(current.code)
document(current.code)

check(current.code)
install(pkg=current.code,local=TRUE)
build(current.code,path=getwd())


cov.mat<-matrix(c(1,4,1,2,6,7,8,2,4,6,2,7,0,7,2,7,3,7,1,6), ncol=4)
 out.vec<-cov.mat[,2]*3 + sample(1:5,5,replace=TRUE)
temp <- fitBMA(x=cov.mat, y=out.vec, g=3)
summBMA(temp)
