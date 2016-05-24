mean.impute.missing <- function(m) {
 mu <- rowMeans(m,na.rm=T)
 isna <- is.na(m)
 m[isna] <- matrix(mu,nrow(m),ncol(m))[isna]
 return(list(m=m,mu=mu))
}
if(!require("irlba")){
	install.packages("irlba",repos="http://cran.r-project.org",lib="./Rlib/")
	library(irlba,lib="./Rlib/")
}else{
	library(irlba)
}
##argv1 hapmap.dat 
##argv2 resource dir
args <- commandArgs(trailingOnly = TRUE)
r <- mean.impute.missing( as.matrix(read.table(args[1])) )
geno.cent <- r$m - matrix(r$mu,nrow(r$m),ncol(r$m))

k <- 2
ir <- irlba(geno.cent,nu=k,nv=k)
U <- matrix(ir$u,ncol=k)
V <- matrix(ir$v,ncol=k)
d <- ir$d
UD <- U %*% diag(d)

write.table(U %*% diag(d), paste(args[1],".UD",sep=""),row.names=FALSE,col.names=FALSE)
write.table(V, paste(args[1],".V",sep=""),row.names=colnames(r$m), col.names=FALSE, quote=FALSE)
write.table(matrix(cbind(row.names(r$m),r$mu),nrow(r$m),2),paste(args[1],".mu",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)

