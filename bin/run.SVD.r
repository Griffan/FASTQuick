mean.impute.missing <- function(m) {
    mu <- rowMeans(m,na.rm=T)
    isna <- is.na(m)
    m[isna] <- matrix(mu,nrow(m),ncol(m))[isna]
    return(list(m=m,mu=mu))
}

## pcs : coordinates in PC space
pc.lik <- function(pcs, gls, UD, mu, n) {
    afs <- (UD %*% pcs + mu)/2
    afs[afs < 0.5/n] <- 0.5/n
    afs[afs > (n-0.5)/n] <- (n-0.5)/n
    gfs <- cbind((1-afs)*(1-afs), 2*afs*(1-afs), afs*afs)
    llk <- sum(log(rowSums(gls * gfs)))
    return(0-llk)
}

library(irlba)
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

n <- ncol(r$m)
GL <- as.matrix(read.table(paste(args[2],'/1kg.phase1.selected.GLs.dat',sep=""))) ##  M * (3Nk)
l <- ncol(GL)/3
ids <- substr(colnames(GL)[(1:l)*3],1,7)
pops <- substr(colnames(GL)[(1:l)*3],9,11)

out <- matrix(NA,l,k)

for(i in 1:l) {
   print(paste(ids[i],pops[i]))
  randstart <- V[as.integer(runif(1,0,n))+1,]
#    #randstart <- V[1,]    
    res <- optim( randstart, pc.lik, gls=10^GL[,3*(i-1)+1:3], UD=UD, mu=r$mu, n=n, method="Nelder-Mead")
    out[i,] <- res$par
#                                        #print(paste(ids[i],pops[i],res$par))
#    
}
write.table(out,paste((argv[2],'/1kg.phase1.selected.pc2.out',sep=""),row.names=ids,col.names=F,quote=F))
