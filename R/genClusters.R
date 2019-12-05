# Generate a dataset consists of two predictor variables and a class variable
genClusters <- function(k=2, n=100, mu=0, sd=5, seed=NULL, plot=TRUE){
   if(!is.null(seed)) set.seed(seed)
   x1 <- matrix(rnorm(2*n),n,2)
   samp <- sample(1:k, n, replace=TRUE)
   x2 <- matrix(rnorm(2*k, mean=mu, sd=sd), k, 2)
   dset <- cbind(x1+x2[samp,], samp)
   colnames(dset) <- c("V1","V2","Class")
   if(plot) plot(dset[,1],dset[,2], col=dset[,3], xlab="V1", ylab="V2", pch=19)
   return(dset)
}
