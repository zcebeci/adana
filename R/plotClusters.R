plotClusters <- function(x, idx=c(1,2,3)){
   plot(x[,idx[1]], x[,idx[2]], col=x[,idx[3]], xlab="V1", ylab="V2", pch=19)
}
