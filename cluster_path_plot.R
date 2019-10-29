rm(list=ls())
load("cluster_path.RData")
minmax=range(beta_out)
plot(lambdas,beta_out[1,],type="l",
     ylim=c(minmax[1]-1,minmax[2]+1),log="x",lwd=2,
     ylab="Group means",xlab = "lambda")
for(j in 2:9){
  lines(lambdas,beta_out[j,],type="l",col=j,lwd=2)
}