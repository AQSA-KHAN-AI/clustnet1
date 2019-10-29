#https://cvxr.rbind.io/cvxr_examples/cvxr_l1-trend-filtering/

rm(list=ls())
library(mixtools)
library(CVXR)
library(l1tf)
library(tidyr)
library(ggplot2)

# 
# #example for fused lasson for smoothing gaussian
# ## lambda = 50
# y <- sp500$log
# lambda_1 <- 50 
# beta <- Variable(length(y))
# objective <- Minimize(0.5 * p_norm(y - beta) +
#                           lambda_1 * p_norm(diff(x = beta, differences = 2), 1))
# p1 <- Problem(objective)
# betaHat <- solve(p1)$getValue(beta)
# 
# 
# #example for fused lasson clustering gaussian-add pairwise penalty
# 
# penalty=function(beta,lambda){
#   out=0
#   n=beta@rows  ### need this
#   for(i in 1:n){
#     for(j in i:n)
#       out=out+abs(beta[i]-beta[j])
#   }
#   lambda * out
# }
# 
# 
# n_lam=20
# lambdas=exp(seq(log(0.008),log(0.04),length.out = n_lam))
# K=10
# y <- c(rnorm(K,1,1),rnorm(K,10,1))
# betasHat=matrix(NA,2*K,n_lam)
# for(i in 1:length(lambdas))
# {
#   lambda <- lambdas[i]
#   beta <- Variable(length(y))
#   obj <- 0.5* p_norm(y - beta) + penalty(beta,lambda)
#   problem<- Problem(Minimize(obj))
#   betasHat[,i] <- solve(problem)$getValue(beta)
# }
# 
# # save(betasHat,file="gaussian_2 grps")
# load(file="gaussian_2_grps")
# minmax=range(betasHat)
# plot(lambdas,betasHat[1,],type="l",
#      ylim=c(minmax[1]-1,minmax[2]+1),lwd=2,
#      ylab="Group means",xlab = "lambda")
# for(j in 2:(2*K)){
#   lines(lambdas,betasHat[j,],type="l",col=j,lwd=2)
# }
# 
# 
# 
# #example for fused lasson clustering poisson 1d
# #beta = log lambda(from poisson)
# 
# 
# n_lam=10
# lambdas=exp(seq(log(0.01),log(0.6),length.out = n_lam))
# K=50
# y <- cbind(c(rpois(K,1),rpois(K,10)),c(rpois(K,3),rpois(K,5)))
# betasHat=array(NA,c(2*K,n_lam,2))
# for(i in 1:length(lambdas))
# {
#   lambda <-lambdas[i]
#   beta <- Variable(nrow(y),ncol(y))
#   obj <-  sum_entries(-y * beta + exp(beta)) + penalty(beta,lambda)
#   problem<- Problem(Minimize(obj))
#   A=solve(problem)
#   betasHat[,i] <- A$getValue(beta)
# }
# # save(betasHat,file="poisson_2_grps")
# # plot(y,exp(A$getValue(beta)))
# # abline(0,1)
# 
# 
# load(file="poisson_2_grps")
# 
# minmax=range(betasHat)
# plot(lambdas,betasHat[1,],type="l",
#      ylim=c(minmax[1]-1,minmax[2]+1),lwd=2,
#      ylab="Group means",xlab = "lambda")
# for(j in 2:(2*K)){
#   lines(lambdas,betasHat[j,],type="l",col=j,lwd=2)
# }
# 
# 
# 
# 
# #example bivariate poisson fused lasson clustering
# #beta = log lambda(from poisson)
# 
penalty=function(beta,lambda){
  out=0
  n=beta@rows  ### need this
  for(i in 1:n){
    for(j in i:n)
      out=out+norm(beta[i,]-beta[j,],"F")
  }
  lambda * out
}
# 
# 
# 
# 
n_lam=10
lambdas=exp(seq(log(0.01),log(0.6),length.out = n_lam))
K=10
load("y_bivariate")
# y <- cbind(c(rpois(K,1),rpois(K,10)),c(rpois(K,3),rpois(K,5)))
# betasHat=array(NA,c(n_lam,2*K,2))
# for(i in 1:length(lambdas))
# {
#   lambda <-lambdas[i]
#   beta <- Variable(nrow(y),ncol(y))
#   obj <-  sum_entries(-y * beta + exp(beta)) + penalty(beta,lambda)
#   problem<- Problem(Minimize(obj))
#   A=solve(problem)
#   betasHat[i,,] <- A$getValue(beta)
# }
# 
# 
# save(betasHat,file="poisson_bivariate_2_grps")

load(file="poisson_bivariate_2_grps")

dat1=data.frame(betasHat[,,1],lambda=lambdas[])
colnames(dat1)[1:(2*K)]=1:(2*K)
plot_dat=gather(dat1,obs,lmu1,1:(2*K))
plot_dat$lmu2=gather(data.frame(betasHat[,,2]))$value
plot_dat$mu1=exp(plot_dat$lmu1)
plot_dat$mu2=exp(plot_dat$lmu2)


Y=data.frame(y,obs=1:(2*K))
colnames(Y)[1:2]=c("mu1","mu2")
true_means=data.frame(mu1=c(1,10),mu2=c(3,5),obs="True mean")
ggplot(plot_dat,aes(mu1,mu2,color=obs))+geom_path()+theme_bw()+  geom_point()+
  geom_point(data=true_means,shape=15,size=5)+
  theme(legend.position = "none")+
  ggtitle("poisson loss")+
  ylim(c(2.5,7))+xlim(c(1,14))
# 
# save(y,file="y_bivariate")
# 
#least squares for comparison

# n_lam=10
# lambdas=exp(seq(log(0.02),log(1.2),length.out = n_lam))
# load("y_bivariate")
# betasHat=array(NA,c(n_lam,2*K,2))
# for(i in 1:length(lambdas))
# {
#   lambda <-lambdas[i]
#   beta <- Variable(nrow(y),ncol(y))
#   obj <-  sum_entries((y-beta)^2) + penalty(beta,lambda)
#   problem<- Problem(Minimize(obj))
#   A=solve(problem)
#   betasHat[i,,] <- A$getValue(beta)
# }
# 
# 
# save(betasHat,file="poisson_bivariate_2_grps_least_squares")




load(file="poisson_bivariate_2_grps_least_squares")

dat1=data.frame(betasHat[,,1],lambda=lambdas)
colnames(dat1)[1:(2*K)]=1:(2*K)
plot_dat=gather(dat1,obs,mu1,1:(2*K))
plot_dat$mu2=gather(data.frame(betasHat[,,2]))$value
Y=data.frame(y,obs=1:(2*K))
colnames(Y)[1:2]=c("mu1","mu2")
true_means=data.frame(mu1=c(1,10),mu2=c(3,5),obs="True mean")
ggplot(plot_dat,aes(mu1,mu2,color=obs))+geom_path()+theme_bw()+
  geom_point(aes(mu1,mu2,color=factor(lambda)))+
  geom_point(data=true_means,shape=15,size=5)+
  theme(legend.position = "none")+
  ggtitle("least squares loss")+
  ylim(c(2.5,7))+xlim(c(1,14))
