library(prclust)
library(tidyr)
library(mvabund)
library(ggplot2)
library(dplyr)
# data(spider)
# 
# n_lam=100
# lambdas=exp(seq(log(0.001),log(2000),length.out = n_lam))
# Y=log(t(spider$abund)+1)
# betasHat=array(NA,c(n_lam,nrow(Y),ncol(Y)))
# for(i in 1:length(lambdas))
# {
#   betasHat[i,,] <- PRclust(Y,lambda1=1,lambda2=lambdas[i],tau=4,loss.method="lasso",algorithm = "ADMM")$mu
# }
# 
# dat1=data.frame(betasHat[,,1])
# plot_dat=gather(dat1,obs,mu2)
# plot_dat$mu5=gather(data.frame(betasHat[,,8]))$value
# ggplot(plot_dat,aes(mu2,mu5,color=obs))+geom_path()+theme_bw()+geom_point()+
#   theme(legend.position = "none")+
#   ggtitle("PRclust_spider")


# 

load(file="poisson_bivariate_4_grps")
# 
# set.seed(2)
# ni=20
# mus=matrix(c(2,6,6,10,6,2,10,6),nrow = 4)
# Y=simulate_biv(rep(ni,nrow(mus)), mus)
# dat=list(Y=Y,mus=mus)

J=matrix(rnorm(prod(dim(dat$Y))),nrow = nrow(dat$Y))/100
Y=t(log((dat$Y+1+J)))
n_lam=1000
lambdas=exp(seq(log(0.1),log(100),length.out = n_lam))

betasHat=array(NA,c(n_lam,nrow(Y),ncol(Y)))
for(i in 1:length(lambdas))
{
  betasHat[i,,] <- PRclust(Y,lambda1=0.5,lambda2=lambdas[i],tau=0.7,loss.method="lasso",algorithm = "ADMM")$mu
}

dat1=data.frame(betasHat[,1,])
plot(dat1[,4])
plot_dat=gather(dat1,obs,mu1)
plot_dat$mu1=exp(plot_dat$mu1)-1
plot_dat$mu2=exp(gather(data.frame(betasHat[,2,]))$value)-1
plot_dat$obs=factor(rep(1:ncol(Y),each=length(lambdas)))
plot_dat$lambda=rep(lambdas,ncol(Y))
true_means=data.frame(dat$mus,obs="True mean")
colnames(true_means)[1:2]=c("mu1","mu2")
y=data.frame(t(exp(Y)-1))
y$obs=factor(1:nrow(y))
colnames(y)[1:2]=c("mu1","mu2")
ggplot(plot_dat,aes(mu1,mu2,group=obs))+theme_bw()+
  geom_path()+
  geom_point(aes(mu1,mu2,color=factor(lambda)))+
  geom_point(data=true_means,shape=15,size=5)+
  geom_point(data=y)+
  theme(legend.position = "none")+
  ggtitle("PRclust")

