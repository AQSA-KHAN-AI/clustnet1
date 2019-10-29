
rm(list=ls())
library(mixtools)
library(CVXR)
library(l1tf)
library(tidyr)
library(ggplot2)
source("functions_clust.R")
# set.seed(1)
# ni=5
# mus=matrix(c(2,6,6,10,6,2,10,6),nrow = 4)
# Y=simulate_biv(rep(ni,nrow(mus)), mus)

# 
# dat=list(Y=Y,mus=mus)
# save(dat,file="poisson_bivariate_4_grps")

plot(jitter(dat$Y))
points(dat$mus,col="blue",pch=16)


load(file="poisson_bivariate_4_grps")
y=dat$Y
n_lam=10

lambdas=exp(seq(log(0.02),log(1.2),length.out = n_lam))
betasHat=array(NA,c(n_lam,nrow(y),2))
for(i in 1:length(lambdas))
{
    lambda <-lambdas[i]
    beta <- Variable(nrow(y),ncol(y))
    obj <-  sum_entries(-y * beta + exp(beta)) + penalty_ebeta(beta,lambda)
    problem<- Problem(Minimize(obj))
    A=solve(problem)
    betasHat[i,,] <- A$getValue(beta)
}


save(betasHat,file="poisson_bivariate_2_grps_poiss_ebeta")

load(file="poisson_bivariate_2_grps_poiss_ebeta")




K=10
dat1=data.frame(betasHat[,,1],lambda=lambdas)
colnames(dat1)[1:(2*K)]=1:(2*K)
plot_dat=gather(dat1,obs,mu1,1:(2*K))
plot_dat$mu2=gather(data.frame(betasHat[,,2]))$value
plot_dat$mu1=exp(plot_dat$mu1)
plot_dat$mu2=exp(plot_dat$mu2)
Y=data.frame(y,obs=factor(1:(2*K)))
colnames(Y)[1:2]=c("mu1","mu2")
true_means=data.frame(mus,obs="True mean")
colnames(true_means)[1:2]=c("mu1","mu2")
ggplot(plot_dat,aes(mu1,mu2,group=obs))+geom_path()+theme_bw()+
  geom_point(aes(mu1,mu2,color=factor(lambda)))+
  geom_point(data=true_means,shape=15,size=5)+
  geom_point(data=Y)+
  theme(legend.position = "none")+
  ggtitle("poisson lik loss")
