
library(boot)
library(mvtnorm)
library(glmnet)
library(bazar)
setwd("C:/Users/Gordana/Dropbox/Projects/Cluster Lasso/Code")
source("functions.R")
N=100
K=6
#simulate data
x=scale(rmvnorm(N,mean = rep(0,K),sigma = diag(K)))
X=cbind(1,x)
beta_true=matrix(c(1,runif(K)))
mu_true=X%*%beta_true
y=rbinom(N,1,inv.logit(mu_true))

#estimate with glm
use_glm=glm(y~x,family = binomial)

#weighted least squares from scratch




eps=1e-7
diff=eps+1

#initialize
beta_current=matrix(c(logit(mean(y)),rep(0,K)))

while(diff>eps){
  mu=X%*%beta_current
  p=inv.logit(mu)
  w=p*(1-p)
  z=mu+(y-p)/(p*(1-p))
  W=diag(as.numeric(w))
  beta_new=qr.solve(t(X)%*%W%*%X, t(X)%*%W%*%z)
  diff=abs(lik.bin(y,X,beta_new)-lik.bin(y,X,beta_current))
  beta_current=beta_new
}

cbind(beta_true,beta_current,use_glm$coefficients,round(beta_current-use_glm$coefficients,4))

wls=beta_current

#RMSE
sqrt(mean((predict(use_glm)-y)^2))

##coordinate descent
#initialize
eps=1e-7
diff=eps+1
beta_current=matrix(c(logit(mean(y)),rep(0,K)))

while(diff>eps){
  
  #do the quadratic approximation on the current beta estimate
  mu=X%*%beta_current
  p=inv.logit(mu)
  w=p*(1-p)
  z=mu+(y-p)/(p*(1-p))
  
  #initialize a local beta for coordinate descent
  beta_loc_current=beta_current

  #initialize beta_old for convergence check
  beta_old=beta_loc_current
  diff2=eps+1
  
  #do coordinate descent till convegence of weighted least squares problem
  #exclude intercept
  while(diff2>(eps)){
    for(j in 1:(K+1)){
      rj=z-X[,-c(j)]%*%beta_loc_current[-c(j)]
      beta_loc_current[j]=sum(X[,j]*w*rj)/sum(X[,j]*w*X[,j])
      
    }
    print(cbind(wls,beta_loc_current))
    #is there convergence on the inner loop
    diff2=abs(lik.bin(y,X,beta_loc_current)-lik.bin(y,X,beta_old))
    beta_old=beta_loc_current
  }
  #is there convergence of outer loop?
  diff=abs(lik.bin(y,X,beta_current)-lik.bin(y,X,beta_loc_current))
  
  # print(diff)
  beta_current=beta_loc_current
}

cbind(beta_true,beta_current,wls,round(beta_current-wls,4))


glmnet_est=glmnet(X[,-1],y,family="binomial",standardize=FALSE)
lambdas=glmnet_est$lambda
my_est=apply(as.matrix(lambdas),1,llasso,y=y,Xmat=X)

par(mfrow=c(2,1))
plot(lambdas,my_est[2,],type = "l",ylim=c(-0.3,1.2),log="x")
for(i in 3:7){
  lines(lambdas,my_est[i,],col=i-1,log="x")
  
}
plot(glmnet_est,xvar="lambda",ylim=c(-0.3,1.2))
#yup, works




#em algorithm for mixture modes 1d
#simulate

mu=c(2,4,-1)
s=c(0.6,1,0.3)
y <- c(rnorm(100, mean = mu[1], sd = s[1]),
       rnorm(50, mean = mu[2], sd = s[2]),
       rnorm(50, mean = mu[3], sd = s[3]))
hist(y,10)
N=length(y)

#initialise
#number of clusters
K=3
init_clust <- kmeans(y, K)

mu_init=init_clust$centers+rnorm(3)/2
mu_t=mu_init
# mu_t=rep(0,K)

sig_t=rep(1,K)#use kmens for this
pi_t=matrix(init_clust$size,K,1)
pi_t=pi_t/sum(pi_t)
phi_t=NA #overdispersion, forget for now



eps=1e-7
diff=eps+1

while(diff>eps){
  #E-Step (eq 3.9 Francis thesis)
  
  tau_t=pi_t[,rep(1,N)]*f_mix(y,mu_t,sig_t,phi_t)
  sum_tau=matrix(apply(tau_t,2,sum),1,N)
  sum_tau=sum_tau[rep(1,K),]
  tau_t=tau_t/sum_tau
  
  #M-Step
  pi_tp1=matrix(apply(tau_t,1,mean),K,1)
  lm_tp1=apply(tau_t,1,gl_alg,y=y)
  mu_tp1=matrix(lm_tp1[1,],K,1)
  sig_tp1=matrix(lm_tp1[2,],K,1)
  
  
  #check convergence
  diff=logl_mix(y,pi_tp1,mu_tp1,sig_tp1,phi)-logl_mix(y,pi_t,mu_t,sig_t,phi)
  pi_t=pi_tp1
  mu_t=mu_tp1
  sig_t=sig_tp1
}


#sigmas are wierd?

cbind(sort(mu),sort(mu_t),sort(mu_init))

cbind(sort(s),sort(sig_t),sort(sqrt(sig_t)))




# install.packages("cvxclustr")
library(cvxclustr)


## Clusterpaths for Mammal Dentition
data(mammals)
X <- as.matrix(mammals[,-1])
X <- t(scale(X,center=TRUE,scale=FALSE))
n <- ncol(X)
## Pick some weights and a sequence of regularization parameters.
k <- 5
phi <- 0.5
w <- kernel_weights(X,phi)
w <- knn_weights(w,k,n)
w[]=1


library(mvabund)
library(mvtnorm)
data(spider)
# X=t(log(as.matrix(spider$abund)+1))
# X=t(matrix(rnorm(prod(dim(spider$abund))),dim(spider$abund)))
X=rmvnorm(12,rep(0,28),cor(t(spider$abund)))
nlambda=20
gamma <- seq(0.0,.25, length.out=nlambda)
## Perform clustering
sol <- cvxclust(X,w,gamma,type=2)




U=array(unlist(sol$U),c(dim(sol$U[[1]]),nlambda))


#plot number of unique "betas

nunique=matrix(NA,12,nlambda)
nclust=rep(NA,nlambda)
for(nl in 1:nlambda){
  A <- create_adjacency(sol$V[[nl]],w,n)
  nclust[nl]=length(unique(find_clusters(A)$cluster))
  for(j in 1:12){
    nunique[j,nl]=length(almost.unique(U[j,,nl],tolerance = 0.01))
  }
}

nunique

plot(nunique[1,],ylim=c(min(nunique),max(nunique)),type="l")
for (i in 2:12) {
  lines(nunique[i,],col=i) 
}
lines(nclust,lwd=3)




plot(U[j,1,],ylim=c(min(U[j,,]),max(U[j,,])),type="l")
for (i in 2:28) {
  lines(U[j,i,],col=i) 
}



plot(U[j,1,],ylim=c(min(U[j,,]),max(U[j,,])),type="l")
for (i in 2:28) {
  lines(U[j,i,],col=i) 
}








nGamma <- sol$nGamma
beta1=matrix(NA,28,nGamma)

for (j in 1:nGamma) {
  beta1[,j] <-sol$U[[j]][3,]
}

plot(beta1[1,],ylim=c(min(beta1),max(beta1)),type="l")
for (i in 2:27) {
  lines(beta1[i,],col=i) 
}


## Output Cluster Assignment at 10th gamma
A <- create_adjacency(sol$V[[60]],w,n,method='admm')
find_clusters(A)
## Visualize Cluster Assignment
G <- graph.adjacency(A, mode = 'upper')
plot(G,vertex.label=as.character(colnames(spider$abund)),vertex.label.cex=0.65,vertex.label.font=2)


n=28
## Plot the cluster path
library(ggplot2)
svdX <- svd(X)
pc <- svdX$u[,1:2,drop=FALSE]
pc.df <- as.data.frame(t(pc)%*%X)
nGamma <- sol$nGamma
df.paths <- data.frame(x=c(),y=c(), group=c())
for (j in 1:nGamma) {
  pcs <- t(pc)%*%sol$U[[j]]
  x <- pcs[1,]
  y <- pcs[2,]
  df <- data.frame(x=pcs[1,], y=pcs[2,], group=1:n)
  df.paths <- rbind(df.paths,df)
}
X_data <- as.data.frame(t(X)%*%pc)
colnames(X_data) <- c("x","y")
X_data$Name <- 1:28
data_plot <- ggplot(data=df.paths,aes(x=x,y=y))
data_plot <- data_plot + geom_path(aes(group=group),colour='grey30',alpha=0.5)
data_plot <- data_plot + geom_text(data=X_data,aes(x=x,y=y,label=Name),
                                   position=position_jitter(h=0.125,w=0.125))
data_plot <- data_plot + geom_point(data=X_data,aes(x=x,y=y),size=1.5)
data_plot <- data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
data_plot + theme_bw()



