

library(CVXR)

library(glmnet)
lasso_reg <- function(beta, lambda = 0) {
  # lasso <- p_norm(beta[-1], 1)
  lasso <- sum(abs(beta[-1]))
  lambda * lasso
}





n=200
K=5
beta.true=c(1,0,1,0,1)
Xmi=matrix(rnorm((K-1)*n),n,(K-1))
X=cbind(1,Xmi)
y=X %*% beta.true
beta <- Variable(K)

lambda=0.1

loss <- sum((y - X %*% beta)^2)/(2*n)

obj <- loss + lasso_reg(beta, lambda=lambda)
prob <- Problem(Minimize(obj))
result <- solve(prob)

result$status
beta.cvxr=round(result$getValue(beta[-1]),5)

plot(result$getValue(beta))


glmnet.mod=glmnet(Xmi,y,lambda = lambda)
data.frame(true=beta.true[-1],cvxr=beta.cvxr,glmnet=as.numeric(glmnet.mod$beta))





library(mvabund)
data(spider)
Y=spider$abund
P=dim(Y)[2]
N=dim(Y)[1]
K=10

# Y=matrix(seq(1:(N*P)),N,P)
# mu=matrix(seq(1:(P*K)),P,K)
pi=runif(K)
pi=pi/sum(pi)


loss_fun=function(beta,Y,pi){
  N=dim(Y)[1]
  K=beta@cols
  outer=0
  for(i in 1:N){
    inner=0
    for(k in 1:K){
      inner=inner+pi[k]*sum((t(Y[i,])-beta[,k])^2)
    }
    outer=outer+log(inner)
  }
  outer
}




  
lasso_reg <- function(beta, lambda = 0) {
  K=beta@cols
  P=beta@rows
  outer=0
  for(k in 1:(K-1)){
    inner=0
    for(p in 1:P){
      inner=inner+(beta[p,k]-beta[p,(k+1)])^2
    }
    outer=outer+inner
  }
  lambda * outer
}

value(beta)<-matrix(seq(1:(P*K)),P,K)
# lasso_reg(beta,0.1)
beta <- Variable(P,K)
loss=loss_fun(beta,Y,pi)


# lambda_1 * p_norm(diff(x = beta, differences = 2), 1))

obj <- loss +lasso_reg(beta, lambda=lambda)
prob <- Problem(Minimize(obj))
result <- solve(prob,ignore_dcp = TRUE)

result$status

result$value

result$getValue(result)




#1d

library(mixtools)
library(CVXR)

mu=c(2,4,-1)
s=c(1,1,1)
Y <- c(rnorm(40, mean = mu[1], sd = s[1]),
       rnorm(20, mean = mu[2], sd = s[2]),
       rnorm(20, mean = mu[3], sd = s[3]))
hist(Y,20)
N=length(Y)

#initialise
#number of clusters
K=N/10

#k_mens for initial centres
init_clust <- kmeans(Y, K)
outEM=normalmixEM(Y, mu = init_clust$centers)

tau=outEM$posterior
beta <- Variable(K,name = "beta")


#As far as I understand the 

loss_fun=function(beta,Y,tau,sig){
  N=length(Y)
  out=0
  for(i in 1:N){
      out=out+sum(tau[i,]*(Y[i]-beta)^2*(1/(2*(sig^2))))
      #funny specification is so it will work with cvxr, e.g. dnorm will not
  }
  out
}

#so if I plug the loss function into optim, it shouldn't move the betas
#but it does

optim_res=optim(out$mu,loss_fun,Y=Y,tau=tau,sig=out$sigma)

cbind(sort(optim_res$par),sort(out$mu))

obj <- loss_fun(beta,Y,tau,as.numeric(out$sigma))
prob <- Problem(Minimize(obj))
result <- solve(prob)

cbind(sort(result$getValue(beta)),sort(optim_res$par),sort(out$mu))













lasso_reg <- function(beta, lambda = 0,w=NULL) {
  # outer=p_norm(diff(x = beta, differences = 1), 1)
  K=beta@rows
  out=0
  if(is.null(w)){
    w=matrix(1,K,K)
  }
  for(k in 1:(K-1)){
    out=out+p_norm(w[k,(k+1):K]*(beta[k]-beta[(k+1):K]),1)
  }
  lambda * out
}


# value(lasso_reg(beta,lambda=1))
# sum(diff(value(beta)))
# 
# 
# value(lasso_reg(beta,0.1))


lambdas=exp((seq(0,15,length.out=6))/6)-1
J=length(lambdas)
beta <- Variable(K)
loss=loss_fun(beta,Y,tau,as.numeric(out$sigma))
beta_out=matrix(NA,K,J)
w=exp(-dist(init_clust$centers))
w=as.matrix(w)

value(beta)=out$mu
value(loss_fun(beta,Y,tau,as.numeric(out$sigma)))

for(j in 1){
  # obj <- loss +lasso_reg(beta, lambda=lambdas[j],w)
  obj <- loss_fun(beta,Y,tau,as.numeric(out$sigma))
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  result$status
  result$value
  beta_out[,j]=result$getValue(beta)
}


cbind(out$mu,beta_out)
plot(lambdas+1,beta_out[1,],type="l",ylim=c(min(mu),max(mu)),log="x")
for(j in 2:J){
  lines(lambdas+1,beta_out[j,],type="l",col=j)
}
abline(h=mu,lty=2,col="red")



