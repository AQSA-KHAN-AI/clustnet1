##soft threshold function
soft_thr=function(beta,lambda){
  out=0
  if(beta>lambda) out=beta-lambda
  if(beta<(-lambda)) out=beta+lambda
  out
}


lik.bin=function(y,X,beta){
  mu=X%*%beta
  logl=sum(dbinom(y,1,inv.logit(mu),log = TRUE))
  -2*logl
}




llasso<-function(lambda,y,Xmat,eps=1e-7){
  K=dim(X)[2]-1
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
        if(j!=1){
          beta_loc_current[j]=soft_thr(sum(X[,j]*w*rj)/N,lambda)/(sum(X[,j]*w*X[,j])/N)
        }
        
      }
      #is there convergence on the inner loop
      diff2=abs(lik.bin(y,X,beta_loc_current)-lik.bin(y,X,beta_old))
      beta_old=beta_loc_current
    }
    #is there convergence of outer loop?
    diff=abs(lik.bin(y,X,beta_current)-lik.bin(y,X,beta_loc_current))
    
    # print(diff)
    beta_current=beta_loc_current
  }
  
  beta_current
}



# ##coordinate descent
# #initialize
# eps=1e-7
# diff=eps+1
# beta_current=matrix(c(logit(mean(y)),rep(0,K)))
# 
# while(diff>eps){
#   
#   #do the quadratic approximation on the current beta estimate
#   mu=X%*%beta_current
#   p=inv.logit(mu)
#   w=p*(1-p)
#   z=mu+(y-p)/(p*(1-p))
#   
#   #initialize a local beta for coordinate descent
#   beta_loc_current=beta_current
#   
#   #initialize beta_old for convergence check
#   beta_old=beta_loc_current
#   diff2=eps+1
#   
#   #do coordinate descent till convegence of weighted least squares problem
#   #exclude intercept
#   while(diff2>(eps)){
#     for(j in 1:(K+1)){
#       rj=z-X[,-c(j)]%*%beta_loc_current[-c(j)]
#       beta_loc_current[j]=sum(X[,j]*w*rj)/sum(X[,j]*w*X[,j])
#       
#     }
#     print(cbind(wls,beta_loc_current))
#     #is there convergence on the inner loop
#     diff2=abs(lik.bin(y,X,beta_loc_current)-lik.bin(y,X,beta_old))
#     beta_old=beta_loc_current
#   }
#   #is there convergence of outer loop?
#   diff=abs(lik.bin(y,X,beta_current)-lik.bin(y,X,beta_loc_current))
#   
#   # print(diff)
#   beta_current=beta_loc_current
# }
# 
# cbind(beta_true,beta_current,wls,round(beta_current-wls,4))


#EM functions
logl_mix<-function(y,pi,mu,sig,phi){
  inner=apply(f_mix(y,mu,sig,phi)*pi[,rep(1,N)],2,sum)
  sum(log(inner))
}

gl_alg<-function(tau,y){
  obj=lm(y~1,weights =tau )
  c(coef(obj),sigma(obj))
}



f_mix<-function(y,mu,sig,phi=NA){
  K=length(mu)
  N=length(y)
  apply(matrix(y),1,dnorm,mean = mu,sd=sig)
}
