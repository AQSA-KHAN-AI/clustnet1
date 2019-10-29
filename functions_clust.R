


simulate_biv<-function(Ns,means, poisson=TRUE){
  out=NULL
  for(i in 1:length(Ns)){
    this_out=cbind(rpois(Ns[i],means[i,1]),rpois(Ns[i],means[i,2]))
    out=rbind(out,this_out)
  }
  # out[sample(1:nrow(out)),]
  out
}


penalty_beta=function(beta,lambda){
  out=0
  n=beta@rows  ### need this
  for(i in 1:(n-1)){
    for(j in (i+1):n)
      out=out+norm(beta[i,]-beta[j,],"F")
  }
  lambda * out
}

