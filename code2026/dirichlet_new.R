source("sampling.R")
library(gtools)

Gibbs = function(data, a0=0, b0=0, N=1000, N_truncations=300, Burn=200)
{
  S = dim(data)[1]
  K = dim(data)[2]
  Ns = apply(data,1,sum)
  
  tau = matrix(0, S, K)
  omega = matrix(0, S, K)
  eta = rep(0, S)
  alpha = rep(1/K, K)
  alpha_samples = matrix(0,N,K)
  
  b_pre = S*(-digamma(1))-b0
  
  for(i in 1:(N+Burn))
  {
    tau_para = data + matrix(rep(alpha,each=S),S,K)
    tau = pmax(matrix(rgamma(S*K, shape = tau_para), S, K),1e-200)
    alpha0 = sum(alpha)
    eta = sapply(Ns, function(x) rbeta(1,shape1=alpha0, shape2=x))
    
    for(k in 1:K)
    {
      omega[,k] = rPIG_ERGamma(S, alpha[k], N_truncations)
      a = sum(omega[,k])
      b = b_pre + sum(log(eta)) + sum(log(tau[,k]))
      alpha[k] = rH(S+a0, a, b)
    }
      
    if(i > Burn)
      alpha_samples[i-Burn,] = alpha
  }
  
  return(alpha_samples)
}

Gibbs_CLT = function(data, a0=0, b0=0, N=1000, N_truncations=300, Burn=200)
{
  S = dim(data)[1]
  K = dim(data)[2]
  Ns = apply(data,1,sum)
  
  alpha = rep(1/K, K)
  alpha_samples = matrix(0,N,K)
  
  b_pre = S*(-digamma(1))-b0
  
  for(i in 1:(N+Burn))
  {
    alpha0 = sum(alpha)
    mu_a = S*(digamma(1+alpha)-digamma(1))/(2*alpha)
    sigma_a = mu_a/(2*alpha^2) - S*trigamma(1+alpha)/(4*alpha^2)
    sigma_a = sqrt(sigma_a)
    
    tau_para = data + matrix(rep(alpha,each=S),S,K)
    mu_b = b_pre + S*digamma(alpha0) + apply(digamma(tau_para),2,sum) - sum(digamma(alpha0+Ns))
    sigma_b = S*trigamma(alpha0) + apply(trigamma(tau_para),2,sum) - sum(trigamma(alpha0+Ns))
    sigma_b = sqrt(sigma_b)
    
    a = rnorm(K)*sigma_a + mu_a
    b = rnorm(K)*sigma_b + mu_b
    
    for(k in 1:K)
      alpha[k] = rH(S+a0, a[k], b[k])
    
    if(i > Burn)
      alpha_samples[i-Burn,] = alpha
  }
  
  return(alpha_samples)
}

Gibbs_E = function(data, a0=0, b0=0, N=1000, N_truncations=300, Burn=200)
{
  S = dim(data)[1]
  K = dim(data)[2]
  Ns = apply(data,1,sum)
  
  alpha = rep(1/K, K)
  alpha_samples = matrix(0,N,K)
  
  b_pre = S*(-digamma(1))-b0
  
  for(i in 1:(N+Burn))
  {
    alpha0 = sum(alpha)
    mu_a = S*(digamma(1+alpha)-digamma(1))/(2*alpha)
    
    tau_para = data + matrix(rep(alpha,each=S),S,K)
    mu_b = b_pre + S*digamma(alpha0) + apply(digamma(tau_para),2,sum) - sum(digamma(alpha0+Ns))
    
    for(k in 1:K)
      alpha[k] = rH(S+a0, mu_a[k], mu_b[k])
    
    if(i > Burn)
      alpha_samples[i-Burn,] = alpha
  }
  
  return(alpha_samples)
}

ltarget = function(data,alpha,a0=0, b0=0)
{

  res = sum(apply(data, 1, function(x){
    lgamma(sum(alpha))-lgamma(sum(alpha+x)) + sum(lgamma(alpha+x)-lgamma(alpha))
  }))
  
  res = res + (a0-1)*sum(log(alpha)) - b0*sum(alpha)
  return(res)
}

MH = function(data,ltarget,a0=0,b0=0,step,init=NULL,N=1000, Burn=200)
{
  K = dim(data)[2]
  x_samples = matrix(0, N, K)
  if(is.null(init))
  {
    x = rep(1/K,K)
  }else{
    x = init
  }
  
  for(i in 1:(N+Burn))
  {
    propose = x*exp(rnorm(K, sd=step))
    ratio = exp(ltarget(data,propose,a0,b0)-ltarget(data,x,a0,b0))
    if(runif(1)<ratio)
    {
      x = propose
    }
    
    if(i >Burn)
    {
      x_samples[i-Burn,] = x
    }
  }
  return(x_samples)
}

# 
# set.seed(123)
# K = 5 # number of categories
# alpha_true = 1:K
# M = 10
# data = matrix(0, M,K)
# p_data = matrix(0, M,K)
# for(i in 1:M)
# {
#   p_true = rdirichlet(1, alpha_true)
#   data[i,] = rmultinom(1, rpois(1,500), prob = p_true)
#   p_data[i,] = p_true
# }
# 
# gibbs_alpha = Gibbs(data)
# gibbs_alpha_c = Gibbs_CLT(data)
# gibbs_alpha_e = Gibbs_E(data)
# MH_alpha = MH(data,ltarget, step=0.1)
