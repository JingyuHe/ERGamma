library(gtools)
library(GIGrvg)
library(truncnorm)

rPIG = function(d, c, N) {
  # d is a vector
  # c is a scalar
  # N is number of truncation
  K = length(d)
  output = 0
  if (N <= K) {
    c_abs = abs(c)
    for (i in 1:N) {
      output = output + rgig(n = 1, lambda = -3 / 2, chi = 1 / 2 / d[i] ^ 2, c_abs ^ 2)
    }
  } else {
    print("N should not be bigger than length of d. \n")
  }
  return(output)
}

rPIG_ERGamma = function(cc, N) {
  # set dk = k
  # c is a scalar
  # N is number of truncation
  output = 0
  c_abs = abs(cc)
  for (i in 1:N) {
    output = output + rgig(n = 1, lambda = -3 / 2, chi = 1 / 2 / i ^ 2, c_abs ^ 2)
    # output = output + rgig(n = 1, lambda = c_abs, chi = 1 / 2 / i^2, psi = 9 / 4)
  }
  return(output)
}


# Gibbs sampler ---
Gibbs = function(data, tau, niter=1000, N_truncations=300, lambda=0)
{
  M = dim(data)[1]
  K = dim(data)[2]
  data_colsum = apply(data,2,sum)
  data_rowsum = apply(data,1,sum)
  w = array(0, dim = c(niter, M, K))
  alpha = matrix(1, niter, K)
  eta = matrix(1, niter, M)
  p = array(1 / K, dim = c(niter, M, K))
  b_vec = matrix(0, niter, K)
  a_vec = matrix(0, niter, K)
  
  for (i in 2:niter) {
    cat("sampling ", i, "\n")
    
    # sampling eta
    for (m in 1:M) {
      eta[i, m] = rgamma(1, shape = sum(alpha[i - 1,]) + data_rowsum[m])
    }
    
    
    # sampling w
    for (m in 1:M) {
      for (k in 1:K) {
        w[i, m, k] = rPIG_ERGamma(sqrt(2 * (alpha[i - 1, k] + data[m, k] - 1)^2), N_truncations)
      }
    }
    
    # sampling p
    for (m in 1:M) {
      p[i, m,] = pmax(rdirichlet(1, alpha[i - 1,] + data[m,]),1e-7) # prevent p=0 case
    }
    
    
    # sampling alpha
    a = apply(w[i,,], 2, sum) + ifelse(tau==0, 0, 1/(2*tau^2))
    b = -2 * apply((data - 1) * as.matrix(w[i,,]), 2, sum) + sum(log(eta[i,])) + M * 0.577216 + apply(log(p[i,,]), 2, sum) - lambda
    b_vec[i,] = b
    a_vec[i,] = a
    
    for (k in 1:K) {
      alpha[i,k] = rtruncnorm(1, a = 0, mean = b[k] / 2 / a[k], sd = sqrt(1 / 2 / a[k]))
    }
    
  }
  
  return(alpha)
}

Gibbs_E = function(data, tau, niter=1000, lambda=0)
{
  M = dim(data)[1]
  K = dim(data)[2]
  data_colsum = apply(data,2,sum)
  data_rowsum = apply(data,1,sum)
  w = array(0, dim = c(niter, M, K))
  alpha = matrix(0.5, niter, K)
  eta = matrix(1, niter, M)
  p = array(1 / K, dim = c(niter, M, K))
  b_vec = matrix(0, niter, K)
  a_vec = matrix(0, niter, K)
  
  for (i in 2:niter) {
    cat("sampling ", i, "\n")
    
    # sampling eta
    for (m in 1:M) {
      eta[i, m] = rgamma(1, shape = sum(alpha[i - 1,]) + data_rowsum[m])
    }
    
    
    # sampling w
    c_array = sqrt(2) * abs(data + rep(alpha[i-1,],each=M) - 1)
    w[i,,] = (0.577216 + digamma(c_array/sqrt(2) + 1))/(sqrt(2)*c_array)
    
    # sampling p
    for (m in 1:M) {
      p[i, m,] = pmax(rdirichlet(1, alpha[i - 1,] + data[m,]),1e-7) # prevent p=0 case
    }
    # p_e = data + rep(alpha[i-1,],each=M)
    # p[i,,] = p_e / apply(p_e,1,sum)
    
    # sampling alpha
    a = apply(w[i,,], 2, sum) + ifelse(tau==0, 0, 1/(2*tau^2))
    b = -2 * apply((data - 1) * as.matrix(w[i,,]), 2, sum) + sum(log(eta[i,])) + M * 0.577216 + apply(log(p[i,,]), 2, sum) - lambda
    # b = -2 * apply((data - 1) * as.matrix(w[i,,]), 2, sum) + sum(log(eta[i,])) + M * 0.577216 + apply(log(p), 2, sum) - lambda
    b_vec[i,] = b
    a_vec[i,] = a
    
    for (k in 1:K) {
      alpha[i,k] = rtruncnorm(1, a = 0, mean = b[k] / 2 / a[k], sd = sqrt(1 / 2 / a[k]))
    }
    
  }
  
  return(alpha)
}

# MH -------
ltarget = function(data,alpha,tau)
{
  res = sum(apply(data, 1, function(x){
    lgamma(sum(alpha))-lgamma(sum(alpha+x)) + sum(lgamma(alpha+x)-lgamma(alpha))
  }))
  res = res+sum(log(2)+log(dnorm(alpha,sd=tau)))
  return(res)
}

MH = function(data,ltarget,tau,step,init=NULL,niter=1000)
{
  K = dim(data)[2]
  x = matrix(0, niter, K)
  if(is.null(init))
  {
    x[1,] = abs(rnorm(K,sd=tau))
  }else{
    x[1,] = init
  }
  
  for(i in 2:niter)
  {
    propose = x[i-1,]*exp(rnorm(K, sd=step))
    ratio = exp(ltarget(data,propose,tau)-ltarget(data,x[i-1,],tau))
    if(runif(1)<ratio)
    {
      x[i,] = propose
    }else{
      x[i,] = x[i-1,]
    }
  }
  return(x)
}


# simulation ----

# data generating process
set.seed(123)
K = 5 # number of categories
alpha_true = 1:K
M = 10
data = matrix(0, M,K)
p_data = matrix(0, M,K)
for(i in 1:M)
{
  p_true = rdirichlet(1, alpha_true)
  data[i,] = rmultinom(1, rpois(1,500), prob = p_true)
  p_data[i,] = p_true
}

gibbs_alpha = Gibbs_E(data, tau = 5, niter=2000)
MH_alpha = MH(data,ltarget, tau = 5, step=0.1, niter=2000)

colMeans(gibbs_alpha[1001:2000,])
colMeans(MH_alpha[1001:2000,])




# opioid data ---
opioid_data = c(118, 96, 298, 58, 170, 18, 171, 240, 403, 72, 180, 24, 250, 527, 470, 67, 209, 23, 280, 682, 405, 97, 180, 39, 109, 237, 327, 150, 400, 17, 154, 386, 418, 179, 394, 21)
opioid_data = matrix(opioid_data, 6, 6, byrow = TRUE)



