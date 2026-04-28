library(GIGrvg)
library(gtools)
gamma_const = 0.577215664901532

# draw k samples from PIG(,cc)
rPIG_ERGamma = function(k,cc,N){
  c_abs = abs(cc)
  chi = 1 / 2 / (1:N) ^ 2
  output = matrix(sapply(chi, function(i) rgig(n = k, lambda = -3 / 2, chi = i, c_abs ^ 2)),
                  k,N)
  return(apply(output,1,sum))
}

# draw k samples from PIG(,0)
rPIG_noc = function(k,N=100){
  rg = 1/matrix(rgamma(k*N,shape = 3/2),N,k)
  rg = rg*(1/4/(1:N)^2)
  return(apply(rg,2,sum))
}

# draw one alpha sample
rh = function(p,a,b,Tol=5000)
{
  tau = sqrt(1/4 + 2*a*p/b^2) + 1/2*ifelse(b>0, 1, -1)
  gamma_rate = tau*abs(b)-b
  for(t in 1:Tol)
  {
    gs = rgamma(1,shape = p, rate=gamma_rate)
    thresh = exp(-a*(gs-tau*abs(b)/(2*a))^2 )
    if(runif(1) < thresh)
      return(gs)
  }
  return(-1)
}

# multinomial_dirichlet
Gibbs_MD = function(data, tau=0, beta=0, niter=1000, N_truncations=300)
{
  M = dim(data)[1]
  K = dim(data)[2]
  p = matrix(1, M, K)
  w = matrix(1, M, K)
  alpha = matrix(1, niter, K)
  
  for(i in 2:niter)
  {
    if(i %% 100 ==0)
      cat("sampling ", i, "\n")
    

    # sample p
    for (m in 1:M)
    {
      p[m,] = pmax(rdirichlet(1, alpha[i-1,] + data[m,]),1e-9) # prevent p=0 case
    }

    # sample eta
    eta = rgamma(M, shape = sum(alpha[i-1,]))
    
    # sample w
    for(k in 1:K)
      w[,k] = rPIG_ERGamma(M,sqrt(2)*alpha[i-1,k],N_truncations)
    
    # sample alpha
    a = apply(w,2,sum)
    b = -M*digamma(1) -beta + sum(log(eta)) + apply(log(p),2,sum)
    for(k in 1:K)
      alpha[i, k] = rh(M+tau, a[k], b[k])
  }
  
  return(alpha)
}

Gibbs_MD_E = function(data, tau=0, beta=0, niter=1000)
{
  M = dim(data)[1]
  K = dim(data)[2]
  alpha = matrix(1, niter, K)
  for(i in 2:niter)
  {
    if(i %% 100 ==0)
      cat("sampling ", i, "\n")
    
    # sample p
    p_e = t(t(data)+alpha[i-1,])
    log_p_e = digamma(p_e) - digamma(apply(p_e,1,sum))
    
    # sample eta
    log_eta_e = digamma(sum(alpha[i-1,]))
    
    # sample w
    w_e = (-digamma(1) + digamma(alpha[i-1,]+1))/(2*alpha[i-1,])
    
    # sample alpha
    a = M*w_e
    b = -M*digamma(1) - beta + M*log_eta_e + apply(log_p_e,2,sum)
    for(k in 1:K)
      alpha[i, k] = rh(M+tau, a[k], b[k])
  }
  
  return(alpha)
}

ltarget = function(data,alpha)
{
  res = sum(apply(data, 1, function(x){
    lgamma(sum(alpha))-lgamma(sum(alpha+x)) + sum(lgamma(alpha+x)-lgamma(alpha))
  }))
  res = res - sum(alpha)
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
    ratio = exp(ltarget(data,propose)-ltarget(data,x[i-1,]))
    if(runif(1)<ratio)
    {
      x[i,] = propose
    }else{
      x[i,] = x[i-1,]
    }
  }
  return(x)
}

