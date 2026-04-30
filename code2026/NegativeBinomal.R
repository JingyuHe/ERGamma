local_lib <- file.path(getwd(), "r-lib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

library(GIGrvg)
library(BayesLogit) 
library(mvtnorm)
source("sampling.R")



log1pexp = function(x)
{
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}


# MCMC

NBR_mcmc = function(X, y, B=1000, burn=200, N_truncations=300, sigma=1, alpha=2, kappa=1)
{
  N = dim(X)[1]
  k = dim(X)[2]
  
  beta_array = matrix(0, B, k)
  psi_array = rep(0, B)
  
  # parameters
  beta = rep(0, k)
  psi = rgamma(1,shape=2,scale = 1) # gamma prior for negative binomial number of "failures"
  
  
  for(i in 1:(B+burn))
  {
    if(i %% 200 == 0) 
      cat("sampling ", i, "\n")
    
    # Partially collapsed update for psi/r.  The Polya-Gamma variables are
    # integrated out in this step; conditioning on them would add a
    # p_PG(xi | y + psi, 0) factor that depends on psi.
    tau = rgamma(N, shape=y+psi, rate=1)
    xi = rPIG_ERGamma(N, psi, N_truncations)
    Xbeta = drop(X %*% beta)
    a = sum(xi)
    b = N*(-digamma(1)) + sum(log(tau)) - sum(log1pexp(Xbeta)) - kappa
    psi_new = rH(alpha+N, a, b)
    if(!is.finite(psi_new) || psi_new <= 0)
      stop("PTN update for psi failed")
    psi = psi_new
    
    # Conditional beta update with fresh Polya-Gamma variables at the new psi.
    omega = sapply(1:N, function(j) rpg(1, y[j]+psi, Xbeta[j]))
    XomegaX = t(X) %*% (X * omega)
    
    V = solve(diag(k)/sigma^2 + XomegaX) 
    m = V %*% t(X) %*% ((y-psi)/2)
    beta = drop(rmvnorm(1,mean=drop(m),sigma = V))
    
    if(i > burn)
    {
      beta_array[i-burn,] = beta
      psi_array[i-burn] = psi
    }
  }
  
  return(list(beta=beta_array, psi=psi_array))
  
}


NBR2_mcmc = function(X, y, beta, B=1000,burn=200, N_truncations=300, sigma=1, alpha=2, kappa=1)
{
  N = dim(X)[1]
  k = dim(X)[2]
  
  psi_array = rep(0, B)
  
  # parameters
  psi = rgamma(1,shape=2,scale = 1) # gamma prior for negative binomial number of "failures"
  Xbeta = X %*% beta
  
  for(i in 1:(B+burn))
  {
    if(i %% 200 == 0) 
      cat("sampling ", i, "\n")
    
    tau = rgamma(N, shape=y+psi, rate=1)
    xi = rPIG_ERGamma(N, psi, N_truncations)
    a = sum(xi)
    b = N*(-digamma(1)) + sum(log(tau)) - sum(log1pexp(Xbeta)) - kappa
    psi_new = rH(alpha+N, a, b)
    if(!is.finite(psi_new) || psi_new <= 0)
      stop("PTN update for psi failed")
    psi = psi_new
    
    if(i > burn)
    {
      psi_array[i-burn] = psi
    }
  }
  
  return(list(psi=psi_array))
  
}





# DGP
# set.seed(5)
k = 5
N = 100
sigma = 1
beta = rnorm(k,s=sigma)
X = matrix(rnorm(N*k,s=0.5),N,k)
Xbeta = X %*% beta
p = 1/(1+exp(Xbeta))
psi = 5
y = sapply(p, function(x) rnbinom(1,size=psi,prob=x))

dalpha=function(r)
{
  logL = sapply(r, function(x){-log(x)-N*lgamma(x) + sum(lgamma(y+x))-sum((y+x)*log1pexp(Xbeta))-500})
  return(exp(logL))
}



res = NBR_mcmc(X, y,B=2000,burn=2000,alpha=0, kappa=0, N_truncations=1000, sigma = 2000)

boxplot(res$beta, main="Posterior Beta")
points(1:k, beta, col="red",cex=1.5, pch=19)

hist(res$psi, breaks = 25, freq = FALSE, xlab = "",
     main="Histogram of Posterior r", col="cornflowerblue")

res2 = NBR2_mcmc(X, y, beta, B=2000,  burn=2000, alpha=0, kappa=0, N_truncations=2000)

hist(res2$psi, breaks = 25, freq = FALSE, xlab = "",
     main="Histogram of Posterior r (beta is known)", col="cornflowerblue")

x = seq(3, 6.5, 0.01)
Z = integrate(dalpha,lower = 4,upper = 6.5)
lines(x,sapply(x, dalpha)/Z$value,col="red",lwd=3)





