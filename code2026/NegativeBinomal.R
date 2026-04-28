library(GIGrvg)
library(BayesLogit) 
library(mvtnorm)
source("sampling.R")



# MCMC

NBR_mcmc = function(X, y, B=1000, burn=200, N_truncations=300, sigma=1, alpha=2, kappa=1)
{
  N = dim(X)[1]
  k = dim(X)[2]
  
  beta_array = matrix(0, B, k)
  psi_array = rep(0, B)
  
  # parameters
  beta = rnorm(k,s=sigma)
  psi = rgamma(1,shape=2,scale = 1) # gamma prior for negative binomial number of "failures"
  
  
  for(i in 1:(B+burn))
  {
    if(i %% 200 == 0) 
      cat("sampling ", i, "\n")
    
    tau = sapply(y+psi, function(x) rgamma(1,shape=x))
    xi = rPIG_ERGamma(N, psi, N_truncations)
    Xbeta = X %*% beta
    
    omega = diag(N)
    for(j in 1:N) omega[j,j] = rpg(1, (y+psi)[j], Xbeta[j])
    Xomega = t(X) %*% omega
    
    V = solve(diag(k)/sigma^2 + Xomega %*% X) 
    m = V %*% Xomega %*% ((y-psi)/2/diag(omega))
    beta = t(rmvnorm(1,mean=m,sigma = V))
    
    a = sum(xi)
    b = N*(-0.1159315) + sum(log(tau) - Xbeta/2)
    psi = rH(alpha+N, a, b+kappa)
    
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
    
    tau = sapply(y+psi, function(x) rgamma(1,shape=x))
    xi = rPIG_ERGamma(N, psi, N_truncations)
    a = sum(xi)
    b = N*(-digamma(1)) + sum(log(tau) - log(1+exp(Xbeta)))
    psi = rH(alpha+N, a, b+kappa)
    
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
  logL = sapply(r, function(x){-log(x)-N*lgamma(x) + sum(lgamma(y+x))-sum((y+x)*log(1+exp(Xbeta)))-500})
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





