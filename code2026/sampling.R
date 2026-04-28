library(GIGrvg)
rPIG_ERGamma = function(n, cc, N) {
  # set dk = k
  # n is number of rvs
  # c is a scalar
  # N is number of truncation
  
  output = rep(0,n)
  
  if(N > 0)
  {
    for (i in 1:N) {
      if(cc == 0)
      {
        output = output + 1/rgamma(n, shape=3/2)/(4*i^2)
      }else
      {
        output = output + rgig(n = n, lambda = -3 / 2, chi = 1 / 2 / i ^ 2, 2*cc ^ 2)
      }
      
    }
  }
  if(cc == 0)
  {
    temp1 = trigamma(1+N)
    temp2 = psigamma(1+N,deriv = 2)
    b = -4*temp1/temp2
    a = 0.5*temp1*b
  }else{
    temp1 = digamma(1+N+cc) - digamma(1+N)
    temp2 = temp1 - cc*trigamma(1+N+cc)
    b = 2*cc^2*temp1 / temp2
    a = b/(2*cc)*temp1
  }

  output = output + rgamma(n,a,rate=b)
  return(output)
}



rH = function(m,a,b,max_try = 500)
{
  tau = ifelse(b>0, 0.5+sqrt(1/4+2*a*m/b^2), -0.5+sqrt(1/4+2*a*m/b^2))
  # cat("m=",m," a=", a, " b=", b, " tau=",tau, "\n")
  v1 = tau*abs(b)-b
  v2 = tau*abs(b)/(2*a)
  
  for(i in 1:max_try)
  {
    x = rgamma(1,m,rate=v1)
    # cat("x=",x," v1=", v1, " v2=", v2, "\n")
    if(runif(1) <= exp(-a*(x-v2)^2))
      return(x)
  }
}


a = rPIG_ERGamma(500,cc=1,N=200)
b = rPIG_ERGamma_2(500,cc=1,N=200)
hist(a,breaks = 30)
hist(b,breaks = 30)


calc_MGF = function(s,sample)
{
  # calculate sample expectation of exp(-s^2*X)
  sapply(s, FUN = function(x) mean(exp(-x^2*sample)))
}


par(mfrow=c(2,3))
for(c in c(0.01,0.05,0.1,0.5,1,5))
{
  s = seq(0,5,0.01)
  mgf = exp(lgamma(1+c)-lgamma(1+sqrt(s^2+c^2))+digamma(1)*(sqrt(s^2+c^2)-c))
  plot(s,mgf,"n",lwd=2,ylab="E[e^(-s^2*X)]",main=paste("c =",c))
  N_list = c(0,1,5,10,50,100,300)
  for(i in 1:7)
  {
    sample = rPIG_ERGamma(1000,cc=c,N=N_list[i])
    res = calc_MGF(s,sample)
    lines(s, res, col=i+1,lwd=2)
  }
  lines(s, mgf, lwd=2)
  legend("topright",legend=paste("N =",N_list),lwd=2,col=2:8)
}







# rPIG_ERGamma_2 = function(n, cc, N, max_try = 500){ # Not better than the original
#   # first sample ERG(0)
#   output = rep(0, n)
#   for(i in 1:n)
#   {
#     for(b in 1:max_try)
#     {
#       U = runif(1)
#       W = sum(1/rgamma(N,shape = 3/2)/(4*(1:N)^2))
#       
#       if(U < exp(-cc^2*W))
#       {
#         break
#       }
#     }
#     if(b == max_try)
#     {
#       print("no acceptance")
#     }
#     output[i] = W
#   }
#   return(output)
# }



# rPIG_ERGamma_saddle = function(n, cc, tol = 1e-6, max.T=500) # Not working correctly
# {
#   
#   gamma_const = -digamma(1)
#   const = gamma_const*cc + lgamma(1+cc)
#   output = rep(0,n)
#   for(i in 1:n)
#   {
#     U = runif(1)
#     t = cc/2
#     for(b in 1:max.T)
#     {
#       print(t)
#       ct = sqrt(cc^2-t)
#       k = -gamma_const*ct - lgamma(1+ct) + const
#       numerator = gamma_const+digamma(1+ct)
#       k1 = numerator/2/ct
#       k2 = (numerator-ct*trigamma(1+ct))/4/ct^3
#       
#       w = sign(t)*sqrt(2*(t*k1 - k))
#       u = t*sqrt(k2)
#       Ft = pnorm(w) + dnorm(w)*((1/w) - 1/(u))
#       
#       delta = (Ft-U)/(sqrt(k2)*dnorm(w))
#       
#       if(abs(delta) < tol)
#       {
#         t = t-delta
#         break 
#       }
#       if(t - delta >= cc)
#       {
#         delta = (t-cc)/2
#       }
#       
#       t = t-delta
#     }
#     ct = sqrt(cc^2-t)
#     output[i] = (gamma_const+digamma(1+ct))/2/ct
#   }
#   
#   return(output)
# }
# 
