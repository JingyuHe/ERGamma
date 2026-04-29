local_lib <- file.path(getwd(), "r-lib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

library(doParallel)
source("dirichlet_new.R")


simulation = function(S,K,b0,N_sample,alpha_type="homo")
{
  if(alpha_type=="homo")
    alpha = rep(1/K,K)
  else
    alpha = (1:K)/K
  
  # simulate data
  data = matrix(0, S, K)
  pdata = matrix(0, S,K)
  for(s in 1:S)
  {
    pdata[s,] = rdirichlet(1, alpha)
    data[s,] = rmultinom(1, 500, prob = pdata[s,])
  }
  
  Time = rep(0,0,0,0)
  
  # gibbs
  Time[1] = system.time({gibbs = Gibbs(data, a0=b0/K, b0=b0, N=N_sample)})[3]
  
  # gibbs clt
  Time[2] = system.time({gibbs_clt = Gibbs_CLT(data, a0=b0/K, b0=b0, N=N_sample)})[3]

  
  # gibbs e
  Time[3] = system.time({gibbs_e = Gibbs_E(data, a0=b0/K, b0=b0, N=N_sample)})[3]
  
  # mh
  Time[4] = system.time({mh_alpha = MH(data, ltarget, a0=b0/K, b0=b0, step=0.2, N=N_sample)})[3]
  
  return(list("time"=Time, "gibbs"=gibbs, "clt"=gibbs_clt, "e"=gibbs_e, "mh"=mh_alpha))
}


# simulation(100, 50, 1, N_sample=500, alpha_type = "homo")
