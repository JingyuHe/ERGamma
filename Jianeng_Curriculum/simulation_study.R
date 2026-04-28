source("./sampling.R")
library(gtools)
library(coda)
simulation = function(M,K,set="A",niter=500)
{
  # DGP
  data = matrix(0, M,K)
  p_data = matrix(0, M,K)
  if(set=="A")
  {
    alpha_true = 1:K/K
  }
  if(set=="B")
  {
    alpha_true = rep(1/K,K)
  }
  
  for(i in 1:M)
  {
    p_true = rdirichlet(1, alpha_true)
    data[i,] = rmultinom(1, rpois(1,500), prob = p_true)
    p_data[i,] = p_true
  }
  
  res = matrix(0,3,3)
  
  ptm <- proc.time()
  gibbs_alpha_1 = Gibbs_MD(data, niter=niter)
  a = proc.time() - ptm
  res[1,2] = a[3]
  res[1,3] = mean(effectiveSize(gibbs_alpha_1)/niter)
  res[1,1] = mean((alpha_true-colMeans(gibbs_alpha_1[101:niter,]))^2)
  save(gibbs_alpha_1, file = paste0(M,K,set,"_PIG.Rdata"))
  
  ptm <- proc.time()
  gibbs_alpha_2 = Gibbs_MD_E(data, niter=niter)
  b = proc.time() - ptm
  res[2,2] = b[3]
  res[2,3] = mean(effectiveSize(gibbs_alpha_2)/niter)
  res[2,1] = mean((alpha_true-colMeans(gibbs_alpha_2[101:niter,]))^2)
  save(gibbs_alpha_2, file = paste0(M,K,set,"_PIGE.Rdata"))
  
  ptm <- proc.time()
  MH_alpha = MH(data,ltarget, tau = 1, step=0.1, niter=niter)
  c = proc.time() - ptm
  res[3,2] = c[3]
  res[3,3] = mean(effectiveSize(MH_alpha)/niter)
  res[3,1] = mean((alpha_true-colMeans(MH_alpha[101:niter,]))^2)
  save(MH_alpha, file = paste0(M,K,set,"_MH.Rdata"))
  
  return(res)
}

M_list = c(100,1000)
K_list = c(10,50)
res_matrix = matrix(0,length(M_list)*length(K_list),8)
time_matrix <- size_matrix<- res_matrix

count = 1
for(M in M_list)
  for(K in K_list)
  {
    res_matrix[count,1] = M
    time_matrix[count,1] = M
    res_matrix[count,2] = K
    time_matrix[count,2] = K
    
    cat(M,"\t", K, "\t A \n")
    resA = simulation(M,K,"A")
    res_matrix[count,1:3+2] = resA[,1]
    time_matrix[count,1:3+2] = resA[,2]
    size_matrix[count,1:3+2] = resA[,3]
    
    cat(M,"\t", K, "\t B \n")
    resB = simulation(M,K,"B")
    res_matrix[count,4:6+2] = resB[,1]
    time_matrix[count,4:6+2] = resB[,2]
    size_matrix[count,4:6+2] = resB[,3]  
  }

save(res_matrix, time_matrix, size_matrix, file = "simulation_res.Rdata")