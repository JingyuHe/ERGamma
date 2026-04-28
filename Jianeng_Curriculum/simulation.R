source("~/dropbox/Curriculum/sampling.R")

# dirichlet heterogenous multinomial
set.seed(123)
K = 50 # number of categories
alpha_true = 1:K
M = 2000
data = matrix(0, M,K)
p_data = matrix(0, M,K)
for(i in 1:M)
{
  p_true = rdirichlet(1, alpha_true)
  data[i,] = rmultinom(1, rpois(1,500), prob = p_true)
  p_data[i,] = p_true
}

ptm <- proc.time()
gibbs_alpha_1 = Gibbs_MD(data, niter=2000)
proc.time() - ptm

ptm <- proc.time()
gibbs_alpha_2 = Gibbs_MD_E(data, niter=2000)
proc.time() - ptm

ptm <- proc.time()
MH_alpha = MH(data,ltarget, tau = 1, step=0.1, niter=2000)
proc.time() - ptm



mean((alpha_true-colMeans(gibbs_alpha_2[1500:2000,]))^2)
mean((alpha_true-colMeans(MH_alpha[1500:2000,]))^2)


acf(MH_alpha[,3])
acf(gibbs_alpha_2[,3])
