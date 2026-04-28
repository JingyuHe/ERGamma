# todo list:
# 1. calculate constant b
# 2. simulate tau from M
# 3. calculate the tail expectation

Bodesson = function(N, dk, T=50)
{
  b = 
   # simulate poisson points
  T_i = matrix(runif(N*T),N,T)*rpois(N, b*T)
  W_i = rgamma(N*T, N, T)
  tau_i = matrix( , N, T)
  X_matrix = exp(T_i)*W_i/tau_i
  final_sample = apply(X_matrix, 1, sum)
  return(final_sample)
}

