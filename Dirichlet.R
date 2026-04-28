# Simulation of homogeneous alpha



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



# data generating process
alpha_true = 1
K = 10 # number of categories
N = 5000 # number of data observations
p_true = rdirichlet(1, rep(alpha_true, K))
n = rmultinom(1, N, prob = p_true)


gamma_const = 0.577216
N_truncations = 500 # number of sum of GIG for PIG
N_draw = 1000

w = matrix(0, N_draw, K)
alpha = rep(1, N_draw)
eta = rep(0, N_draw)
p = matrix(1 / K, N_draw, K)
b_vec = rep(0, N_draw)
a_vec = rep(0, N_draw)
# Gibbs sampler
alpha2 = rep(1, N_draw)

for (i in 2:N_draw) {
  # sampling eta
  eta[i] = rgamma(1, shape = K * alpha[i - 1] + sum(n), scale = 1)

  # sampling w
  for (j in 1:K) {
    # w[i, j] = rPIG_ERGamma(sqrt(2 * (alpha[i - 1] + n[j] - 1)^2), N_truncations)
    w[i,j] = 0.2
  }

  # sampling p
#   p[i,] = rdirichlet(1, rep(alpha[i - 1] + sum(n), K))
    p[i, ] = rdirichlet(1, alpha[i - 1] + n)
    # p[i, ] = p_true

  # sampling alpha
  a = sum(w[i,])
  b = (-2 * sum((n - 1) * w[i,]) + K * log(eta[i]) + sum(gamma_const + log((p[i,]))))
  b = b
  b_vec[i] = b
  a_vec[i] = a
  # alpha[i] = rnorm(1, mean = b / 2 / a, sd = sqrt(1 / 2 / a))
  alpha[i] = rtruncnorm(1, a=0, mean = b / 2 / a, sd = sqrt(1 / 2 / a))
    # alpha[i] = alpha_true

}


par(mfrow = c(4,4))

for(i in 1:K){
    plot(p[, i])
    abline(h = p_true[i], col = "red", main = "p")
}


plot(alpha)
abline(h = alpha_true, col = "red", main = "alpha")


plot(alpha2)
abline(h = alpha_true, col = "red", main = "alpha2")


plot(b_vec, main = "b")
plot(a_vec, main = "a")
