library(GIGrvg)
library(truncnorm)
library(ggplot2)

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


gamma_const = 0.577216
alpha = 3
beta = 2
# set.seed(123)
set.seed(789)
y = rgamma(200,shape = alpha,rate = beta)

# prior
a = 2
b = 1
c = 1


N = 2000
a_prime = a*prod(y)
b_prime = b+length(y)
c_prime = c+length(y)
beta_y = a_prime*beta^c_prime

mu_numerator = gamma_const*b_prime+log(beta_y)
if(mu_numerator == Inf)
{
  mu_numerator = gamma_const*b_prime+log(a)+sum(log(y))+c_prime*log(beta)
}

alpha_path = 2
for(i in 1:N)
{
  print(i)
  t_alpha = alpha_path[i] - 1
  w = rPIG_ERGamma(b_prime,sqrt(2)*t_alpha,N=200)
  denominator = 2*sum(w)
  t_alpha = rtruncnorm(1, a=-1,
                       mean = mu_numerator/denominator,
                       sd = sqrt(1 / denominator))
  alpha_path = c(alpha_path,t_alpha + 1)
}


##############

dalpha = function(x)
{
  prod(beta^x/gamma(x)*y^(x-1)*exp(-beta*y))
}

dalpha2 = function(x)
{
  1/gamma(x)^b_prime * beta_y^x
}

x = seq(2.5,3.5,0.01)
Z = integrate(dalpha2,lower = 2,upper = 4)


par(mfrow = c(1,2))

# pdf("gamma_shape.pdf",height=5,width = 5)
hist(alpha_path[1501:2001],breaks = 30,freq = FALSE,
     main = "Posterior Gamma Shape (beta=2, a=2, b=1, c=1)",
     xlab = "alpha")
lines(density(alpha_path[1501:2001]),lwd=2, lty = 2)
lines(x,sapply(x, dalpha2)/Z$value,col="red",lwd=2)

alpha_plot = data.frame(alpha=alpha_path[1501:2001])
h <- ggplot(alpha_plot, aes(x=alpha)) + 
  scale_fill_brewer()+
  geom_histogram(aes(y =..density..),
                 color="darkblue",fill="#7da7ca",bins=30) +
  geom_density(size=1,linetype="dashed") + 
  stat_function(fun = function(x) dalpha2(x)/Z$value,
                color="red",size=1) + 
  ggtitle("Posterior Gamma Shape (beta=2, a=2, b=1, c=1)") + 
  theme(plot.title = element_text(hjust = 0.5))
print(h)
# dev.off()










#########################################################
# Random walk metropolis Hastings

# posterior parameters
n = length(y)
b2 = b + n
c2 = c + n * log(beta) + sum(log(y))

log_post = function(alpha){
  output = c2 * alpha - b2 * lgamma(alpha)
  return(output)
}

alpha_vec = rep(1, N)

step_size = 0.05

# Start from correct value!!! easier
alpha_vec[1] = 2

for(i in 2:N){
  alpha_prop = exp( log(alpha_vec[i-1]) + rnorm(1, sd = step_size) )
  MH_ratio = log_post(alpha_prop) - log_post(alpha_vec[i-1])
  MH_ratio = exp(MH_ratio)
  if(runif(1) < MH_ratio){
    alpha_vec[i] = alpha_prop
  }else{
    alpha_vec[i] = alpha_vec[i-1]
  }
}


# hist(alpha_vec[1501:2000],breaks = 30,freq = FALSE,
    #  main = "Metropolis-Hastings (beta=2, a=2, b=1, c=1)",
    #  xlab = "alpha")
# lines(density(alpha_vec[1501:2000]),lwd=2, lty = 2)
# lines(x,sapply(x, dalpha2)/Z$value,col="red",lwd=2)

pdf("gamma_shape_MH.pdf",height=5,width = 5)

alpha_plot = data.frame(alpha=alpha_vec[1501:2000])
h <- ggplot(alpha_plot, aes(x=alpha)) + 
  scale_fill_brewer()+
  geom_histogram(aes(y =..density..),
                 color="darkblue",fill="#7da7ca",bins=30) +
  geom_density(size=1,linetype="dashed") + 
  stat_function(fun = function(x) dalpha2(x)/Z$value,
                color="red",size=1) + 
  ggtitle("MH (beta=2, a=2, b=1, c=1)") + 
  theme(plot.title = element_text(hjust = 0.5))
print(h)
dev.off()


pdf("gamma_shape_trace.pdf",height=5,width =8)
par(mfrow = c(1,2))
plot(alpha_path, main = "Trace plot of PIG samples", bty = "l")
plot(alpha_vec, main = "Trace plot of MH samples", bty = "l")
dev.off()

pdf("gamma_shape_acf.pdf",height=5,width =8)
par(mfrow = c(1,2))
acf(alpha_path, bty = "l", main = "Autocorrelation of P-IG draws")
acf(alpha_vec, bty = "l", main = "Autocorrelation of M-H draws")
dev.off()


library(coda)
# effective sample size
effectiveSize(alpha_path) / length(alpha_path)
effectiveSize(alpha_vec) / length(alpha_vec)

