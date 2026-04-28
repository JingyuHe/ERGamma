source("sampling.R")
library(ggplot2)


alpha = 3
beta = 2
# set.seed(123)
set.seed(789)
y = rgamma(200,shape = alpha,rate = beta)

# prior
b = -1
c = 1

N = 2000

b_prime = b+length(y)
c_prime =c+length(y)*log(beta)+sum(log(y))


alpha_path = rep(2,N)
for(i in 2:N)
{
  print(i)
  w = rPIG_ERGamma(b_prime,alpha_path[i-1]*sqrt(2),N=500)
  
  alpha_path[i] = rh(b_prime+1, sum(w), c_prime-b_prime*digamma(1))
}


##############

dalpha2 = function(x)
{
  exp(c_prime*x - b_prime*lgamma(x))
}

x = seq(2.5,3.5,0.01)
Z = integrate(dalpha2,lower = 2.5,upper = 3.5)


par(mfrow = c(1,2))

# pdf("gamma_shape.pdf",height=5,width = 5)
hist(alpha_path[1501:2000],breaks = 30,freq = FALSE,
     main = "Posterior Gamma Shape (beta=2, b=1, c=1)",
     xlab = "alpha")
lines(density(alpha_path[1501:2000]),lwd=2, lty = 2)
lines(x,sapply(x, dalpha2)/Z$value,col="red",lwd=2)

alpha_plot = data.frame(alpha=alpha_path[1501:2000])
h <- ggplot(alpha_plot, aes(x=alpha)) + 
  scale_fill_brewer()+
  geom_histogram(aes(y =..density..),
                 color="darkblue",fill="#7da7ca",bins=30) +
  geom_density(size=1,linetype="dashed") + 
  stat_function(fun = function(x) dalpha2(x)/Z$value,
                color="red",size=1) + 
  ggtitle("Posterior Gamma Shape (beta=2, b=-1, c=1)") + 
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
h2 <- ggplot(alpha_plot, aes(x=alpha)) + 
  scale_fill_brewer()+
  geom_histogram(aes(y =..density..),
                 color="darkblue",fill="#7da7ca",bins=30) +
  geom_density(size=1,linetype="dashed") + 
  stat_function(fun = function(x) dalpha2(x)/Z$value,
                color="red",size=1) + 
  ggtitle("MH (beta=2, b=-1, c=1)") + 
  theme(plot.title = element_text(hjust = 0.5))
print(h2)
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

