source("sampling.R")
library(ggplot2)
library(reshape2)
library(moments)

mm = function(x)
{
  cat("mean = ",mean(x),"\n")
  cat("variance = ",var(x),"\n")
  cat("skewness = ",skewness(x),"\n")
  cat("kurtosis = ",kurtosis(x), "\n")
}

# Miller's derivative matching strategy
calc_AB = function(delta, mu, init.x=1, tol=1e-3,max.T=500)
{
  x = init.x
  for(i in 1:max.T)
  {
    A = 1 - x^2*delta*(delta*trigamma(delta*x+1)-trigamma(x))
    B = (A-1)/x - delta*(digamma(delta*x+1) - digamma(x) -log(mu*delta))
    x_new = A/B
    if(abs(x/x_new - 1) < tol)
    {
      return(c(A,B))
    }
    x = x_new
  }
  print("not convergent")
  return(c(A,B))
}

calc_gamma_mom = function(A,B)
{
  m1 = A/B
  m2 = A/B^2
  m3 = 2/sqrt(A)
  m4 = 6/A
  return(c(m1,m2,m3,m4))
}

# calculate normalization constants
Zs = c()
for(delta in c(5,10,30))
{
  if(delta == 30)
  {
    mu = 5.09/4.26
  }
  if(delta == 10)
  {
    mu = 5.57/5.01
  }
  if(delta == 5)
  {
    mu = 7.19/6.05
  }
  b_pre = delta*gamma_const - delta*log(delta*mu)
  dalpha = function(x)
  {
    exp(lgamma(delta*x+1)-delta*lgamma(x)) * (delta*mu)^(-delta*x)
  }
  
  if(delta<30)
  {
    Z = integrate(dalpha,lower = 0,upper = 20)
  }else{
    Z = integrate(dalpha,lower = 1,upper = 6)
  }
  
  Zs = c(Zs,Z$value)
}


dalpha = function(x,delta,mu,z)
{
  exp(lgamma(delta*x+1)-delta*lgamma(x)) * (delta*mu)^(-delta*x) / z
}

N=5000
B=5000
gamma_const = -digamma(1)
xi_data = data.frame(x=rep(0,N))
xi_init = 5 # initial start

for(delta in c(5,10,30))
{
  if(delta == 30)
  {
    mu = 5.09/4.26
  }
  if(delta == 10)
  {
    mu = 5.57/5.01
  }
  if(delta == 5)
  {
    mu = 7.19/6.05
  }
  b_pre = delta*gamma_const - delta*log(delta*mu)
  xi = xi_init
  xi_sample = rep(0,N)
  for(i in 1:(N+B))
  {
    if(i %% 200 ==0) print(i)
    tau = rgamma(1,shape = delta*xi+1)
    omega = rPIG_ERGamma(delta, xi, 1000)
    b = b_pre + delta*log(tau)
    xi = rH(delta+1, sum(omega), b)
    if(i > B)
      xi_sample[i-B] = xi
  }
  xi_data[,paste0("n=",delta)] = xi_sample
}

# Miller
for(delta in c(5,10,30))
{
  if(delta == 30)
  {
    mu = 5.09/4.26
  }
  if(delta == 10)
  {
    mu = 5.57/5.01
  }
  if(delta == 5)
  {
    mu = 7.19/6.05
  }
  
  ab = calc_AB(delta,mu)
  print(mm(rgamma(5000,shape=ab[1],rate = ab[2])))

}



xi_data = xi_data[,-1]
save(file="gamma_shape_data.Rdata",xi_data)
p_data = melt(xi_data)
colnames(p_data) = c("sample_size","alpha")


p = ggplot(p_data, aes(x = alpha, fill = sample_size)) +
  geom_histogram(aes(y=..density..), position = "identity", alpha = 0.6, bins = 70)+
  geom_function(fun = dalpha, size=1, colour = "red", args=list(5,7.19/6.05,Zs[1])) + 
  geom_function(fun = dalpha, size=1, colour = "green", args=list(10,5.57/5.01,Zs[2])) + 
  geom_function(fun = dalpha, size=1, colour = "blue", args=list(30,5.09/4.26,Zs[3]))+
  geom_vline(xintercept = 5, linetype="dotted", size=1.5)+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()
pdf("gamma_shape_unknown_scale.pdf",width = 12, height = 6)
p
dev.off()

pdf("gamma_trace.pdf",width = 12, height = 5)
par(mfrow=c(1,3))
plot(tail(xi_data[,1],2000), main="n = 5",
     xlab = "iterate", ylab = "alpha", pch=19, cex=0.35,ylim = c(0,18))
abline(h=5,col="red",lwd=3,lty=2)
plot(tail(xi_data[,2],2000),main="n = 10",
     xlab = "iterate", ylab = "alpha", pch=19, cex=0.35,ylim = c(0,18))
abline(h=5,col="red",lwd=3,lty=2)
plot(tail(xi_data[,3],2000),main="n = 30",
     xlab = "iterate", ylab = "alpha", pch=19, cex=0.35,ylim = c(0,18))
abline(h=5,col="red",lwd=3,lty=2)
dev.off()

mm(xi_data[,1])

# EM
set.seed(123)
n = 500; alpha = 3; beta = 1.5
B = 100
xpath = rep(NA, B)
plot(1:B, xpath, type="n", xlab="step", ylab="alpha", main="Path of alpha",
     ylim=c(0,8))
y = rgamma(n,shape=alpha,rate=beta)
ya = mean(y)
yg = (prod(y))^(1/n)
# prior
delta = 0
eta = 2
mu = 1
# posterior
delta_new = delta + n
eta_new = delta/(delta+n)*eta + n/(delta+n)*ya
mu_new = mu^(delta/(delta+n))*yg*(n/(delta+n))

for(j in 1:30)
{
  # algo
  B = 100
  x = rnorm(1,mean=3,sd=2)
  xpath = rep(0,B)
  kappa_1 = delta_new
  k2_pre = delta_new*gamma_const
  k3_pre = k2_pre - delta_new*log(delta_new*eta_new/mu_new)
  
  for(i in 1:B)
  {
    kappa_2 = (k2_pre + delta_new*digamma(x+1))/(2*x)
    kappa_3 = k3_pre + delta_new*digamma(delta_new*x+1)
    x = (kappa_3 + sqrt(kappa_3^2+8*kappa_1*kappa_2))/(4*kappa_2)
    xpath[i] = x
  }
  lines(1:B, xpath, col=j)
}


# EM2
B = 100
xpath = rep(NA, B)
plot(1:B, xpath, type="n", xlab="step", ylab="alpha", main="Path of alpha",
     ylim=c(0,8))
mu = 5.09/4.26
delta = 30
post_mode = c()

for(j in 1:30)
{
  # algo
  B = 100
  x = rnorm(1,mean=3,sd=2)
  xpath = rep(0,B)
  kappa_1 = delta
  k2_pre = delta*gamma_const
  k3_pre = k2_pre - delta*log(delta*mu)
  
  for(i in 1:B)
  {
    kappa_2 = (k2_pre + delta*digamma(x+1))/(2*x)
    kappa_3 = k3_pre + delta*digamma(delta*x+1)
    x = (kappa_3 + sqrt(kappa_3^2+8*kappa_1*kappa_2))/(4*kappa_2)
    xpath[i] = x
  }
  post_mode = c(post_mode,xpath[B])
  lines(1:B, xpath, col=j)
}



# Miller
delta = 20
mu = 20

ab = calc_AB(delta,mu)
dalpha = function(x)
{
  exp(lgamma(delta*x+1)-delta*lgamma(x) - delta*x*log(delta*mu))
}
upper = 1
x_grid = seq(0.001,upper,length.out=1000)
Z = integrate(dalpha,0.001,upper)$value
plot(x_grid, dalpha(x_grid)/Z,"l")
lines(x_grid, dgamma(x_grid,shape=ab[1],rate=ab[2]),lty=2,col="red",lwd=2)

# counter example
y = function(x)
{
  gamma(x)/gamma(0.04*x)^2*exp(-0.06*x^2-1.28*x)
}
x = 10
tol = 1e-3
for(i in 1:10000)
{
  A = 1 - x^2*(-0.12-0.0032*trigamma(0.04*x)+trigamma(x))
  B = (A-1)/x - (-0.12*x-0.08*digamma(0.04*x)+digamma(x)-1.28)
  x_new = A/B
  if(abs(x/x_new - 1) < tol)
  {
    break
  }
  x = x_new
}


upper = 25
x_grid = seq(0.001,upper,length.out=1000)
Z = integrate(y,0.001,upper)$value
plot(x_grid, y(x_grid)/Z,"l",ylim=c(0,0.12))
lines(x_grid, dgamma(x_grid,shape=A,rate=B),lty=2,col="red",lwd=2)

N = 20000
sample = rep(0,N)
x = 10
for(i in 1:N)
{
  eta = rgamma(1,shape=x)
  w = rPIG_ERGamma(2,0.04*x,N=500)
  x = rH(3,0.04^2*sum(w)+0.06,0.08*gamma_const-1.28+log(eta))
  sample[i] = x
  if(i %% 100 == 0) print(i)
}

pd = data.frame(value=tail(sample,5000))
p <- ggplot(pd, aes(x=value)) + 
  geom_histogram(aes(y=..density..),position = "identity", alpha = 0.6, binwidth = 0.5)
p

pdf("counter.pdf",height=5,width = 7)
hist(tail(sample,20000),breaks=30,freq = FALSE,main="",xlab="",
     ylim=c(0,0.1),col="lightblue")
lines(x_grid, y(x_grid)/Z,lwd=2)
lines(x_grid, dgamma(x_grid,shape=A,rate=B),lty=2,col=2,lwd=2)
legend("topright",legend=c("true density","Miller method"),lwd=2,col=1:2,lty=1:2)
dev.off()
