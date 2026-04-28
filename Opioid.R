rm(list = ls())
library(GIGrvg)
library(truncnorm)
library(gtools)
library(ggplot2)
library(reshape2)

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


# opioid data ---------------
OverdoseDeaths = read.csv("OverdoseDeaths.csv", header=T, stringsAsFactors = F, colClasses = c("character", "character", "integer", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
data = as.matrix(OverdoseDeaths[,5:10])

OverdoseDeaths = OverdoseDeaths[,2:10]
rates = OverdoseDeaths[,4:9]/apply(OverdoseDeaths[,4:9],1,sum)

df = data.frame(value = c(rates[,1], rates[,2], rates[,3], rates[,4], rates[,5], rates[,6]), Drug = rep(c("Cocaine", "Synthetic", "Heroin", "Methadone", "NatOpioids", "PsyStim"), each = 76) )
p0 = ggplot(df, aes(x = value, fill = Drug)) + 
  geom_histogram(color="darkblue")+
  xlab("Proportion of overdose deaths")

pdf("figure2.pdf",width = 6,height = 6)
print(p0)
dev.off()

data_colsum = apply(data,2,sum)
data_rowsum = apply(data,1,sum)

## **********
data[62,6] = 1
## **********

# gibbs sampler ---------
gamma = 0.577216
N_truncations = 500 # number of sum of GIG for PIG
N_draw = 1000

# number of probability vectors
M = dim(data)[1]
# number of alphas (category)
K = dim(data)[2]

opioid_alpha = function(prior_sigma2 = 1, lambda=0)
{
  w = array(0, dim = c(N_draw, M, K))
  alpha = matrix(1, N_draw, K)
  eta = matrix(1, N_draw, M)
  p = array(1 / K, dim = c(N_draw, M, K))
  b_vec = matrix(0, N_draw, K)
  a_vec = matrix(0, N_draw, K)
  
  # Gibbs sampler
  for (i in 2:N_draw) {
    cat("sampling ", i, "\n")
    
    # sampling eta
    for (m in 1:M) {
      eta[i, m] = rgamma(1, shape = sum(alpha[i - 1,]) + data_rowsum[m])
    }

    
    # sampling w
    for (m in 1:M) {
      for (k in 1:K) {
        w[i, m, k] = rPIG_ERGamma(sqrt(2 * (alpha[i - 1, k] + data[m, k] - 1)^2), N_truncations)
      }
    }

    # sampling p
    for (m in 1:M) {
      p[i, m,] = rdirichlet(1, alpha[i - 1,] + data[m,])
    }
    
    
    # sampling alpha
    a = apply(w[i,,], 2, sum) + ifelse(prior_sigma2==0, 0, 1/(2*prior_sigma2))
    b = -2 * apply((data - 1) * as.matrix(w[i,,]), 2, sum) + sum(log(eta[i,])) + M * gamma + apply(log(p[i,,]), 2, sum) - lambda
    b_vec[i,] = b
    a_vec[i,] = a
    
    for (k in 1:K) {
      alpha[i,k] = rtruncnorm(1, a = 0, mean = b[k] / 2 / a[k], sd = sqrt(1 / 2 / a[k]))
    }

  }
  
  return(alpha[701:N_draw,])
}

# sampling alpha ---------
alpha_1 = opioid_alpha(prior_sigma2=0.20888569^2)
alpha_2 = opioid_alpha(prior_sigma2=0.10444284^2)
alpha_3 = opioid_alpha(prior_sigma2=0.06962856^2)

names = colnames(data)
alpha_sample = rbind(data.frame(alpha = c(alpha_1), tau2 = 0.044, drug=rep(names,each=300)),
                     data.frame(alpha = c(alpha_2), tau2 = 0.011, drug=rep(names,each=300)),
                     data.frame(alpha = c(alpha_3), tau2 = 0.005, drug=rep(names,each=300)))
alpha_sample$tau2 = paste0("tau2=",alpha_sample$tau2)

p = ggplot(alpha_sample,aes(x=alpha,fill=drug)) + 
  geom_histogram(color="darkblue",bins=50) +
  facet_grid( drug ~ tau2, scales = "free") + 
  ggtitle("Posterior Alpha Histogram (Normal Prior)") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")

pdf("opioid_normal_new.pdf", height = 9, width = 9)
print(p)
dev.off()


alpha2_5 = opioid_alpha(prior_sigma2=0,lambda=5) # doesn't work when M = 76
alpha2_10 = opioid_alpha(prior_sigma2=0,lambda=10) # doesn't work when M = 76
alpha2_20 = opioid_alpha(prior_sigma2=0,lambda=20) # doesn't work when M = 76
alpha2_48 = opioid_alpha(prior_sigma2=0,lambda=48)
alpha2_45 = opioid_alpha(prior_sigma2=0,lambda=45)
alpha2_75 = opioid_alpha(prior_sigma2=0,lambda=75)
alpha2_100 = opioid_alpha(prior_sigma2=0,lambda=100)

alpha_sample2 = rbind(data.frame(alpha = c(alpha2_50), lambda = 50, drug=rep(names,each=300)),
                     data.frame(alpha = c(alpha2_75), lambda = 75, drug=rep(names,each=300)),
                     data.frame(alpha = c(alpha2_100), lambda = 100, drug=rep(names,each=300)))
alpha_sample2$lambda = as.factor(paste0("lambda=",alpha_sample2$lambda))
alpha_sample2$lambda = factor(alpha_sample2$lambda, levels=c("lambda=50","lambda=75","lambda=100"))

p2 = ggplot(alpha_sample2,aes(x=alpha,fill=drug)) + 
  geom_histogram(color="darkblue",,bins=50) +
  facet_grid( drug ~ lambda, scales = "free") + 
  ggtitle("Posterior Alpha Histogram (Exponential Prior)") + 
  theme(plot.title = element_text(hjust = 0.5),,legend.position = "none")

pdf("opioid_exponential.pdf", height = 9, width = 9)
print(p2)
dev.off()

# posterior predictive ------------------------
post_p_1 = NULL
post_p_2 = NULL
post_p_3 = NULL

post_p2_50 = NULL
post_p2_75 = NULL
post_p2_100 = NULL

N_post_alpha = dim(alpha_1)[1]
for(i in 1:1000)
{
  post_p_1 = rbind(post_p_1, rdirichlet(1, alpha_1[sample(N_post_alpha,1),]))
  post_p_2 = rbind(post_p_2, rdirichlet(1, alpha_2[sample(N_post_alpha,1),]))
  post_p_3 = rbind(post_p_3, rdirichlet(1, alpha_3[sample(N_post_alpha,1),]))
  # post_p2_50 = rbind(post_p2_50, rdirichlet(2, alpha2_50[sample(N_post_alpha,1),]))
  # post_p2_75 = rbind(post_p2_75, rdirichlet(2, alpha2_75[sample(N_post_alpha,1),]))
  # post_p2_100 = rbind(post_p2_100, rdirichlet(2, alpha2_100[sample(N_post_alpha,1),]))
}

post_p_sample = rbind(data.frame(p = c(post_p_1), tau2 = 0.044, drug=rep(names,each=1000)),
                     data.frame(p = c(post_p_2), tau2 = 0.011, drug=rep(names,each=1000)),
                     data.frame(p = c(post_p_3), tau2 = 0.005, drug=rep(names,each=1000)))
post_p_sample$tau2 = paste0("tau2=",post_p_sample$tau2)
  
p3 = ggplot(post_p_sample,aes(x=p,fill=drug)) + 
  geom_histogram(color="darkblue",bins=50) +
  facet_grid( drug ~ tau2, scales = "free") + 
  ggtitle("Posterior Predictive Histogram (Normal Prior)") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")

pdf("opioid_normal_predictive_new.pdf", height = 9, width = 9)
print(p3)
dev.off()

post_p2_sample = rbind(data.frame(p = c(post_p2_50), lambda = 50, drug=rep(names,each=1000)),
                      data.frame(p = c(post_p2_75), lambda = 75, drug=rep(names,each=1000)),
                      data.frame(p = c(post_p2_100), lambda = 100, drug=rep(names,each=1000)))
post_p2_sample$lambda = as.factor(paste0("lambda=",post_p2_sample$lambda))
post_p2_sample$lambda = factor(post_p2_sample$lambda, levels=c("lambda=50","lambda=75","lambda=100"))


p4 = ggplot(post_p2_sample,aes(x=p,fill=drug)) + 
  geom_histogram(color="darkblue",bins=50) +
  facet_grid(drug~lambda, scales = "free") + 
  ggtitle("Posterior Predictive Histogram (Exponential Prior)") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")

pdf("opioid_exponential_predictive.pdf", height = 9, width = 9)
print(p4)
dev.off()

test = NULL
for(i in 1:1000)
  test = rbind(test, rdirichlet(1,alpha=rexp(6,rate=100)))