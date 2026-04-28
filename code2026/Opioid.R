rm(list = ls())
source("sampling.R")
library(ggplot2)
library(reshape2)


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


# wrapper
opioid_alpha =function(data, e=1/6, v=1, N_draw = 1000, N_truncations=500)
{
  # v: prior variance, i.e. 

  M = dim(data)[1]
  K = dim(data)[2]
  beta = e/v
  tau = e^2/v
  
  alpha_sample = Gibbs_MD(data, tau, beta, niter=2*N_draw)
  return(alpha_sample[1:N_draw + N_draw, ])
}



# sampling alpha ---------
alpha_1 = opioid_alpha(data,v=1)
alpha_5 = opioid_alpha(data,v=5)
alpha_50 = opioid_alpha(data,v=50)

names = colnames(data)
alpha_sample = rbind(data.frame(alpha = c(alpha_1), v = 1, drug=rep(names,each=1000)),
                     data.frame(alpha = c(alpha_5), v = 5, drug=rep(names,each=1000)),
                     data.frame(alpha = c(alpha_50), v = 50, drug=rep(names,each=1000)))
alpha_sample$v = paste0("prior_var=",alpha_sample$v)

p = ggplot(alpha_sample,aes(x=alpha,fill=drug)) + 
  geom_histogram(color="darkblue",bins=50) +
  facet_grid( drug ~ v, scales = "free") + 
  ggtitle("Posterior Alpha Histogram") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")

pdf("opioid_var.pdf", height = 9, width = 9)
print(p)
dev.off()


alpha_e1 = opioid_alpha(data,e=1)
alpha_e2 = opioid_alpha(data,e=2)
alpha_e3 = opioid_alpha(data,e=3)
alpha_e5 = opioid_alpha(data,e=5)
alpha_e10 = opioid_alpha(data,e=10)
alpha_e50 = opioid_alpha(data,e=50)

alpha_sample2 = rbind(data.frame(alpha = c(alpha_e1), e=1, drug=rep(names,each=1000)),
                     data.frame(alpha = c(alpha_e5), e=5, drug=rep(names,each=1000)),
                     data.frame(alpha = c(alpha_e50), e=50, drug=rep(names,each=1000)))
alpha_sample2$e = paste0("prior_mean=",alpha_sample2$e)

p2 = ggplot(alpha_sample2,aes(x=alpha,fill=drug)) + 
  geom_histogram(color="darkblue",,bins=50) +
  facet_grid( drug ~ e, scales = "free") + 
  ggtitle("Posterior Alpha Histogram") + 
  theme(plot.title = element_text(hjust = 0.5),,legend.position = "none")

pdf("opioid_mean.pdf", height = 9, width = 9)
print(p2)
dev.off()

# posterior predictive ------------------------
post_p_1 = NULL
post_p_5 = NULL
post_p_50 = NULL


N_post_alpha = dim(alpha_e1)[1]
for(i in 1:1000)
{
  post_p_1 = rbind(post_p_1, rdirichlet(1, alpha_e1[sample(N_post_alpha,1),]))
  post_p_5 = rbind(post_p_5, rdirichlet(1, alpha_e5[sample(N_post_alpha,1),]))
  post_p_50 = rbind(post_p_50, rdirichlet(1, alpha_e50[sample(N_post_alpha,1),]))
}

post_p_sample = rbind(data.frame(p = c(post_p_1), e = 1, drug=rep(names,each=1000)),
                     data.frame(p = c(post_p_5), e = 5, drug=rep(names,each=1000)),
                     data.frame(p = c(post_p_50), e = 50, drug=rep(names,each=1000)))
post_p_sample$e = paste0("prior_mean=",post_p_sample$e)
  
p3 = ggplot(post_p_sample,aes(x=p,fill=drug)) + 
  geom_histogram(color="darkblue",bins=50) +
  facet_grid( drug ~ e, scales = "free") + 
  ggtitle("Posterior Predictive Histogram") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")

pdf("opioid_predictive_new.pdf", height = 9, width = 9)
print(p3)
dev.off()

